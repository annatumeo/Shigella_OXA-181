############################################
## Packages
############################################
library(tidyverse)
library(lubridate)
library(stringr)
library(ggridges)


############################################
## 1. Load SNP distance matrix (clean IDs)
############################################
snp_df <- read_tsv(
  "~/Desktop/core_alignment.matrix.clean_names.tsv",
  col_types = cols()
)

snp_mat <- snp_df %>%
  column_to_rownames(var = colnames(.)[1]) %>%
  as.matrix()

ref <- "H64788"

if (!ref %in% rownames(snp_mat)) stop("H64788 not found in matrix")

ref_dists <- snp_mat[ref, ]

dist_df <- tibble(
  accession = names(ref_dists),
  snp_dist  = as.numeric(ref_dists)
) %>%
  filter(accession != ref)

############################################
## 2. Load & clean metadata
############################################
meta <- read_tsv(
  "~/Desktop/accessions_metadata.tsv",
  col_types = cols()
)
meta_clean <- meta %>%
  rename(
    accession = `Assembly Accession`,
    organism  = `ANI Best ANI match Organism`,
    geo       = `Assembly BioSample Geographic location`,
    coll_date = `Assembly BioSample Collection date`
  ) %>%
  mutate(
    accession = as.character(accession),
    coll_date_parsed = suppressWarnings(
      parse_date_time(coll_date, orders = c("d/m/Y", "Y"))
    ),
    coll_year = year(coll_date_parsed),
    country_raw = case_when(
      is.na(geo) | geo == "" ~ NA_character_,
      TRUE ~ str_trim(str_split_fixed(geo, "[:;,]", 2)[, 1])
    ),
    country = case_when(
      country_raw %in% c("missing", "not applicable", "not provided", "N/A", "NA") ~ NA_character_,
      country_raw == "Viet Nam" ~ "Vietnam",
      TRUE ~ country_raw
    ),
    region = case_when(
      # Asia
      country %in% c("China", "India", "Bangladesh", "Pakistan", "Vietnam",
                     "Hong Kong", "Singapore") ~ "Asia",
      # Europe
      country %in% c("Spain", "France", "Portugal", "Czech Republic",
                     "United Kingdom", "UK", "Netherlands", "The Netherlands") ~ "Europe",
      # Americas
      country %in% c("USA", "United States", "Canada", "Haiti") ~ "Americas",
      # Africa
      country %in% c("Kenya", "Nigeria", "Senegal") ~ "Africa",
      # Australia as its own thing
      country %in% c("Australia") ~ "Australia",
      is.na(country) ~ NA_character_,
      TRUE ~ "Other"
    )
  )
############################################
## 2b. Define SAME region colours as barplot
############################################
region_palette <- c(
  "Asia"      = "#1f78b4",
  "Europe"    = "#33a02c",
  "Americas"  = "#e31a1c",
  "Africa"    = "#ff7f00",
  "Australia" = "#6a3d9a",
  "Other"     = "#b15928"
)
############################################
## 3. Join SNP distances + metadata
############################################
plot_df <- dist_df %>%
  left_join(meta_clean, by = "accession")

############################################
## 4. Keep only countries with enough genomes
############################################
min_n <- 5

country_counts <- plot_df %>%
  filter(!is.na(country)) %>%             # drop missing / junk
  count(country, name = "n") %>%
  arrange(desc(n))

top_countries <- country_counts %>%
  filter(n >= min_n) %>%
  pull(country)
ridge_df <- plot_df %>%
  filter(!is.na(country)) %>%
  group_by(country) %>%
  filter(n() >= min_n) %>%
  ungroup() %>%
  mutate(
    country = fct_reorder(country, snp_dist, .fun = median, .desc = TRUE)
  )
p_ridge <- ggplot(
  ridge_df,
  aes(x = snp_dist, y = country, fill = region)
) +
  geom_density_ridges(
    alpha          = 0.8,
    scale          = 3,
    rel_min_height = 0.01,
    color          = "white",
    size           = 0.3
  ) +
  scale_fill_manual(
    values = region_palette,
    na.translate = FALSE    # optional: drop NA from legend
  ) +
  labs(
    x = "SNP distance from H64788",
    y = NULL,
    title    = "Distribution of SNP distances from H64788 by country",
    subtitle = paste0("Countries with ≥", min_n, " genomes"),
    fill     = "Region"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold")
  )

p_ridge

library(svglite)
svglite::svglite(
  "~/Desktop/H64788_SNP_ridgeline_by_country.svg",
  width = 12,
  height = 12
)
print(p_ridge)
dev.off()



