############################################
## PACKAGES
############################################
library(tidyverse)
library(lubridate)
library(stringr)
library(ggplot2)


############################################
## 1. LOAD CLEAN SNP MATRIX
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
## 2. LOAD + CLEAN METADATA
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
    coll_date = `Assembly BioSample Collection date`,
    host      = `Assembly BioSample Host`
  ) %>%
  mutate(
    accession = as.character(accession),
    
    ## Date handling
    coll_date_parsed = suppressWarnings(
      parse_date_time(coll_date, orders = c("d/m/Y", "Y"))
    ),
    coll_year = year(coll_date_parsed),
    
    ## Extract country (first part before colon, comma, etc.)
    country_raw = case_when(
      is.na(geo) | geo == "" ~ NA_character_,
      TRUE ~ str_trim(str_split_fixed(geo, "[:;,]", 2)[, 1])
    ),
    
    ## Clean country names
    country = case_when(
      country_raw %in% c("missing","not applicable","not provided",
                         "N/A","NA","not collected","unknown",
                         "not determined") ~ NA_character_,
      country_raw == "Viet Nam" ~ "Vietnam",
      TRUE ~ country_raw
    ),
    
    ## Assign region
    region = case_when(
      ## Asia
      country %in% c("China","India","Bangladesh","Pakistan","Vietnam",
                     "Hong Kong","Singapore","Thailand","Japan",
                     "South Korea","Taiwan") ~ "Asia",
      
      ## Europe
      country %in% c("Spain","France","Portugal","Czech Republic",
                     "United Kingdom","UK","Netherlands",
                     "The Netherlands","Switzerland","Romania") ~ "Europe",
      
      ## Americas
      country %in% c("USA","United States","Canada","Haiti",
                     "Argentina","Chile") ~ "Americas",
      
      ## Africa
      country %in% c("Kenya","Nigeria","Senegal","Ghana","Ethiopia",
                     "South Africa","Mozambique","Madagascar",
                     "Gambia","Tanzania") ~ "Africa",
      
      ## Australia
      country %in% c("Australia") ~ "Australia",
      
      ## Missing
      is.na(country) ~ NA_character_,
      
      ## Everything else stays here
      TRUE ~ "Other"
    )
  )


############################################
## 3. JOIN SNP DISTANCES + METADATA
############################################
plot_df <- dist_df %>%
  left_join(meta_clean, by = "accession")


############################################
## 4. COUNT ISOLATES PER REGION
############################################
region_counts <- plot_df %>%
  filter(!is.na(region)) %>%
  count(region, name = "n") %>%
  dplyr::arrange(desc(n))

print(region_counts)


############################################
## 5. DEFINE CONSISTENT REGION COLOURS
############################################
region_palette <- c(
  "Asia"      = "#1f78b4",
  "Europe"    = "#33a02c",
  "Americas"  = "#e31a1c",
  "Africa"    = "#ff7f00",
  "Australia" = "#6a3d9a",
  "Other"     = "#b15928"   # will not be used now, but fine to keep
)

p_region_bar <- ggplot(region_counts,
                       aes(x = reorder(region, -n), y = n, fill = region)) +
  geom_col(color = "black", width = 0.75) +
  
  # Big numbers on top of bars
  geom_text(
    aes(label = n),
    vjust = -0.7,
    size  = 5,          # ~ 12pt equivalent
    fontface = "bold"
  ) +
  
  scale_fill_manual(values = region_palette, name = "Region") +
  labs(
    x = "Region",
    y = "Number of isolates",
    title = "Total number of isolates per region"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 20, 10)
  ) +
  coord_cartesian(clip = "off")   # so labels aren’t clipped

p_region_bar


library(svglite)
svglite::svglite(
  "~/Desktop/region_bar.svg",
  width = 12,
  height = 12
)
print(p_region_bar)
dev.off()







