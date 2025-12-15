############################################
## PACKAGES
############################################
library(tidyverse)
library(lubridate)
library(stringr)
library(pheatmap)
library(svglite)

############################################
## 1. LOAD SNP MATRIX (CLEAN IDs)
############################################
snp_df <- read_tsv(
  "~/Desktop/core_alignment.matrix.clean_names.tsv",
  col_types = cols()
)

# First col = rownames
snp_tmp <- snp_df %>%
  column_to_rownames(var = colnames(.)[1])

rn <- rownames(snp_tmp)
cn <- colnames(snp_tmp)

# Force ALL entries to numeric
snp_mat <- apply(as.matrix(snp_tmp), 2, as.numeric)
snp_mat <- as.matrix(snp_mat)
rownames(snp_mat) <- rn
colnames(snp_mat) <- cn

cat("snp_mat:\n")
str(snp_mat)

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
    
    coll_date_parsed = suppressWarnings(
      parse_date_time(coll_date, orders = c("d/m/Y", "Y"))
    ),
    coll_year = year(coll_date_parsed),
    
    country_raw = case_when(
      is.na(geo) | geo == "" ~ NA_character_,
      TRUE ~ str_trim(str_split_fixed(geo, "[:;,]", 2)[, 1])
    ),
    
    country = case_when(
      country_raw %in% c("missing","not applicable","not provided",
                         "N/A","NA","not collected","unknown",
                         "not determined") ~ NA_character_,
      country_raw == "Viet Nam" ~ "Vietnam",
      TRUE ~ country_raw
    ),
    
    region = case_when(
      country %in% c("China","India","Bangladesh","Pakistan","Vietnam",
                     "Hong Kong","Singapore","Thailand","Japan",
                     "South Korea","Taiwan") ~ "Asia",
      country %in% c("Spain","France","Portugal","Czech Republic",
                     "United Kingdom","UK","Netherlands","The Netherlands",
                     "Switzerland","Romania") ~ "Europe",
      country %in% c("USA","United States","Canada","Haiti",
                     "Argentina","Chile") ~ "Americas",
      country %in% c("Kenya","Nigeria","Senegal","Ghana","Ethiopia",
                     "South Africa","Mozambique","Madagascar",
                     "Gambia","Tanzania") ~ "Africa",
      country %in% c("Australia") ~ "Australia",
      is.na(country) ~ NA_character_,
      TRUE ~ "Other"
    )
  )

############################################
## 3. JOIN SNP DISTANCES + METADATA
############################################
plot_df <- dist_df %>%
  left_join(meta_clean, by = "accession")

############################################
## 4. SELECT H64788 + TOP 10 CLOSEST ISOLATES
############################################
top10 <- plot_df %>%
  arrange(snp_dist) %>%
  slice(1:10)

iso_ids <- c("H64788", top10$accession)

############################################
## 5. SUBSET SNP MATRIX TO THESE 11 ISOLATES
############################################
sub_mat <- snp_mat[iso_ids, iso_ids]

# Build a FRESH numeric matrix for pheatmap
sub_mat_num <- matrix(
  as.numeric(sub_mat),
  nrow = nrow(sub_mat),
  ncol = ncol(sub_mat),
  dimnames = dimnames(sub_mat)
)

cat("sub_mat_num:\n")
str(sub_mat_num)

############################################
## 6. BUILD ANNOTATION DF (YEAR, HOST, REGION)
############################################
ann_df <- meta_clean %>%
  filter(accession %in% iso_ids) %>%
  select(accession, coll_year, host, region) %>%
  mutate(
    accession = factor(accession, levels = iso_ids),
    coll_year = as.factor(coll_year),
    host      = as.factor(host),
    region    = as.factor(region)
  ) %>%
  arrange(accession)

# Drop tibble class, go to plain data.frame
ann_df <- as.data.frame(ann_df)

# Set rownames & drop accession col
rownames(ann_df) <- ann_df$accession
ann_df$accession <- NULL

str(ann_df)


############################################
## 7. ALL-vs-ALL SNP DISTANCE HEATMAP
############################################
pheatmap(
  sub_mat_num,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = TRUE,
  number_format   = "%.0f",
  fontsize_number = 14,
  annotation_row  = ann_df,
  color           = colorRampPalette(c("white", "red"))(100),
  border_color    = "grey80",
  main            = "Core SNP distances among H64788 and its 10 closest isolates"
)

############################################
## 8. SAVE TO SVG
############################################
svglite(
  "~/Desktop/H64788_top10_SNP_heatmap.svg",
  width  = 12,
  height = 12
)

pheatmap(
  sub_mat_num,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = TRUE,
  number_format   = "%.0f",
  fontsize_number = 14,
  annotation_row  = ann_df,
  color           = colorRampPalette(c("white", "red"))(100),
  border_color    = "grey80",
  main            = "Core SNP distances among H64788 and its 10 closest isolates"
)

dev.off()



