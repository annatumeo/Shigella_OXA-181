############################################
## Bubble chart: ARG prevalence per lineage
## Facet by antibiotic class
## Removes Unassigned
##
## Inputs (in working dir):
##  - metadata_matched_to_tree_with_lineage.tsv
##  - abricate.arg.summary.cleaned.tsv
##  - non_core_ARGs_with_classes.tsv
##
## Outputs:
##  - bubble_ARG_prevalence_by_lineage.tsv  (plot-ready data)
##  - bubble_ARG_prevalence_by_lineage.svg
############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- Files ----
meta_file   <- "metadata_matched_to_tree_with_lineage.tsv"
arg_file    <- "abricate.arg.summary.cleaned.tsv"
class_file  <- "non_core_ARGs_with_classes.tsv"

# ---- Options ----
DROP_CLASSES <- c("multidrug efflux/regulator")  # optional
MIN_LINEAGE_N <- 10   # drop tiny lineages if you want stability (set 1 to keep all)
MIN_PRESENT   <- 5    # keep genes present in >= this many genomes overall
MIN_ABSENT    <- 5    # and absent in >= this many genomes overall

# ---- Helper: normalize IDs ----
norm_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\\s+", "")
  # optional: convert trailing .1 -> _1 (matches many of your cleaned IDs)
  x <- stringr::str_replace(x, "\\.(\\d+)$", "_\\1")
  x
}

# ---- 1) Metadata: sample + lineage ----
meta <- readr::read_tsv(meta_file, show_col_types = FALSE) %>%
  transmute(
    sample  = norm_id(Tree_label),
    lineage = if_else(is.na(lineage) | lineage == "", "Unassigned", lineage)
  ) %>%
  filter(lineage != "Unassigned") %>%
  add_count(lineage, name = "lineage_n") %>%
  filter(lineage_n >= MIN_LINEAGE_N) %>%
  select(-lineage_n)

# ---- 2) ARG table (ABRicate summary) ----
arg <- readr::read_tsv(arg_file, show_col_types = FALSE) %>%
  mutate(sample = norm_id(`#FILE`)) %>%
  select(-`#FILE`)

# Merge to keep only samples with lineage
dat <- arg %>%
  inner_join(meta, by = "sample")

gene_cols <- setdiff(colnames(dat), c("sample", "lineage"))

# Convert to presence/absence:
# ABRicate summary convention: "." (or blank/NA) = absent; else present
bin <- dat %>%
  mutate(across(all_of(gene_cols), ~ !is.na(.) & str_trim(as.character(.)) != "." & str_trim(as.character(.)) != ""))
# ---- 3) Gene -> antibiotic class mapping ----
gene_class <- readr::read_tsv(class_file, show_col_types = FALSE) %>%
  transmute(
    gene = as.character(gene),
    antibiotic_class_raw = as.character(antibiotic_class)
  ) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(
    antibiotic_class = str_to_lower(str_trim(antibiotic_class_raw)),
    antibiotic_class = case_when(
      # ---- STANDARDISE SYNONYMS (PREVENT DUPLICATE FACETS) ----
      antibiotic_class %in% c("biocide (qac)", "qac", "biocide/qac", "biocide") ~ "biocide",
      antibiotic_class %in% c("polymyxin/colistin (lps modification)",
                              "polymyxin/colistin (lps-modification)",
                              "polymyxin/colistin", "polymyxin", "colistin") ~ "polymyxin/colistin",
      antibiotic_class %in% c("quinolone/fluoroquinolone", "fluoroquinolone", "quinolone") ~ "quinolone/fluoroquinolone",
      antibiotic_class %in% c("macrolide/lincosamide/streptogramin",
                              "macrolide", "lincosamide", "streptogramin") ~ "macrolide/lincosamide/streptogramin",
      TRUE ~ antibiotic_class
    )
  ) %>%
  filter(!(antibiotic_class %in% str_to_lower(DROP_CLASSES))) %>%
  mutate(
    # ---- MERGED FACET GROUPS ----
    facet_group = case_when(
      antibiotic_class %in% c("fosfomycin", "polymyxin/colistin") ~ "Polymyxin/Colistin/Fosfomycin",
      antibiotic_class %in% c("rifamycin", "tetracycline", "biocide") ~ "Biocide/Tetracycline/Rifamycin",
      antibiotic_class %in% c("sulfonamide", "macrolide/lincosamide/streptogramin") ~ "Sulfonamide/Macrolide",
      TRUE ~ str_to_title(antibiotic_class)
    )
  )



# Keep only genes that exist in your ABRicate matrix AND have a class label
genes_keep <- intersect(gene_cols, gene_class$gene)

# ---- 4) Optional global prevalence filters (avoid ultra-rare genes) ----
gene_prev <- tibble(gene = genes_keep) %>%
  mutate(
    n_present = map_int(gene, ~ sum(bin[[.x]], na.rm = TRUE)),
    n_total   = nrow(bin),
    n_absent  = n_total - n_present
  ) %>%
  filter(n_present >= MIN_PRESENT, n_absent >= MIN_ABSENT)

genes_keep2 <- gene_prev$gene

# ---- 5) Build plot-ready data: prevalence per lineage per gene ----
bubble_df <- bin %>%
  select(sample, lineage, all_of(genes_keep2)) %>%
  pivot_longer(cols = all_of(genes_keep2), names_to = "gene", values_to = "present") %>%
  group_by(lineage, gene) %>%
  summarise(
    n_present = sum(present, na.rm = TRUE),
    n_total   = n(),
    prevalence = n_present / n_total,
    .groups = "drop"
  ) %>%
  left_join(gene_class, by = "gene") %>%
  filter(!is.na(antibiotic_class)) %>%
  mutate(
    # nicer ordering: genes with higher max prevalence first within each class
    gene = fct_reorder(gene, prevalence, .fun = max),
    lineage = fct_infreq(lineage),
    facet_group = fct_infreq(facet_group),
    antibiotic_class = fct_infreq(antibiotic_class),
    label = paste0(n_present, "/", n_total)
  )

# ---- 6) Write the exact data you will plot ----
readr::write_tsv(bubble_df, "bubble_ARG_prevalence_by_lineage.tsv")
# ---- 7) Bubble chart (facet by antibiotic class) ----
# Make sure lineage is a factor with stable ordering
bubble_df <- bubble_df %>%
  mutate(
    lineage = factor(lineage, levels = sort(unique(as.character(lineage))))
  )

p <- ggplot(bubble_df, aes(x = lineage, y = gene)) +
  geom_point(aes(size = prevalence, color = lineage), alpha = 0.85) +
  facet_wrap(~ facet_group, scales = "free_y", ncol = 3) +
  scale_size_continuous(
    name   = "Prevalence",
    range  = c(0.5, 9),
    labels = scales::percent_format(accuracy = 1)
  ) +
  # one colour per lineage (discrete)
  scale_color_viridis_d(name = "Lineage", end = 0.9) +
  labs(
    x = "Lineage",
    y = "ARG",
    title = "Accessory ARG prevalence by lineage (faceted by antibiotic class)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    strip.text         = element_text(face = "bold")
  )

print(p)
ggsave("bubble_ARG_prevalence_by_lineage.svg", p, width = 13, height = 8, units = "in", bg = "white")
