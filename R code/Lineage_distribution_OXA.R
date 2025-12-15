############################################################
## Nature-level Figure A–B
## A) Lineage × Region distribution (within-lineage proportions)
## B) OXA prevalence within lineage × region
##
## - Consistent ordering across panels
## - Legend matches gradient (alpha overridden in legend)
## - Alpha encodes support (n)
## - Large fonts
## - White text on dark tiles
##
## Input:  shigella_metadata_with_lineage_year_region_oxa.tsv
## Output: SVG / PDF / PNG
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
})

#-----------------------------------------------------------
# 0) GLOBAL SETTINGS
#-----------------------------------------------------------
infile <- "shigella_metadata_with_lineage_year_region_oxa.tsv"

BASE_SIZE <- 18      # global font size
SMOOTHING <- "laplace"  # "none" or "laplace"

ALPHA_MIN_N <- 3
ALPHA_MAX_N <- 80

#-----------------------------------------------------------
# 1) LOAD DATA + DEFINE any_OXA
#-----------------------------------------------------------
df <- readr::read_tsv(infile, show_col_types = FALSE)

oxa_bool_cols <- names(df)[stringr::str_detect(names(df), "^OXA-")]

df <- df %>%
  mutate(
    any_OXA = if (length(oxa_bool_cols) > 0) {
      if_any(all_of(oxa_bool_cols), ~ .x %in% TRUE)
    } else {
      !is.na(OXA_genes) & OXA_genes != ""
    }
  ) %>%
  filter(!is.na(lineage),
         !is.na(Region),
         Region != "Unknown")

#-----------------------------------------------------------
# 2) SHARED ORDERING (CRITICAL)
#-----------------------------------------------------------
region_order <- df %>%
  count(Region, sort = TRUE) %>%
  pull(Region)

lineage_order <- df %>%
  count(lineage, sort = TRUE) %>%
  pull(lineage)

df <- df %>%
  mutate(
    Region  = factor(Region, levels = region_order),
    lineage = factor(lineage, levels = lineage_order)
  )

#-----------------------------------------------------------
# 3) PANEL A — LINEAGE × REGION DISTRIBUTION
#-----------------------------------------------------------
tab_lr <- with(df, table(lineage, Region))
ct_lr  <- suppressWarnings(chisq.test(tab_lr, correct = FALSE))

N_lr <- sum(tab_lr)
V_lr <- sqrt(as.numeric(ct_lr$statistic) /
               (N_lr * min(nrow(tab_lr)-1, ncol(tab_lr)-1)))

A_df <- as.data.frame(tab_lr) %>%
  as_tibble() %>%
  rename(n = Freq) %>%
  group_by(lineage) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

subA <- sprintf(
  "χ²(%d)=%.2f, p=%s; Cramér’s V=%.2f; N=%d (Unknown excluded)",
  ct_lr$parameter,
  as.numeric(ct_lr$statistic),
  format.pval(ct_lr$p.value, digits = 3, eps = 1e-300),
  V_lr, N_lr
)

pA <- ggplot(A_df, aes(x = Region, y = lineage, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = if_else(n == 0, "", as.character(n))),
            size = 5,
            color = "black") +
  scale_fill_gradient(
    low = "#f7fbff",
    high = "#2171b5",
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    name = "Within-lineage\nproportion"
  ) +
  labs(
    title = "A  Lineage × Region distribution",
    subtitle = subA,
    x = "Region",
    y = "Lineage"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = BASE_SIZE + 4),
    plot.subtitle = element_text(size = BASE_SIZE),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE),
    axis.text.y = element_text(size = BASE_SIZE),
    axis.title  = element_text(size = BASE_SIZE + 2),
    legend.title = element_text(size = BASE_SIZE + 1),
    legend.text  = element_text(size = BASE_SIZE),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

#-----------------------------------------------------------
# 4) PANEL B — OXA PREVALENCE
#-----------------------------------------------------------
B_df <- df %>%
  group_by(lineage, Region) %>%
  summarise(
    n = n(),
    n_oxa = sum(any_OXA),
    .groups = "drop"
  ) %>%
  mutate(
    prop_raw = n_oxa / n,
    prop = if (SMOOTHING == "laplace")
      (n_oxa + 0.5) / (n + 1) else prop_raw,
    
    label = paste0(n_oxa, "/", n),
    
    n_clip = pmax(pmin(n, ALPHA_MAX_N), ALPHA_MIN_N),
    alpha  = (n_clip - ALPHA_MIN_N) / (ALPHA_MAX_N - ALPHA_MIN_N),
    alpha  = alpha * 0.35 + 0.65,
    
    label_col = if_else(prop >= 0.55, "white", "black")
  )

grp <- df %>% unite(group, lineage, Region, remove = FALSE)
tab_g <- table(grp$group, grp$any_OXA)
ct_g  <- suppressWarnings(chisq.test(tab_g, correct = FALSE))

subB <- sprintf(
  "Fill shows %% OXA+ within each lineage×region (smoothing: %s). Transparency scales with n.  χ²(%d)=%.1f, p=%s",
  SMOOTHING,
  ct_g$parameter,
  as.numeric(ct_g$statistic),
  format.pval(ct_g$p.value, digits = 3, eps = 1e-300)
)

pB <- ggplot(B_df, aes(x = Region, y = lineage)) +
  geom_tile(aes(fill = prop, alpha = alpha),
            color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = label_col),
            size = 5.2) +
  scale_color_identity() +
  scale_fill_gradient(
    low  = "#f7fbff",
    high = "#2171b5",
    limits = c(0, 1),
    labels = percent_format(accuracy = 1),
    name = "% OXA+",
    guide = guide_colorbar(
      barheight = grid::unit(60, "pt"),
      override.aes = list(alpha = 1)
    )
  ) +
  scale_alpha(range = c(0.65, 1), guide = "none") +
  labs(
    title = "B  OXA prevalence within lineage × region strata",
    subtitle = subB,
    x = "Region",
    y = "Lineage"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold", size = BASE_SIZE + 4),
    plot.subtitle = element_text(size = BASE_SIZE),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE),
    axis.text.y = element_text(size = BASE_SIZE),
    axis.title  = element_text(size = BASE_SIZE + 2),
    legend.title = element_text(size = BASE_SIZE + 1),
    legend.text  = element_text(size = BASE_SIZE),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

#-----------------------------------------------------------
# 5) COMBINE + SAVE
#-----------------------------------------------------------
fig <- pA / pB + plot_layout(heights = c(1, 1))

print(fig)

ggsave("Figure_AB_LineageRegion_and_OXA.svg", fig, width = 15, height = 15)
ggsave("Figure_AB_LineageRegion_and_OXA.pdf", fig, width = 15, height = 15)
ggsave("Figure_AB_LineageRegion_and_OXA.png", fig, width = 8.5, height = 11, dpi = 450)
