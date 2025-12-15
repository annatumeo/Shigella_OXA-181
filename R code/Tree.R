############################################################
## Shigella: circular bootstrap tree + rings
## Rings: Collection year, Region, OXA_genes (which)
## Uses metadata with lineage already merged
############################################################

suppressPackageStartupMessages({
  library(treeio)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggnewscale)
  library(phangorn)
})

setwd("~/Desktop")

tree_file <- "RAxML_bipartitionsBranchLabels.all200.renamed.newick"
meta_file <- "metadata_matched_to_tree_with_lineage.tsv"

out_png <- "Shigella_tree_bootstrap_rings_OXAgenes_year_region.png"
out_pdf <- "Shigella_tree_bootstrap_rings_OXAgenes_year_region.pdf"
out_svg <- "Shigella_tree_bootstrap_rings_OXAgenes_year_region.svg"

FOCAL_IDS    <- c("H64788")
hi_ids       <- unique(c(OUTGROUP_IDS, FOCAL_IDS))

norm_id <- function(x) str_replace_all(as.character(x), "\\s+", "")

############################################################
## 1) Load metadata
############################################################
meta <- read.delim(meta_file, sep = "\t", header = TRUE,
                   stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(all(c("Tree_label","Collection Year","Region","OXA_genes","lineage") %in% colnames(meta)))

meta <- meta %>%
  mutate(
    Tree_label = norm_id(Tree_label),
    YEAR       = suppressWarnings(as.integer(`Collection Year`)),
    Region     = ifelse(is.na(Region) | Region == "", "Unknown", Region),
    lineage    = ifelse(is.na(lineage) | lineage == "", "Unassigned", lineage),
    OXA_genes  = ifelse(is.na(OXA_genes), "", as.character(OXA_genes))
  )

# Parse OXA_genes into a single "which" category per isolate
# - If multiple genes, keep them joined with "+"
# - If blank, set to "None"
split_oxa <- function(x) {
  x <- str_trim(x)
  if (is.na(x) || x == "") return("None")
  parts <- unlist(str_split(x, "[,;\\s]+"))
  parts <- parts[parts != ""]
  parts <- sort(unique(parts))
  if (length(parts) == 0) "None" else paste(parts, collapse = "+")
}

meta$OXA_which <- vapply(meta$OXA_genes, split_oxa, character(1))

# Keep OXA levels stable (None first, then alphabetical)
oxa_levels <- c("None", sort(unique(meta$OXA_which[meta$OXA_which != "None"])))
meta$OXA_which <- factor(meta$OXA_which, levels = unique(oxa_levels))

############################################################
## 2) Read tree + midpoint root
############################################################
raxml_td <- read.raxml(tree_file)
raxml_td@phylo <- phangorn::midpoint(raxml_td@phylo)
raxml_td@phylo$tip.label <- norm_id(raxml_td@phylo$tip.label)

tips <- raxml_td@phylo$tip.label

############################################################
## 3) Join tips to metadata
############################################################
tip_meta <- tibble(label = tips) %>%
  left_join(meta, by = c("label" = "Tree_label")) %>%
  mutate(
    lineage   = ifelse(is.na(lineage) | lineage == "", "Unassigned", lineage),
    Region    = ifelse(is.na(Region)  | Region  == "", "Unknown", Region),
    OXA_which = ifelse(is.na(OXA_which), "None", as.character(OXA_which)),
    OXA_which = factor(OXA_which, levels = levels(meta$OXA_which))
  )

############################################################
## 4) Base tree plot
############################################################
p <- ggtree(
  raxml_td,
  layout        = "circular",
  branch.length = "none",
  aes(color = bootstrap)
) +
  scale_color_viridis_c(name = "Bootstrap", na.value = "grey70") +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 55, 5, 5, unit = "mm")
  )

# tip points by lineage
p <- p + ggnewscale::new_scale_color()
p <- p %<+% tip_meta +
  geom_tippoint(aes(color = lineage), size = 1.6, alpha = 0.95) +
  guides(color = guide_legend(title = "Lineage", override.aes = list(size = 3)))
# Highlight + print ONLY these two labels (strings must be quoted!)
OUTGROUP_IDS <- c("GCF_002950395_1")
FOCAL_IDS    <- c("H64788")
hi_ids       <- c(OUTGROUP_IDS, FOCAL_IDS)


############################################################
## 5) Rings via geom_tile at radial x offsets (stable)
############################################################
pd <- p$data
tip_pos <- pd %>% filter(isTip) %>% select(label, y, x)

base_x <- max(pd$x, na.rm = TRUE)

ring_thick <- 4.0
gap        <- 1.8

# (optional) move rings slightly further out so they don't crowd the tree
x_year   <- base_x + 1.5   # was +0.7
x_oxa    <- x_year + ring_thick + gap
x_region <- x_oxa  + ring_thick + gap


# YEAR ring
dat_year <- tip_meta %>%
  select(label, YEAR) %>%
  left_join(tip_pos %>% select(label, y), by = "label") %>%
  mutate(x = x_year) %>%
  filter(!is.na(y))

# OXA WHICH ring
dat_oxa <- tip_meta %>%
  select(label, OXA_which) %>%
  left_join(tip_pos %>% select(label, y), by = "label") %>%
  mutate(x = x_oxa) %>%
  filter(!is.na(y))

# REGION ring
region_levels <- c("Africa","Americas","Asia","Australia","Europe","Unknown")
dat_region <- tip_meta %>%
  transmute(label,
            Region = factor(ifelse(is.na(Region) | Region == "", "Unknown", Region),
                            levels = region_levels)) %>%
  left_join(tip_pos %>% select(label, y), by = "label") %>%
  mutate(x = x_region) %>%
  filter(!is.na(y))

# expand x-limits so rings are not clipped
x_max_needed <- x_region + ring_thick + 3.0
p <- p + expand_limits(x = x_max_needed)
p <- p +
  theme(plot.margin = margin(5, 120, 5, 5, unit = "mm"))
# Add YEAR ring
p <- p +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = dat_year,
    aes(x = x, y = y, fill = YEAR),
    inherit.aes = FALSE,
    width = ring_thick,
    height = 1
  ) +
scale_fill_viridis_c(
  name = "Collection year",
  option = "magma",   # try: "plasma", "inferno", "cividis", "magma", "turbo"
  direction = 1,
  na.value = "grey95"
)
# Add OXA WHICH ring (categorical)
p <- p +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = dat_oxa,
    aes(x = x, y = y, fill = OXA_which),
    inherit.aes = FALSE,
    width = ring_thick,
    height = 1
  ) +
  scale_fill_discrete(name = "OXA genes (which)")

# Add REGION ring
p <- p +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = dat_region,
    aes(x = x, y = y, fill = Region),
    inherit.aes = FALSE,
    width = ring_thick,
    height = 1
  ) +
  scale_fill_manual(
    name = "Region",
    values = c(
      Africa    = "#1b9e77",
      Americas  = "#d95f02",
      Asia      = "#7570b3",
      Europe    = "#e7298a",
      Australia = "#66a61e",
      Unknown   = "grey70"
    )
  ) +
  theme(
    legend.box = "vertical",
    legend.spacing.y = unit(0.15, "cm"),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )
# ---- FORCE highlights on top (add AFTER rings) ----
p <- p +
  geom_tippoint(
    aes(subset = isTip & label %in% FOCAL_IDS),
    shape = 21, fill = "yellow", color = "black", stroke = 1.1, size = 5.0
  ) +
  geom_tiplab(aes(subset = isTip & label %in% FOCAL_IDS, label = label),
    size = 3.4,
    fontface = "bold",
    align = TRUE,
    linetype = 1,
    linesize = 0.25,
    offset = 11,
    show.legend = FALSE
  )


print(p)

############################################################
## 6) Save
ggsave(out_png, p, width = 11, height = 11, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(out_pdf, p, width = 11, height = 11, units = "in", bg = "white", limitsize = FALSE)
ggsave(out_svg, p, width = 11, height = 11, units = "in", bg = "white", limitsize = FALSE)


cat("Saved:\n", out_png, "\n", out_pdf, "\n", out_svg, "\n")

