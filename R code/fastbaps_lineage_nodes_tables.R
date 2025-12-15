#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(fastbaps)
})

# ---------------------------
# Minimal CLI parsing (no optparse dependency)
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}

snp_fasta   <- get_arg("--snp")
tree_file   <- get_arg("--tree")
out_prefix  <- get_arg("--out_prefix", "fastbaps_lineage")
outgroup_f  <- get_arg("--outgroup_file", NA)

if (is.na(snp_fasta) || is.na(tree_file)) {
  cat(
    "Usage:\n",
    "  Rscript fastbaps_lineage_nodes.R --snp core_SNPs.fasta --tree RAxML_tree.newick \\\n",
    "    [--outgroup_file outgroup.txt] [--out_prefix results/lineages]\n\n",
    "Outputs:\n",
    "  <out_prefix>_tips.tsv\n",
    "  <out_prefix>_nodes.tsv\n",
    sep = ""
  )
  quit(status = 2)
}

trim_vec <- function(x) trimws(gsub("^>\\s*", "", x))

# ---------------------------
# 1) Read tree
# ---------------------------
tree <- read.tree(tree_file)
tree$tip.label <- trim_vec(tree$tip.label)

# ---------------------------
# 2) Optionally drop outgroup from tree (and later from lineage tables)
# ---------------------------
outgroup <- character(0)
if (!is.na(outgroup_f)) {
  outgroup <- trim_vec(readLines(outgroup_f, warn = FALSE))
  outgroup <- outgroup[nzchar(outgroup)]
  drop_og <- intersect(tree$tip.label, outgroup)
  if (length(drop_og) > 0) {
    tree <- drop.tip(tree, drop_og)
  }
}

# ---------------------------
# 3) Import SNP FASTA (fastbaps sparse format)
#    NOTE: import_fasta_sparse_nt expects sequences as COLUMNS (samples)
# ---------------------------
snp_data <- import_fasta_sparse_nt(snp_fasta)
colnames(snp_data$snp.matrix) <- trim_vec(colnames(snp_data$snp.matrix))

# Keep only samples present in tree (after dropping outgroup)
common_ids <- intersect(tree$tip.label, colnames(snp_data$snp.matrix))
if (length(common_ids) < 3) {
  stop("Too few overlapping IDs between tree tips and SNP alignment after filtering.")
}

# Subset SNP matrix to common IDs and keep same order as tree tips
common_ids <- tree$tip.label[tree$tip.label %in% common_ids]
snp_data_sub <- snp_data
snp_data_sub$snp.matrix <- snp_data$snp.matrix[, common_ids, drop = FALSE]

# ---------------------------
# 4) Run fastbaps
# ---------------------------
set.seed(42)
baps_hc  <- fast_baps(snp_data_sub)
clusters <- best_baps_partition(snp_data_sub, baps_hc)

if (is.null(names(clusters))) names(clusters) <- colnames(snp_data_sub$snp.matrix)

tip_lineage <- data.frame(
  tip     = names(clusters),
  lineage = paste0("L", as.integer(clusters)),
  stringsAsFactors = FALSE
)

# Ensure outgroup excluded (extra safety in case outgroup wasn't in tree but was in SNPs)
if (length(outgroup) > 0) {
  tip_lineage <- tip_lineage[!(tip_lineage$tip %in% outgroup), , drop = FALSE]
}

# Write tip lineage TSV
tips_out <- paste0(out_prefix, "_tips.tsv")
write.table(tip_lineage[order(tip_lineage$lineage, tip_lineage$tip), ],
            file = tips_out, sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------
# 5) Node lineage attribution
#    For each node, look at descendant tips:
#      - if all descendants share the same lineage => assign that lineage
#      - else => MIXED
# ---------------------------
# Fast lookup: tip -> lineage
lin_map <- setNames(tip_lineage$lineage, tip_lineage$tip)

# Ensure tree tips all have lineage labels (they should, by construction)
missing_lin <- setdiff(tree$tip.label, names(lin_map))
if (length(missing_lin) > 0) {
  stop("Some tree tips have no lineage (check naming consistency):\n  ",
       paste(head(missing_lin, 20), collapse = ", "))
}

n_tips   <- length(tree$tip.label)
n_nodes  <- tree$Nnode
all_nodes <- 1:(n_tips + n_nodes)

# helper: descendant tip labels of a node
desc_tips <- function(node) {
  # node can be a tip itself
  if (node <= n_tips) return(tree$tip.label[node])
  tree$tip.label[phangorn::Descendants(tree, node, type = "tips")[[1]]]
}

# We only need phangorn::Descendants; load phangorn only now
suppressPackageStartupMessages(library(phangorn))

node_df <- lapply(all_nodes, function(nd) {
  tips <- desc_tips(nd)
  lins <- unname(lin_map[tips])
  uniq <- unique(lins)
  
  assigned <- if (length(uniq) == 1) uniq else "MIXED"
  
  data.frame(
    node_id          = nd,
    is_tip           = nd <= n_tips,
    tip_label        = if (nd <= n_tips) tree$tip.label[nd] else NA_character_,
    bootstrap_label  = if (nd > n_tips) {
      # internal node labels in ape are in node.label indexed 1..Nnode
      idx <- nd - n_tips
      if (!is.null(tree$node.label) && idx <= length(tree$node.label)) tree$node.label[idx] else NA_character_
    } else NA_character_,
    assigned_lineage = assigned,
    n_desc_tips      = length(tips),
    stringsAsFactors = FALSE
  )
})

node_df <- do.call(rbind, node_df)

# Write node TSV
nodes_out <- paste0(out_prefix, "_nodes.tsv")
write.table(node_df, file = nodes_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
cat("Wrote:\n  ", tips_out, "\n  ", nodes_out, "\n", sep = "")
 

# SNPs (sites) and samples used by fastBAPS
n_snps    <- nrow(snp_data_sub$snp.matrix)
n_samples <- ncol(snp_data_sub$snp.matrix)

cat("fastBAPS input:", n_snps, "SNP sites across", n_samples, "samples\n")