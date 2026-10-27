# Genuinely held-out CV: fit OUTRIDER (autoencoder+NB) on each fold's TRAIN
# HC only, then score the held-out TEST HC through the frozen model
# (frozen_outrider.R::score_frozen). Test samples never touch E/D/b/theta/mu
# or the size-factor/centering references -- unlike ../insample_comparison/
# run_outrider_cv.R, which fits train+test jointly and only reads out the
# test rows afterward.

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source("/project/cfRNA_NormativeModeling/MixedEffectsModeling/notebookscripts/OutriderComparison/held_out_comparison/frozen_outrider.R")

IN_DIR <- commandArgs(trailingOnly = TRUE)[1]   # config.OUTRIDER_COMPARISON_DIR
OUT_DIR <- commandArgs(trailingOnly = TRUE)[2]  # config.OUTRIDER_HELD_OUT_DIR

hc_counts <- fread(file.path(IN_DIR, "hc_counts.csv"))
hc_genes <- hc_counts[[1]]
hc_mat <- as.matrix(hc_counts[, -1]); rownames(hc_mat) <- hc_genes

gene_len <- fread(file.path(IN_DIR, "gene_lengths.csv"))
setnames(gene_len, c("gene", "length"))
gene_len <- gene_len[match(hc_genes, gene_len$gene)]
stopifnot(identical(gene_len$gene, hc_genes))

folds <- fromJSON(file.path(IN_DIR, "cv_folds.json"))
bp <- MulticoreParam(workers = 4, stop.on.error = FALSE)

for (fi in names(folds)) {
  cat(sprintf("=== fold %s (held-out) ===\n", fi))
  train_names <- folds[[fi]]$train
  test_names <- folds[[fi]]$test

  tf <- fit_train(hc_mat[, train_names, drop = FALSE], gene_len$length, q = 20, iterations = 5, bp = bp)
  cat(sprintf("  train fit: %d genes kept, %d dropped\n", length(tf$genes), tf$n_dropped))
  tf$ods_train <- NULL  # not needed past validation, keeps saved object small

  sc <- score_frozen(tf, hc_mat[, test_names, drop = FALSE])

  saveRDS(list(z = sc$z, mu = sc$mu_nb, y = sc$counts, theta = tf$theta, genes = tf$genes),
          file.path(OUT_DIR, paste0("cv_fold", fi, "_full.rds")))
}
cat("ALL DONE\n")
