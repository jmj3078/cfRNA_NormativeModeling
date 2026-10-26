# Genuinely held-out LOBO: fit OUTRIDER on HC from every OTHER batch, then
# score the held-out batch's HC+disease samples through the frozen
# encoder/decoder/NB params (frozen_outrider.R::score_frozen). The held-out
# batch's disease samples never enter the autoencoder fit at all -- unlike
# ../insample_comparison/run_outrider_lobo.R, which concatenates them into
# the same matrix passed to controlForConfounders().

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source("/project/cfRNA_NormativeModeling/MixedEffectsModeling/OutriderComparison/held_out_comparison/frozen_outrider.R")

IN_DIR <- commandArgs(trailingOnly = TRUE)[1]   # config.OUTRIDER_COMPARISON_DIR
OUT_DIR <- commandArgs(trailingOnly = TRUE)[2]  # config.OUTRIDER_HELD_OUT_DIR

hc_counts <- fread(file.path(IN_DIR, "hc_counts.csv"))
hc_genes <- hc_counts[[1]]
hc_mat <- as.matrix(hc_counts[, -1]); rownames(hc_mat) <- hc_genes
hc_meta <- fread(file.path(IN_DIR, "hc_meta.csv"))
stopifnot(identical(colnames(hc_mat), hc_meta$sample))

dis_counts <- fread(file.path(IN_DIR, "disease_counts.csv"))
dis_genes <- dis_counts[[1]]
dis_mat <- as.matrix(dis_counts[, -1]); rownames(dis_mat) <- dis_genes
stopifnot(identical(hc_genes, dis_genes))

gene_len <- fread(file.path(IN_DIR, "gene_lengths.csv"))
setnames(gene_len, c("gene", "length"))
gene_len <- gene_len[match(hc_genes, gene_len$gene)]
stopifnot(identical(gene_len$gene, hc_genes))

test_meta <- fromJSON(file.path(IN_DIR, "lobo_test_meta.json"))
batches <- names(test_meta)
bp <- MulticoreParam(workers = 4, stop.on.error = FALSE)

out_rows <- list()
for (b in batches) {
  safe <- gsub(" ", "_", b)
  cat(sprintf("=== %s (held-out) ===\n", b))
  test_names <- test_meta[[b]]$test_names
  test_is_hc <- test_meta[[b]]$test_is_hc

  train_hc_names <- hc_meta$sample[hc_meta$batch != b]  # this batch's HC fully excluded too

  test_from_hc <- intersect(test_names, colnames(hc_mat))
  test_from_dis <- intersect(test_names, colnames(dis_mat))
  test_mat <- cbind(hc_mat[, test_from_hc, drop = FALSE], dis_mat[, test_from_dis, drop = FALSE])
  test_mat <- test_mat[, test_names, drop = FALSE]  # restore original order

  tf <- fit_train(hc_mat[, train_hc_names, drop = FALSE], gene_len$length, q = 20, iterations = 5, bp = bp)
  cat(sprintf("  train fit: %d genes kept, %d dropped\n", length(tf$genes), tf$n_dropped))
  tf$ods_train <- NULL

  sc <- score_frozen(tf, test_mat)
  saveRDS(sc$z, file.path(OUT_DIR, paste0("z_test_", safe, ".rds")))
  fwrite(as.data.table(sc$z, keep.rownames = "sample"), file.path(OUT_DIR, paste0("z_test_", safe, ".csv")))
  out_rows[[b]] <- data.table(batch = b, sample = test_names, is_hc = test_is_hc)
}

fwrite(rbindlist(out_rows), file.path(OUT_DIR, "outrider_lobo_membership.csv"))
cat("ALL DONE\n")
