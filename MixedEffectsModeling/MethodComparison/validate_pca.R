# Same contract validate_frozen.R enforces for the autoencoder, applied to the PCA arm:
# scoring the TRAIN samples through the frozen path must reproduce native OUTRIDER output.
# Also confirms controlForConfounders(implementation="pca") runs at all once `iterations`
# is withheld.
suppressPackageStartupMessages({ library(OUTRIDER) })
source("/project/cfRNA_NormativeModeling/MixedEffectsModeling/_temp_method_comparison/frozen_latent.R")

set.seed(42)
ods0 <- makeExampleOutriderDataSet(dataset = "Kremer")
keep <- head(which(rowSums(counts(ods0)) > 0), 400)
mat <- counts(ods0)[keep, ]
gene_len <- rep(2000, nrow(mat))

for (impl in c("pca", "autoencoder")) {
  cat("===", impl, "===\n")
  tf <- fit_train_impl(mat, gene_len, q = 10, implementation = impl, iterations = 3,
                       bp = SerialParam())
  cat("  genes kept:", length(tf$genes), " dropped:", tf$n_dropped,
      " E:", paste(dim(tf$E), collapse = "x"), "\n")

  # native in-sample reference, refit the same way
  ods <- OutriderDataSet(countData = mat)
  mcols(ods)$basepairs <- gene_len
  ods <- filterExpression(ods, filterGenes = TRUE, fpkmCutoff = 1)
  ods <- estimateSizeFactors(ods)
  set.seed(42)
  args <- list(q = 10, implementation = impl, BPPARAM = SerialParam())
  if (impl == "autoencoder") args$iterations <- 3
  r1 <- do.call(retry_drop, c(list(controlForConfounders, ods), args))
  ods <- retry_drop(fit, r1$ods)$ods
  ods <- computePvalues(computeZscores(ods), method = "None")

  common <- intersect(tf$genes, rownames(ods))
  sc <- score_frozen(tf, mat)
  zd <- max(abs(sc$z[, common] - t(as.matrix(assay(ods, "zScore")))[, common]), na.rm = TRUE)
  pd <- max(abs(sc$pval[, common] - t(pValue(ods))[, common]), na.rm = TRUE)
  cat(sprintf("  max |Z diff| = %.3e, max |p diff| = %.3e\n", zd, pd))
  stopifnot(zd < 1e-6, pd < 1e-6)
  cat("  PASS\n")
}
cat("ALL PASS\n")
