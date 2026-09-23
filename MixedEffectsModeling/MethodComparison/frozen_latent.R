# Delta over held_out_comparison/frozen_outrider.R: lets the frozen fit use either
# OUTRIDER's autoencoder or its PCA confounder control, and persists the whole frozen fit
# so perturbed counts can be re-scored as a matmul instead of a 600s refit.
#
# Everything else -- compute_x, compute_size_factors, score_frozen, retry_drop -- is
# sourced from the original, which validate_frozen.R already proved reproduces native
# OUTRIDER output exactly (max |Z diff| = 0.000e+00).
#
# controlForConfounders(implementation = "pca") dispatches to autoCorrectPCA, which
# accepts only `trim` -- passing `iterations` (a fitAutoencoder argument) errors with
# "unused argument". Hence the conditional below. autoCorrectPCA writes E/D/b into the
# same slots the autoencoder path uses, so score_frozen works unchanged on a PCA fit.

source("/project/cfRNA_NormativeModeling/MixedEffectsModeling/OutriderComparison/held_out_comparison/frozen_outrider.R")

fit_train_impl <- function(train_mat, gene_len, q = 20, implementation = "autoencoder",
                           iterations = 5, bp = SerialParam()) {
  ods <- OutriderDataSet(countData = train_mat)
  mcols(ods)$basepairs <- gene_len
  ods <- filterExpression(ods, filterGenes = TRUE, fpkmCutoff = 1)
  ods <- estimateSizeFactors(ods)
  set.seed(42)

  ctrl_args <- list(q = q, implementation = implementation, BPPARAM = bp)
  if (implementation == "autoencoder") ctrl_args$iterations <- iterations
  r1 <- do.call(retry_drop, c(list(controlForConfounders, ods), ctrl_args))
  r2 <- retry_drop(fit, r1$ods)
  ods <- r2$ods

  xc <- compute_x(counts(ods, normalized = FALSE), sizeFactors(ods))
  l2fc_train <- .log2fc(ods)

  list(
    implementation = implementation,
    genes = rownames(ods),
    E = metadata(ods)[["E"]],
    D = metadata(ods)[["D"]],
    b = mcols(ods)[["b"]],
    gene_mu = mcols(ods)[["mu"]],
    theta = theta(ods),
    loggeomeans = mcols(ods)[["loggeomeans"]],
    x_center = xc$center,
    l2fc_mean = rowMeans(l2fc_train),
    l2fc_sd = matrixStats::rowSds(as.matrix(l2fc_train)),
    n_dropped = r1$n_dropped + r2$n_dropped
  )
}
