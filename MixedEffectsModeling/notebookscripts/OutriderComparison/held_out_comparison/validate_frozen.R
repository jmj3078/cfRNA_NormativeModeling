# Self-check: if score_frozen() applied to the TRAIN samples themselves
# reproduces OUTRIDER's own native in-sample Z-scores/p-values, the frozen
# encode/decode math in frozen_outrider.R is a faithful reimplementation --
# only then is applying it to genuinely unseen samples trustworthy.

suppressPackageStartupMessages({
  library(OUTRIDER)
})
source("/project/cfRNA_NormativeModeling/MixedEffectsModeling/notebookscripts/OutriderComparison/held_out_comparison/frozen_outrider.R")

set.seed(1)
ods0 <- makeExampleOutriderDataSet(dataset = "Kremer")
mat <- counts(ods0)
mat <- mat[rowSums(mat) > 0, ]
mat <- mat[seq_len(min(400, nrow(mat))), ]  # small for a fast smoke test
gene_len <- rep(2000L, nrow(mat))

tf <- fit_train(mat, gene_len, q = 10, iterations = 3)
cat(sprintf("train fit: %d genes kept (%d dropped)\n", length(tf$genes), tf$n_dropped))

# native in-sample outputs from the very same fitted ods
ods_native <- computeZscores(tf$ods_train)
ods_native <- computePvalues(ods_native, method = "None")
native_z <- t(as.matrix(assay(ods_native, "zScore")))
native_p <- t(pValue(ods_native))

sc <- score_frozen(tf, mat)  # score the TRAIN samples through the "frozen" path

common_genes <- tf$genes
z_diff <- abs(sc$z[, common_genes] - native_z[, common_genes])
p_diff <- abs(sc$pval[, common_genes] - native_p[, common_genes])

cat(sprintf("max |Z diff| = %.3e, max |p diff| = %.3e\n", max(z_diff, na.rm = TRUE), max(p_diff, na.rm = TRUE)))
stopifnot(max(z_diff, na.rm = TRUE) < 1e-6)
stopifnot(max(p_diff, na.rm = TRUE) < 1e-6)
cat("PASS: frozen scoring reproduces native OUTRIDER output on in-sample data.\n")
