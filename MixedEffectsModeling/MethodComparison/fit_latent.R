# Fit the frozen latent reference (autoencoder and PCA) on each fold's TRAIN HC only, and
# export every parameter needed to score unseen samples later.
#
# run_held_out_cv.R discards these (it keeps only z/mu/y/theta for the unperturbed test
# counts), so perturbed counts cannot be re-scored without a refit. Exporting E/D/b here
# makes every later re-score a matmul in numpy -- normF = exp((x %*% E) %*% t(D) + b) * sf,
# verified against OUTRIDER:::predictMatC to 2e-13 relative.
#
# Each (implementation, fold) is cached, so a rerun only fits what is missing.
#
# Usage: Rscript fit_latent.R [cv|lobo] [fold]      (no fold = all folds of that split)
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
source("/project/cfRNA_NormativeModeling/MixedEffectsModeling/_temp_method_comparison/frozen_latent.R")

BASE <- "/project/cfRNA_NormativeModeling/MixedEffectsModeling/_temp_method_comparison"
IN_DIR <- "/project/cfRNA_NormativeModeling/MixedEffectsModeling/OutriderComparison/insample_comparison"
OUT_DIR <- file.path(BASE, "cache", "latent_fits")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
split <- if (length(args) >= 1) args[1] else "cv"

folds_path <- if (split == "lobo") file.path(BASE, "cache", "lobo_folds.json") else
                                   file.path(IN_DIR, "cv_folds.json")
folds <- fromJSON(folds_path)
want <- if (length(args) >= 2) args[2] else names(folds)

hc_counts <- fread(file.path(IN_DIR, "hc_counts.csv"))
hc_genes <- hc_counts[[1]]
hc_mat <- as.matrix(hc_counts[, -1]); rownames(hc_mat) <- hc_genes

gene_len <- fread(file.path(IN_DIR, "gene_lengths.csv"))
setnames(gene_len, c("gene", "length"))
gene_len <- gene_len[match(hc_genes, gene_len$gene)]
stopifnot(identical(gene_len$gene, hc_genes))

bp <- MulticoreParam(workers = 4, stop.on.error = FALSE)

# fold key -> filename-safe tag; LOBO keys are batch names with spaces and dots
tag_of <- function(impl, fi) {
  key <- if (split == "lobo") gsub("[^A-Za-z0-9_.-]+", "_", fi) else fi
  sprintf("%s_%s%s", impl, if (split == "lobo") "lobo_" else "fold", key)
}

export <- function(tf, tag) {
  fwrite(data.table(gene = tf$genes, b = tf$b, gene_mu = tf$gene_mu, theta = tf$theta,
                    loggeomeans = tf$loggeomeans, x_center = tf$x_center,
                    l2fc_mean = tf$l2fc_mean, l2fc_sd = tf$l2fc_sd),
         file.path(OUT_DIR, paste0(tag, "_genes.csv")))
  fwrite(as.data.table(tf$E), file.path(OUT_DIR, paste0(tag, "_E.csv")))
  fwrite(as.data.table(tf$D), file.path(OUT_DIR, paste0(tag, "_D.csv")))
}

for (fi in want) {
  train_names <- folds[[fi]]$train
  for (impl in c("pca", "autoencoder")) {
    tag <- tag_of(impl, fi)
    if (file.exists(file.path(OUT_DIR, paste0(tag, "_E.csv")))) {
      cat(sprintf("skip %s (cached)\n", tag)); next
    }
    t0 <- Sys.time()
    tf <- fit_train_impl(hc_mat[, train_names, drop = FALSE], gene_len$length,
                         q = 20, implementation = impl, iterations = 5, bp = bp)
    export(tf, tag)
    cat(sprintf("%s: %d genes, %d dropped, %.0fs\n", tag, length(tf$genes), tf$n_dropped,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
cat("FIT_LATENT_DONE\n")
