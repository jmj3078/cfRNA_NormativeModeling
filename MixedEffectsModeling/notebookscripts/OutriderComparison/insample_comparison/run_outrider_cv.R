suppressPackageStartupMessages({
  library(OUTRIDER)
  library(data.table)
  library(jsonlite)
})

dir <- commandArgs(trailingOnly = TRUE)[1]  # config.OUTRIDER_COMPARISON_DIR

hc_counts <- fread(file.path(dir, "hc_counts.csv"))
hc_genes <- hc_counts[[1]]
hc_mat <- as.matrix(hc_counts[, -1]); rownames(hc_mat) <- hc_genes

gene_len <- fread(file.path(dir, "gene_lengths.csv"))
setnames(gene_len, c("gene", "length"))
gene_len <- gene_len[match(hc_genes, gene_len$gene)]
stopifnot(identical(gene_len$gene, hc_genes))

folds <- fromJSON(file.path(dir, "cv_folds.json"))
q <- 20
bp <- MulticoreParam(workers = 4, stop.on.error = FALSE)

retry_drop <- function(f, ods, ...) {
  n_dropped <- 0
  repeat {
    msg <- NULL
    res <- tryCatch(list(ok = TRUE, val = f(ods, ...)),
                    error = function(e) { msg <<- conditionMessage(e); list(ok = FALSE) })
    if (res$ok) return(list(ods = res$val, n_dropped = n_dropped))
    m <- regmatches(msg, regexpr("element index: [0-9, ]+", msg))
    idx_str <- sub("element index: ", "", m)
    bad <- as.integer(trimws(strsplit(idx_str, ",")[[1]]))
    bad <- bad[!is.na(bad)]
    if (length(bad) == 0) stop(msg)
    cat(sprintf("  dropping %d gene(s) that fail NB optimization: %s\n", length(bad),
               paste(rownames(ods)[bad], collapse = ", ")))
    n_dropped <- n_dropped + length(bad)
    ods <- ods[-bad, ]
  }
}

for (fi in names(folds)) {
  cat(sprintf("=== fold %s ===\n", fi))
  train_names <- folds[[fi]]$train
  test_names <- folds[[fi]]$test
  mat <- hc_mat[, c(train_names, test_names), drop = FALSE]

  ods <- OutriderDataSet(countData = mat)
  mcols(ods)$basepairs <- gene_len$length
  ods <- filterExpression(ods, filterGenes = TRUE, fpkmCutoff = 1)
  ods <- estimateSizeFactors(ods)
  set.seed(42)
  r1 <- retry_drop(controlForConfounders, ods, q = q, iterations = 5, BPPARAM = bp)
  r2 <- retry_drop(fit, r1$ods, BPPARAM = bp)
  ods <- r2$ods
  cat(sprintf("  n_genes final=%d, dropped=%d\n", nrow(ods), r1$n_dropped + r2$n_dropped))
  ods <- computeZscores(ods)

  z <- t(as.matrix(assay(ods, "zScore")))
  mu <- t(normalizationFactors(ods))
  th <- theta(ods)
  y <- t(counts(ods))

  z_test <- z[test_names, , drop = FALSE]
  mu_test <- mu[test_names, , drop = FALSE]
  y_test <- y[test_names, , drop = FALSE]
  saveRDS(list(z = z_test, mu = mu_test, y = y_test, theta = th, genes = rownames(ods)),
          file.path(dir, paste0("cv_fold", fi, "_full.rds")))
}
cat("ALL DONE\n")
