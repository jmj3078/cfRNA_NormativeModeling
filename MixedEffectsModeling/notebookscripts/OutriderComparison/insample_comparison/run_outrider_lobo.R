suppressPackageStartupMessages({
  library(OUTRIDER)
  library(data.table)
  library(jsonlite)
})

dir <- commandArgs(trailingOnly = TRUE)[1]  # config.OUTRIDER_COMPARISON_DIR

hc_counts <- fread(file.path(dir, "hc_counts.csv"))
hc_genes <- hc_counts[[1]]
hc_mat <- as.matrix(hc_counts[, -1]); rownames(hc_mat) <- hc_genes
hc_meta <- fread(file.path(dir, "hc_meta.csv"))
stopifnot(identical(colnames(hc_mat), hc_meta$sample))

dis_counts <- fread(file.path(dir, "disease_counts.csv"))
dis_genes <- dis_counts[[1]]
dis_mat <- as.matrix(dis_counts[, -1]); rownames(dis_mat) <- dis_genes
stopifnot(identical(hc_genes, dis_genes))

gene_len <- fread(file.path(dir, "gene_lengths.csv"))
setnames(gene_len, c("gene", "length"))
gene_len <- gene_len[match(hc_genes, gene_len$gene)]
stopifnot(identical(gene_len$gene, hc_genes))

test_meta <- fromJSON(file.path(dir, "lobo_test_meta.json"))
batches <- names(test_meta)

q <- 20
bp <- MulticoreParam(workers = 4, stop.on.error = FALSE)

# Same retry-drop as run_outrider.R: neither controlForConfounders() nor fit() has a
# partial-success mode -- ANY gene's L-BFGS-B divergence aborts the whole call.
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

out_rows <- list()

for (b in batches) {
  safe <- gsub(" ", "_", b)
  cat(sprintf("=== %s ===\n", b))
  test_names <- test_meta[[b]]$test_names
  test_is_hc <- test_meta[[b]]$test_is_hc

  train_hc_names <- hc_meta$sample[hc_meta$batch != b]
  train_mat <- hc_mat[, train_hc_names, drop = FALSE]

  test_from_hc <- intersect(test_names, colnames(hc_mat))
  test_from_dis <- intersect(test_names, colnames(dis_mat))
  test_mat <- cbind(hc_mat[, test_from_hc, drop = FALSE], dis_mat[, test_from_dis, drop = FALSE])
  test_mat <- test_mat[, test_names, drop = FALSE]  # restore original order

  mat <- cbind(train_mat, test_mat)
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
  z <- t(as.matrix(assay(ods, "zScore")))  # samples x genes

  z_test <- z[test_names, , drop = FALSE]
  saveRDS(z_test, file.path(dir, paste0("z_test_", safe, ".rds")))
  out_rows[[b]] <- data.table(batch = b, sample = test_names, is_hc = test_is_hc)
}

fwrite(rbindlist(out_rows), file.path(dir, "outrider_lobo_membership.csv"))
cat("ALL DONE\n")
