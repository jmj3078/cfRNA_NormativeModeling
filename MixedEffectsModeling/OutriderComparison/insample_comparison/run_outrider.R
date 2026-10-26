suppressPackageStartupMessages({
  library(OUTRIDER)
  library(data.table)
})

dir <- commandArgs(trailingOnly = TRUE)[1]  # config.OUTRIDER_COMPARISON_DIR
counts <- fread(file.path(dir, "hc_counts.csv"))
genes <- counts[[1]]
mat <- as.matrix(counts[, -1])
rownames(mat) <- genes
mode(mat) <- "integer"

meta <- fread(file.path(dir, "hc_meta.csv"))
stopifnot(identical(colnames(mat), meta$sample))

gene_len <- fread(file.path(dir, "gene_lengths.csv"))
setnames(gene_len, c("gene", "length"))
gene_len <- gene_len[match(genes, gene_len$gene)]
stopifnot(identical(gene_len$gene, genes))

ods <- OutriderDataSet(countData = mat)
mcols(ods)$basepairs <- gene_len$length
ods <- filterExpression(ods, filterGenes = TRUE, fpkmCutoff = 1)
ods <- estimateSizeFactors(ods)

set.seed(42)
q <- 20
bp <- MulticoreParam(workers = 4, stop.on.error = FALSE)

# Neither controlForConfounders() nor fit() has a partial-success mode -- ANY gene's
# per-gene NB optimization failing (L-BFGS-B non-finite) aborts the WHOLE call. Parse the
# raised bplist error for the failing gene indices, drop those genes, and retry (looping
# in case a new gene fails once the offending ones are gone).
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
    cat(sprintf("dropping %d gene(s) that fail NB optimization: %s\n", length(bad),
               paste(rownames(ods)[bad], collapse = ", ")))
    n_dropped <- n_dropped + length(bad)
    ods <- ods[-bad, ]
  }
}

r1 <- retry_drop(controlForConfounders, ods, q = q, iterations = 5, BPPARAM = bp)
ods <- r1$ods
r2 <- retry_drop(fit, ods, BPPARAM = bp)
ods <- r2$ods
n_dropped <- r1$n_dropped + r2$n_dropped
cat(sprintf("total genes dropped for NB fit failure: %d (out of filtered set)\n", n_dropped))
ods <- computePvalues(ods, alternative = "two.sided")
ods <- computeZscores(ods)

pv <- as.matrix(assay(ods, "pValue"))
padj <- as.matrix(assay(ods, "padjust"))
z <- as.matrix(assay(ods, "zScore"))

aberrant <- padj < 0.05
per_sample_aberrant <- colSums(aberrant, na.rm = TRUE)
per_sample_mean_absz <- colMeans(abs(z), na.rm = TRUE)
per_sample_pos_frac <- colMeans(z > 0, na.rm = TRUE)

out <- data.table(sample = colnames(mat), batch = meta$batch,
                   n_aberrant = per_sample_aberrant,
                   mean_absz = per_sample_mean_absz,
                   pos_frac = per_sample_pos_frac)
fwrite(out, file.path(dir, "outrider_per_sample.csv"))
saveRDS(ods, file.path(dir, "ods_fit.rds"))
cat("DONE q=", q, "\n")
