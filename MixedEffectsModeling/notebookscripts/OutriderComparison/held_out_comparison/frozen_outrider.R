# Genuine held-out OUTRIDER scoring: fit autoencoder+NB (E, D, b, theta, mu,
# size-factor reference) on a TRAIN-only sample set, then project NEW test
# samples through the frozen model with zero test-data leakage into any fitted
# parameter. Standard OUTRIDER (controlForConfounders/fit/computeZscores) has
# no such train/apply split -- test samples always sit inside the same fit
# matrix (see run_outrider_lobo.R / run_outrider_cv.R in ../insample_comparison).
# All formulas below are copied 1:1 from OUTRIDER's own internal (non-exported)
# functions (see scratchpad/outrider_src*.txt for the dumped source), so a
# frozen-model scoring of the TRAIN samples themselves must reproduce
# OUTRIDER's native in-sample output exactly -- see validate_frozen.R.

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(data.table)
})

.ns <- getNamespace("OUTRIDER")
.predictMatC <- get("predictMatC", envir = .ns)
.log2fc <- get("log2fc", envir = .ns)

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

# Replicates OUTRIDER's internal x(ods): centered log-normalized counts.
# counts_mat: genes x samples. sf: per-sample size factor. center: per-gene
# offset to subtract; NULL derives it from this call's own samples (do this
# only for the TRAIN fit -- test scoring must always pass the frozen
# train center).
compute_x <- function(counts_mat, sf, center = NULL) {
  k <- t(counts_mat)                    # samples x genes
  x0 <- log((1 + k) / sf)                # sf recycled per-row (per-sample)
  if (is.null(center)) center <- colMeans(x0)
  list(x = t(t(x0) - center), center = center)
}

# Replicates estimateSizeFactors()'s median-of-ratios against a frozen
# per-gene log-geometric-mean reference (loggeomeans), so test size factors
# never depend on other test samples.
compute_size_factors <- function(counts_mat, loggeomeans) {
  apply(counts_mat, 2, function(cnts) {
    exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
  })
}

# Fit the full OUTRIDER pipeline (autoencoder + per-gene NB) on TRAIN samples
# only, and freeze every parameter needed to score unseen samples later.
fit_train <- function(train_mat, gene_len, q = 20, iterations = 5, bp = SerialParam()) {
  ods <- OutriderDataSet(countData = train_mat)
  mcols(ods)$basepairs <- gene_len
  ods <- filterExpression(ods, filterGenes = TRUE, fpkmCutoff = 1)
  ods <- estimateSizeFactors(ods)
  set.seed(42)
  r1 <- retry_drop(controlForConfounders, ods, q = q, iterations = iterations, BPPARAM = bp)
  r2 <- retry_drop(fit, r1$ods)
  ods <- r2$ods

  xc <- compute_x(counts(ods, normalized = FALSE), sizeFactors(ods))
  l2fc_train <- .log2fc(ods)

  list(
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
    n_dropped = r1$n_dropped + r2$n_dropped,
    ods_train = ods  # kept only for the in-sample self-check; drop before saving to disk
  )
}

# Score NEW samples through the frozen train_fit. test_mat must contain at
# least train_fit$genes (extra genes are ignored). Returns z/pval/normF/l2fc
# matrices, samples x genes.
score_frozen <- function(train_fit, test_mat) {
  genes <- train_fit$genes
  k_test <- test_mat[genes, , drop = FALSE]

  sf_test <- compute_size_factors(k_test, train_fit$loggeomeans)
  xc_test <- compute_x(k_test, sf_test, center = train_fit$x_center)

  # samples x genes, frozen decoder -- test data never touched E/D/b fitting
  normF_test_sxg <- .predictMatC(xc_test$x, train_fit$E, train_fit$D, train_fit$b, sf_test)
  normF_test <- t(normF_test_sxg)  # genes x samples, matches counts orientation
  dimnames(normF_test) <- dimnames(k_test)

  l2fc_test <- log2(k_test + 1) - log2(normF_test + 1)
  z_test <- (l2fc_test - train_fit$l2fc_mean) / train_fit$l2fc_sd

  mu_nb <- train_fit$gene_mu * normF_test  # genes x samples, NB mean per pVal()
  theta_mat <- matrix(train_fit$theta, nrow = length(genes), ncol = ncol(k_test))
  pless <- pnbinom(k_test, size = theta_mat, mu = mu_nb)
  dval <- dnbinom(k_test, size = theta_mat, mu = mu_nb)
  pval_test <- 2 * pmin(0.5, pless, 1 - pless + dval)
  pval_test <- matrix(pval_test, nrow = nrow(k_test), ncol = ncol(k_test), dimnames = dimnames(k_test))

  list(z = t(z_test), pval = t(pval_test), normF = t(normF_test), l2fc = t(l2fc_test),
       counts = t(k_test), mu_nb = t(mu_nb))
}
