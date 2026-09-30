suppressPackageStartupMessages({
  library(optparse); library(parallel); library(jsonlite)
})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", grep("--file=", .args, value = TRUE))
source(file.path(dirname(normalizePath(.script_path)), "glmm_helpers.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--x", type = "character"), make_option("--y", type = "character"),
  make_option("--batch", type = "character"), make_option("--genes", type = "character"),
  make_option("--trend", type = "character"), make_option("--mode", type = "character", default = "cascade"),
  make_option("--disp-prior", type = "character", default = ""),
  make_option("--fit-params", type = "character", default = ""),
  make_option("--out", type = "character"), make_option("--chunk-size", type = "integer", default = 200),
  make_option("--cores", type = "integer", default = min(parallel::detectCores() - 1, 12))
)))

X <- as.matrix(read.csv(opt$x, row.names = 1))
Y <- read.csv(opt$y, row.names = 1)
batch <- read.csv(opt$batch, row.names = 1)[[1]]
gene_meta <- read.csv(opt$genes)  # columns: gene, [stage] (stage only needed for fixed_stage mode)
# The trend is itself estimated FROM a calibration run, so calib mode has none yet:
# alpha_of returns NA there, trend_alpha is reported as NA and pcis_outliers
# no-ops (a calibration run exists only to supply hyperparameters, not deployed fits).
alpha_of <- function(mean_y) NA_real_
if (opt$mode != "calib" && nzchar(opt$trend) && file.exists(opt$trend)) {
  trend <- fromJSON(opt$trend)
  alpha_of <- function(mean_y) {
    lm <- log(max(mean_y, 1e-8))
    s <- exp(approx(trend$lowess_logmu, trend$lowess_logsigma, xout = lm, rule = 2)$y)
    min(max(s, trend$alpha_floor), trend$alpha_cap)
  }
}
safe_names <- sanitize_names(colnames(X))
colnames(X) <- safe_names
# Fallback defaults, identical to config.FIT_PARAMS. They exist only so the
# script still runs standalone without --fit-params; config.py is the source.
FP <- list(beta_explode_thr = 3.0, tau2_max = 9.0, disp_intercept_max = 10.0,
           pcis_cut = 3.62, max_outlier_frac = 0.05)
if (nzchar(opt$`fit-params`) && file.exists(opt$`fit-params`)) {
  loaded <- fromJSON(opt$`fit-params`)
  for (k in names(FP)) if (!is.null(loaded[[k]])) FP[[k]] <- loaded[[k]]
  if (!is.null(loaded$chunk_size)) opt$`chunk-size` <- as.integer(loaded$chunk_size)
  if (!is.null(loaded$cores)) opt$cores <- as.integer(loaded$cores)
}
BETA_EXPLODE_THR <- FP$beta_explode_thr
TAU2_MAX <- FP$tau2_max
DISP_INTERCEPT_MAX <- FP$disp_intercept_max
PCIS_CUT <- FP$pcis_cut
MAX_OUTLIER_FRAC <- FP$max_outlier_frac

# EB prior sd for the dispersion SLOPES, one value per covariate, estimated by a
# --mode calib run (eb_shrinkage.estimate_slope_prior). NULL in calib mode and
# whenever the json is absent, so the calibration run itself fits with no
# dispersion prior at all.
TAU_SLOPE <- NULL
if (opt$mode != "calib" && nzchar(opt$`disp-prior`) && file.exists(opt$`disp-prior`)) {
  TAU_SLOPE <- as.numeric(fromJSON(opt$`disp-prior`)$tau_slope)
  stopifnot(length(TAU_SLOPE) == length(safe_names))
}

done_genes <- character(0)
if (file.exists(opt$out)) done_genes <- read.csv(opt$out)$gene

fit_stage <- function(y, stage, trend_alpha) fit_stage_gene(
  y, safe_names, X, batch, stage, TAU_SLOPE, trend_alpha,
  BETA_EXPLODE_THR, TAU2_MAX, DISP_INTERCEPT_MAX, PCIS_CUT, MAX_OUTLIER_FRAC)

# v3: 2-stage cascade (nbi_full_eb -> nbi_intercept_eb). force-return at nbi_intercept_eb regardless of ok --
# a gene that fails there is genuinely unmodelable and comes back ok=FALSE
# (caller routes it to "excluded"). Every stage's own reject reason is preserved.
fit_one_cascade <- function(g) {
  y <- as.numeric(Y[[g]])
  alpha_g <- alpha_of(mean(y))
  reasons <- list(nbi_full_eb = "", nbi_intercept_eb = "")
  for (stage in c("nbi_full_eb", "nbi_intercept_eb")) {
    r <- fit_stage(y, stage, alpha_g)
    reasons[[stage]] <- if (isTRUE(r$ok)) "" else r$fail_reason
    if (isTRUE(r$ok) || stage == "nbi_intercept_eb") {
      gc()
      return(c(list(gene = g, trend_alpha = alpha_g,
                    nbi_full_eb_reject_reason = reasons$nbi_full_eb,
                    nbi_intercept_eb_reject_reason = reasons$nbi_intercept_eb), r))
    }
  }
}

# calib: stage nbi_full_eb only, no dispersion prior, used to estimate TAU_SLOPE.
# fixed_stage: refit the stage the full run already chose (CV folds).
fit_one_single <- function(g, stage) {
  y <- as.numeric(Y[[g]])
  alpha_g <- alpha_of(mean(y))
  r <- fit_stage(y, stage, alpha_g)
  gc()
  c(list(gene = g, trend_alpha = alpha_g,
        nbi_full_eb_reject_reason = "", nbi_intercept_eb_reject_reason = ""), r)
}

worker <- switch(opt$mode,
  cascade = fit_one_cascade,
  calib = function(g) fit_one_single(g, "nbi_full_eb"),
  fixed_stage = function(g) fit_one_single(g, gene_meta$stage[gene_meta$gene == g]),
  stop(sprintf("unknown mode '%s'", opt$mode)))

genes_todo <- setdiff(gene_meta$gene, done_genes)
chunks <- split(genes_todo, ceiling(seq_along(genes_todo) / opt$`chunk-size`))

t0 <- Sys.time()
n_ok_cum <- length(done_genes)
n_total_cum <- length(done_genes)
for (i in seq_along(chunks)) {
  results <- mclapply(chunks[[i]], worker, mc.cores = opt$cores)
  rows <- lapply(results, function(r) {
    p <- 11  # 1 intercept + 10 covariates
    pad <- function(v) c(v, rep(NA, p - length(v)))[1:p]
    mu_padded <- pad(r$mu_coef)
    disp_padded <- pad(r$disp_coef)
    se_padded <- pad(r$disp_se)
    row <- c(list(gene = r$gene, stage = r$stage, ok = r$ok, singular = r$singular,
                 tau2 = r$tau2, trend_alpha = r$trend_alpha,
                 n_outliers = r$n_outliers, outlier_refit_failed = r$outlier_refit_failed,
                 fail_reason = r$fail_reason,
                 nbi_full_eb_reject_reason = r$nbi_full_eb_reject_reason,
                 nbi_intercept_eb_reject_reason = r$nbi_intercept_eb_reject_reason))
    for (j in seq_len(p)) {
      row[[paste0("mu_coef_", j-1)]] <- mu_padded[j]
      row[[paste0("disp_coef_", j-1)]] <- disp_padded[j]
      row[[paste0("disp_se_", j-1)]] <- se_padded[j]
    }
    row
  })
  df <- do.call(rbind, lapply(rows, as.data.frame))
  write.table(df, opt$out, sep = ",", append = file.exists(opt$out), col.names = !file.exists(opt$out), row.names = FALSE)
  gc()

  n_ok_cum <- n_ok_cum + sum(df$ok, na.rm = TRUE)
  n_total_cum <- n_total_cum + nrow(df)
  elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  eta_min <- (length(chunks) - i) * (elapsed_min / i)
  stage_counts <- paste(sprintf("%s=%d", names(table(df$stage)), table(df$stage)), collapse = ",")
  cat(sprintf("[%s] chunk %d/%d done (%d genes, ok_rate=%.2f, %s, out_genes=%d) | elapsed=%.1fmin eta=%.1fmin | cum_ok_rate=%.3f (%d/%d)\n",
             format(Sys.time(), "%H:%M:%S"), i, length(chunks), nrow(df), mean(df$ok, na.rm = TRUE), stage_counts,
             sum(df$n_outliers > 0, na.rm = TRUE), elapsed_min, eta_min, n_ok_cum / n_total_cum, n_ok_cum, n_total_cum))
}
cat("DONE\n")
