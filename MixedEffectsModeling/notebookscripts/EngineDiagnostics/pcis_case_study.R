suppressPackageStartupMessages({library(glmmTMB); library(parallel); library(jsonlite)})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", grep("--file=", .args, value = TRUE))
source(file.path(dirname(normalizePath(.script_path)), "..", "..", "core", "glmm_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
WORK <- args[1]; PRIOR <- args[2]; FP_PATH <- args[3]; OUT_GENE <- args[4]; OUT_OBS <- args[5]

X <- as.matrix(read.csv(file.path(WORK, "X.csv.gz"), row.names = 1))
batch <- read.csv(file.path(WORK, "batch.csv.gz"), row.names = 1)[[1]]
Y <- read.csv(file.path(WORK, "Y_cases.csv.gz"), row.names = 1, check.names = FALSE)
meta <- read.csv(file.path(WORK, "genes_cases.csv"))
fp <- fromJSON(FP_PATH)
tau_slope <- as.numeric(fromJSON(PRIOR)$tau_slope)
sn <- sanitize_names(colnames(X)); colnames(X) <- sn
Xa <- cbind(1, X); n <- nrow(Xa); p <- ncol(Xa)
Z <- model.matrix(~ 0 + factor(batch))

fit_gene <- function(d, stage) {
  mu_fml <- as.formula(paste("y__ ~", paste(sn, collapse = " + "), "+ (1 | batch__)"))
  if (stage == "nbi_full_eb") {
    pr <- data.frame(prior = sprintf("normal(0, %.6g)", tau_slope), class = "betad", coef = sn)
    glmmTMB(mu_fml, dispformula = as.formula(paste("~", paste(sn, collapse = " + "))),
            family = nbinom2(), data = d, priors = pr)
  } else {
    glmmTMB(mu_fml, dispformula = ~ 1, family = nbinom2(), data = d)
  }
}

pcis_with <- function(mu, y, alpha, tau2) {
  r <- (y - mu) / sqrt(mu + alpha * mu^2)
  sw <- sqrt(mu / (1 + alpha * mu))
  Mw <- cbind(Xa, Z) * sw
  A <- solve(crossprod(Mw) + diag(c(rep(0, p), rep(1 / max(tau2, 1e-6), ncol(Z)))))
  h <- pmin(pmax(rowSums((Mw %*% A) * Mw), 0), 1 - 1e-8)
  list(pcis = r^2 * h / (sum(h) * (1 - h)^2), h = h, p_eff = sum(h))
}

cooks_classical <- function(mu, y, alpha) {
  r <- (y - mu) / sqrt(mu + alpha * mu^2)
  Xw <- Xa * sqrt(mu / (1 + alpha * mu))
  h <- pmin(pmax(rowSums((Xw %*% solve(crossprod(Xw))) * Xw), 0), 1 - 1e-8)
  r^2 * h / (p * (1 - h)^2)
}

own_alpha <- function(fit, Xd) {
  a <- exp(-as.numeric(Xd %*% fixef(fit)$disp))
  if (length(a) == 1) rep(a, nrow(Xd)) else a
}

work <- function(i) {
  g <- meta$gene[i]; st <- meta$stage[i]; ta <- meta$trend_alpha[i]
  y <- as.integer(round(Y[[g]]))
  d <- as.data.frame(X); d$y__ <- y; d$batch__ <- factor(batch)
  fit0 <- tryCatch(fit_gene(d, st), error = function(e) NULL)
  if (is.null(fit0)) return(NULL)
  mu0 <- as.numeric(predict(fit0, type = "response"))
  tau2_0 <- as.numeric(VarCorr(fit0)$cond[[1]][1, 1])
  Xd <- if (st == "nbi_full_eb") Xa else matrix(1, n, 1)
  a_own0 <- own_alpha(fit0, Xd)
  pt <- pcis_with(mu0, y, ta, tau2_0)
  po <- pcis_with(mu0, y, a_own0, tau2_0)
  cd <- cooks_classical(mu0, y, a_own0)
  drop <- pcis_outliers(fit0, Xa, y, batch, ta, fp$pcis_cut, fp$max_outlier_frac)
  a_own1 <- rep(NA_real_, n); ok1 <- NA; tau2_1 <- NA_real_; beta1 <- rep(NA_real_, p)
  if (length(drop) > 0) {
    fit1 <- tryCatch(fit_gene(droplevels(d[-drop, , drop = FALSE]), st), error = function(e) NULL)
    if (!is.null(fit1)) {
      ok1 <- isTRUE(is_converged(fit1, fp$beta_explode_thr, fp$tau2_max, fp$disp_intercept_max)$ok)
      a_own1 <- own_alpha(fit1, Xd); tau2_1 <- as.numeric(VarCorr(fit1)$cond[[1]][1, 1])
      beta1 <- as.numeric(fixef(fit1)$cond)
    }
  }
  beta0 <- as.numeric(fixef(fit0)$cond)
  gene_row <- data.frame(gene = g, stage = st, trend_alpha = ta, n_drop = length(drop),
    refit_ok = ok1, alpha_own_with = exp(mean(log(a_own0))), alpha_own_without = exp(mean(log(a_own1))),
    tau2_with = tau2_0, tau2_without = tau2_1, p_eff = pt$p_eff,
    max_abs_dbeta = if (length(drop) > 0) max(abs(beta1[-1] - beta0[-1])) else NA_real_,
    dintercept = if (length(drop) > 0) beta1[1] - beta0[1] else NA_real_)
  keep <- if (meta$full_obs[i]) seq_len(n) else drop
  obs <- if (length(keep)) data.frame(gene = g, idx = keep - 1L, batch = batch[keep], y = y[keep], mu = mu0[keep],
    pcis_trend = pt$pcis[keep], pcis_own = po$pcis[keep], cooks = cd[keep], h = pt$h[keep],
    alpha_own = a_own0[keep], removed = keep %in% drop) else NULL
  list(gene = gene_row, obs = obs)
}

res <- mclapply(seq_len(nrow(meta)), work, mc.cores = as.integer(fp$cores))
res <- res[!vapply(res, is.null, logical(1))]
write.csv(do.call(rbind, lapply(res, `[[`, "gene")), OUT_GENE, row.names = FALSE)
write.csv(do.call(rbind, lapply(res, `[[`, "obs")), gzfile(OUT_OBS), row.names = FALSE)
cat(sprintf("done: %d genes\n", length(res)))
