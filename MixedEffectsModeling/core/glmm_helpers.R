suppressPackageStartupMessages(library(glmmTMB))

sanitize_names <- function(names) {
  safe <- gsub("[^A-Za-z0-9_]", "_", names)
  bad <- grepl("^[^A-Za-z.]", safe)
  safe[bad] <- paste0("v", safe[bad])
  safe
}

safe_max_abs <- function(x) if (length(x) == 0) 0 else max(abs(x))
is_converged <- function(fit, beta_explode_thr, tau2_max, disp_intercept_max) {
  if (inherits(fit, "try-error")) return(list(ok = FALSE, singular = NA, tau2 = NA))
  beta_max <- safe_max_abs(c(fixef(fit)$cond[-1], fixef(fit)$disp[-1]))
  disp0 <- fixef(fit)$disp[1]
  tau2 <- as.numeric(VarCorr(fit)$cond[[1]][1, 1])
  if (isTRUE(beta_max >= beta_explode_thr)) return(list(ok = FALSE, singular = NA, tau2 = tau2))
  if (length(disp0) > 0 && isTRUE(abs(disp0) >= disp_intercept_max)) return(list(ok = FALSE, singular = NA, tau2 = tau2))
  if (isTRUE(fit$sdr$pdHess)) {
    if (isTRUE(tau2 >= tau2_max)) return(list(ok = FALSE, singular = NA, tau2 = tau2))
    return(list(ok = TRUE, singular = FALSE, tau2 = tau2))
  }
  if (isTRUE(tau2 < 1e-5)) return(list(ok = TRUE, singular = TRUE, tau2 = 0.0))
  return(list(ok = FALSE, singular = NA, tau2 = tau2))
}

# Standard errors of the dispersion fixed effects. NA when sdreport is unusable
# -- the downstream EB squeeze reads NA as SE^2=inf and returns exactly the
# lowess trend value, recovering v2's hard-fixed dispersion as a limiting case.
disp_se <- function(fit, ncoef) {
  V <- tryCatch(vcov(fit, full = TRUE), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(V)) return(rep(NA_real_, ncoef))
  idx <- grep("^disp", names(rownames(V)))
  if (length(idx) != ncoef) return(rep(NA_real_, ncoef))
  s <- suppressWarnings(sqrt(diag(V))[idx])
  ifelse(is.finite(s), s, NA_real_)
}

# PCIS -- Prior-Conditioned Impact Score. Cook-shaped but deliberately NOT Cook's
# distance, and named apart from it so the difference cannot be lost. Both of its
# departures are prior-conditioned: the variance is conditioned on the lowess
# trend prior, and the leverage on the prior-penalised mixed design.
#
#   w_i   = mu_i / (1 + alpha_trend * mu_i)          NB2 log-link IRLS weight
#   M     = [Xf Z],  P = blkdiag(0_p, I/tau^2)       fixed design + batch design
#   H     = W^1/2 M (M'WM + P)^-1 M' W^1/2
#   p_eff = tr(H)
#   PCIS_i= r_i^2 / p_eff * h_ii / (1 - h_ii)^2,  r_i = (y_i-mu_i)/sqrt(mu_i+alpha_trend*mu_i^2)
#
# Two deliberate departures from Cook's distance, both measured:
#
# 1. The variance uses the lowess TREND dispersion, not the gene's own fitted
#    one. With a freely estimated dispersion an outlier masks itself: three 20x
#    outliers on a synthetic near-Poisson gene inflated alpha 0.004 -> 0.147
#    (36x), dropping their statistic from 4.5-9.4 to 0.6-1.1 -- below threshold,
#    so nothing was flagged at ANY outlier magnitude. This is what breaks the
#    one-step deletion approximation, so PCIS is an influence heuristic, not an
#    estimate of the parameter shift under deletion.
# 2. The leverage includes the batch design Z under a ridge penalty 1/tau^2
#    (Hodges & Sargent effective df), because mu already contains the BLUP: a
#    fixed-effect-only hat matrix mixes two different models inside one
#    statistic. Measured effect: p_eff 18.4 vs p 11 (40% of the model's effective
#    complexity was being ignored) and singleton-batch leverage 0.017 -> 0.165
#    (10x). tau2 -> 0 sends the penalty to infinity, so p_eff -> p automatically.
#
# PCIS has no F reference distribution (see point 1 above), so the cut is a
# fixed constant read off an empirical null: genes regenerated from their own
# fitted (beta, gamma, tau^2) and refit under the same stage/prior. The target
# is a population-level per-observation false-alarm rate of 1e-4 -- not a
# principled level, but the first one just above the point where null-driven
# removals cross below the observed real removal rate.
#
# The null's leverage comes from model.matrix(~0+factor(batch)), so the cut is
# tied to the batch design and must be recalibrated whenever Batch_ID changes.
# 2026-07-28: 19,158 genes x 693 obs -> 2.28 (deployed 2.25).
# 2026-09-30: 1,912 genes x 679 obs, 25 batches -> 3.62. The whole move is the
# design (Moore_1 71 -> 48/7/16); a composition-matched control on the old null
# still gave 2.29.
# Observations are dropped rather than replaced by a trimmed mean, which would
# fabricate counts and bias dispersion downward. Returns indices ordered by
# decreasing influence.
pcis_outliers <- function(fit, Xa, y, batch, trend_alpha, cut, max_frac) {
  p <- ncol(Xa); n <- length(y)
  if (n - p < 4) return(integer(0))
  mu <- tryCatch(as.numeric(predict(fit, type = "response")), error = function(e) NULL)
  if (is.null(mu) || !all(is.finite(mu)) || !is.finite(trend_alpha) || trend_alpha <= 0) return(integer(0))
  alpha <- trend_alpha
  r <- (y - mu) / sqrt(mu + alpha * mu^2)
  sw <- sqrt(mu / (1 + alpha * mu))

  tau2 <- tryCatch(as.numeric(VarCorr(fit)$cond[[1]][1, 1]), error = function(e) NA_real_)
  Z <- tryCatch(model.matrix(~ 0 + factor(batch)), error = function(e) NULL)
  h <- NULL
  if (!is.null(Z) && isTRUE(is.finite(tau2))) {
    Mw <- cbind(Xa, Z) * sw
    P <- diag(c(rep(0, p), rep(1 / max(tau2, 1e-6), ncol(Z))))
    A <- tryCatch(solve(crossprod(Mw) + P), error = function(e) NULL)
    if (!is.null(A)) h <- rowSums((Mw %*% A) * Mw)
  }
  if (is.null(h)) {                       # fall back to the fixed-effect design
    Xw <- Xa * sw
    A <- tryCatch(solve(crossprod(Xw)), error = function(e) NULL)
    if (is.null(A)) return(integer(0))
    h <- rowSums((Xw %*% A) * Xw)
  }
  h <- pmin(pmax(h, 0), 1 - 1e-8)
  p_eff <- sum(h)
  if (!is.finite(p_eff) || p_eff < 1 || n - p_eff < 4) return(integer(0))

  PCIS <- r^2 * h / (p_eff * (1 - h)^2)
  over <- which(is.finite(PCIS) & PCIS > cut)
  if (length(over) == 0) return(integer(0))
  n_keep <- min(length(over), floor(max_frac * n))
  if (n_keep < 1) return(integer(0))
  over[order(-PCIS[over])][seq_len(n_keep)]
}

# Fits ONE stage for ONE gene. Caller (glmm_fit.R) drives the demotion order.
# v3: stages are "nbi_full_eb" (dispersion regressed on covariates) and "nbi_intercept_eb"
# (dispersion intercept only). "nbi_disp_intercept" is gone -- with a
# properly-scaled slope prior it was no longer a distinct model, and "nbi_intercept_eb"
# now occupies its structural position. The dispersion INTERCEPT is left
# unpenalized here; it is squeezed toward the lowess trend analytically
# downstream (eb_shrinkage.squeeze_log_theta). tau_slope is the per-covariate EB
# prior sd for the dispersion SLOPES, estimated by a --mode calib run.
fit_stage_gene <- function(y, safe_names, X, batch, stage, tau_slope, trend_alpha,
                           beta_explode_thr, tau2_max, disp_intercept_max,
                           pcis_cut, max_outlier_frac) {
  df <- as.data.frame(X); colnames(df) <- safe_names
  df$y__ <- as.integer(round(y))
  df$batch__ <- factor(batch)

  mu_fml <- as.formula(paste("y__ ~", paste(safe_names, collapse = " + "), "+ (1 | batch__)"))
  disp_fml <- switch(stage,
    nbi_full_eb = as.formula(paste("~", paste(safe_names, collapse = " + "))),
    nbi_intercept_eb = as.formula("~ 1"),
    stop(sprintf("fit_stage_gene: unknown stage '%s'", stage)))
  n_disp <- if (stage == "nbi_full_eb") length(safe_names) + 1L else 1L

  pr <- NULL
  if (stage == "nbi_full_eb" && !is.null(tau_slope)) {
    pr <- data.frame(prior = sprintf("normal(0, %.6g)", tau_slope),
                     class = "betad", coef = safe_names)
  }
  run <- function(d) tryCatch({
    if (is.null(pr)) glmmTMB(mu_fml, dispformula = disp_fml, family = nbinom2(), data = d)
    else glmmTMB(mu_fml, dispformula = disp_fml, family = nbinom2(), data = d, priors = pr)
  }, error = function(e) structure(conditionMessage(e), class = "try-error"))

  fit <- run(df)
  if (inherits(fit, "try-error")) {
    return(list(stage = stage, ok = FALSE, singular = NA, tau2 = NA, n_outliers = 0L,
               outlier_refit_failed = FALSE, mu_coef = numeric(0), disp_coef = numeric(0),
               disp_se = rep(NA_real_, n_disp), fail_reason = as.character(fit)))
  }

  conv <- is_converged(fit, beta_explode_thr, tau2_max, disp_intercept_max)
  n_out <- 0L; refit_failed <- FALSE
  # PCIS is computed even when the first fit failed the convergence gate --
  # outliers can be the cause of the failure, so the refit gets a chance.
  drop_idx <- pcis_outliers(fit, cbind(1, as.matrix(X)), df$y__, batch,
                            trend_alpha, pcis_cut, max_outlier_frac)
  if (length(drop_idx) > 0) {
    # droplevels: HC has 5 singleton batches, so removing one observation can
    # empty a random-effect level entirely.
    fit2 <- run(droplevels(df[-drop_idx, , drop = FALSE]))
    conv2 <- if (inherits(fit2, "try-error")) list(ok = FALSE) else
      is_converged(fit2, beta_explode_thr, tau2_max, disp_intercept_max)
    if (isTRUE(conv2$ok)) {
      fit <- fit2; conv <- conv2; n_out <- length(drop_idx)
    } else {
      refit_failed <- TRUE
    }
  }

  list(stage = stage, ok = conv$ok, singular = conv$singular, tau2 = conv$tau2,
      n_outliers = n_out, outlier_refit_failed = refit_failed,
      mu_coef = as.numeric(fixef(fit)$cond), disp_coef = as.numeric(fixef(fit)$disp),
      disp_se = disp_se(fit, n_disp),
      fail_reason = if (isTRUE(conv$ok)) "" else "not_converged_or_explosion_or_tau2_bound")
}

# Shared-beta pooled GLM (route "pool") + batch random intercept. Unused this
# round (pooling deferred until nz_a_max is picked) -- kept so the file stays
# a complete, working unit.
fit_pooled_glmm <- function(Y_block, X, batch, mean_hc, eps, rare_overdisp_thr) {
  n_hc <- nrow(X); n_g <- ncol(Y_block)
  safe_names <- sanitize_names(colnames(X))
  sample_idx <- rep(seq_len(n_hc), n_g)
  gene_idx <- rep(seq_len(n_g), each = n_hc)
  df <- as.data.frame(X[sample_idx, , drop = FALSE]); colnames(df) <- safe_names
  df$y__ <- as.integer(round(Y_block[cbind(sample_idx, gene_idx)]))
  df$batch__ <- factor(batch[sample_idx])
  df$off__ <- log(mean_hc[gene_idx] + eps)
  mu_fml <- as.formula(paste("y__ ~ offset(off__) +", paste(safe_names, collapse = " + "), "+ (1 | batch__)"))

  mult_clip <- function(beta) list(
    mult_lo = as.numeric(quantile(exp(X %*% beta[-1]), 0.001)),
    mult_hi = as.numeric(quantile(exp(X %*% beta[-1]), 0.999)))

  pool_tau2 <- function(fit) as.numeric(VarCorr(fit)$cond[[1]][1, 1])

  fit_pois <- tryCatch(glmmTMB(mu_fml, family = poisson(), data = df), error = function(e) NULL)
  if (is.null(fit_pois)) return(list(ok = FALSE))
  ratio <- sum(residuals(fit_pois, type = "pearson")^2) / df.residual(fit_pois)
  if (ratio <= rare_overdisp_thr) {
    beta <- as.numeric(fixef(fit_pois)$cond)
    mc <- mult_clip(beta)
    return(list(ok = TRUE, family = "poisson", beta = beta, alpha = NA, tau2 = pool_tau2(fit_pois),
               overdisp_ratio = ratio, mult_lo = mc$mult_lo, mult_hi = mc$mult_hi))
  }
  fit_nb <- tryCatch(glmmTMB(mu_fml, family = nbinom2(), data = df), error = function(e) NULL)
  if (is.null(fit_nb)) return(list(ok = FALSE))
  beta <- as.numeric(fixef(fit_nb)$cond)
  mc <- mult_clip(beta)
  list(ok = TRUE, family = "negbin", beta = beta,
      alpha = exp(-fixef(fit_nb)$disp[["(Intercept)"]]),
      tau2 = pool_tau2(fit_nb), overdisp_ratio = ratio, mult_lo = mc$mult_lo, mult_hi = mc$mult_hi)
}
