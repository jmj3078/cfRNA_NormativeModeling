suppressPackageStartupMessages(suppressWarnings(library(gamlss)))

# ── Column name sanitization ──────────────────────────────────────
# Covariate names from Python (e.g. "log(Total Reads)") contain special
# characters that R formula parser rejects. We replace them with safe names
# and store the mapping in a global for reuse.

sanitize_names <- function(names) {
  safe <- gsub("[^A-Za-z0-9_]", "_", names)
  # R identifiers must start with a letter or dot, not underscore/digit
  bad <- grepl("^[^A-Za-z.]", safe)
  safe[bad] <- paste0("v", safe[bad])
  safe
}

# ── NBI RQR on new data ───────────────────────────────────────────
# gamlss::qresid() only works on fitted objects (training data).
# For test-set RQR we compute manually via pNBI CDF.

rqr_nbi <- function(y, mu, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- as.integer(round(y))
  a <- ifelse(y > 0, pNBI(y - 1L, mu = mu, sigma = sigma), 0)
  b <- pNBI(y,        mu = mu, sigma = sigma)
  # Clamp to avoid Inf from qnorm
  a <- pmax(pmin(a, 1 - 1e-8), 1e-8)
  b <- pmax(pmin(b, 1 - 1e-8), 1e-8)
  # Guard against floating-point inversion
  lo <- pmin(a, b); hi <- pmax(a, b)
  u  <- runif(length(y), min = lo, max = hi)
  qnorm(u)
}

# ── Main fitting function ─────────────────────────────────────────
#
# Args
#   y_train, y_test : numeric vectors of raw counts
#   X_train, X_test : numeric matrices (n x p), same column order
#   seed            : integer for RQR reproducibility
#   n_cyc           : max gamlss backfitting cycles
#
# Returns a named list:
#   z       : RQR z-scores on test set (numeric vector, may contain NA)
#   mu_test : predicted mu on test set
#   sigma_test : predicted sigma (1/theta) on test set
#   success : logical
#   msg     : error message if success == FALSE

# ── Internal: fit NBI on a data frame and return the gamlss object ──
.fit_nbi <- function(df_tr, mu_fml, sigma_fml, n_cyc) {
  tryCatch(
    gamlss(
      formula       = mu_fml,
      sigma.formula = sigma_fml,
      family        = NBI(),
      data          = df_tr,
      control       = gamlss.control(n.cyc = n_cyc, trace = FALSE)
    ),
    error = function(e) e
  )
}

# ── NBI main fitting function (with iterative outlier removal) ─────
#
# Additional args vs. original:
#   outlier_z : |z_train| threshold to flag and remove outliers (default 5)
#   max_iter  : maximum refinement iterations (default 3)
#
# Returns a named list:
#   z, mu_test, sigma_test, success, msg
#   n_removed : total training samples removed across all iterations

fit_gamlss_gene <- function(y_train, y_test, X_train, X_test,
                             seed = NULL, n_cyc = 50,
                             outlier_z = 5.0, max_iter = 2L,
                             max_remove_frac = 0.05) {
  n_te       <- length(y_test)
  n_tr_orig  <- length(y_train)
  safe_names <- sanitize_names(colnames(X_train))
  mu_rhs     <- paste(safe_names, collapse = " + ")
  mu_fml     <- as.formula(paste("y__ ~", mu_rhs))
  sigma_fml  <- as.formula(paste("~",     mu_rhs))

  df_tr <- as.data.frame(X_train)
  colnames(df_tr) <- safe_names
  df_tr$y__ <- as.integer(round(y_train))

  df_te <- as.data.frame(X_test)
  colnames(df_te) <- safe_names

  na_result <- list(z = rep(NA_real_, n_te), mu_test = rep(NA_real_, n_te),
                    sigma_test = rep(NA_real_, n_te),
                    success = FALSE, msg = "", n_removed = 0L)

  keep      <- rep(TRUE, n_tr_orig)
  n_removed <- 0L

  for (iter in seq_len(max_iter)) {
    fit <- .fit_nbi(df_tr[keep, , drop = FALSE], mu_fml, sigma_fml, n_cyc)
    if (inherits(fit, "error")) {
      na_result$msg <- conditionMessage(fit)
      return(na_result)
    }

    # Training z-scores to identify outliers
    pred_tr <- tryCatch(
      predictAll(fit, newdata = df_tr[keep, , drop = FALSE],
                 type = "response", data = df_tr[keep, , drop = FALSE]),
      error = function(e) e
    )
    if (inherits(pred_tr, "error")) break

    z_tr    <- rqr_nbi(df_tr$y__[keep],
                        mu    = as.numeric(pred_tr$mu),
                        sigma = as.numeric(pred_tr$sigma))
    outlier <- is.finite(z_tr) & (abs(z_tr) > outlier_z)

    if (!any(outlier)) break  # converged — no more outliers

    # Safety cap: if a single iteration would remove too many samples, stop
    if (sum(outlier) / n_tr_orig > max_remove_frac) break

    idx_keep  <- which(keep)
    keep[idx_keep[outlier]] <- FALSE
    n_removed <- n_removed + sum(outlier)
  }

  # Final prediction on test set
  pred_te <- tryCatch(
    predictAll(fit, newdata = df_te, type = "response",
               data = df_tr[keep, , drop = FALSE]),
    error = function(e) e
  )
  if (inherits(pred_te, "error")) {
    na_result$msg <- conditionMessage(pred_te)
    return(na_result)
  }

  mu_te    <- as.numeric(pred_te$mu)
  sigma_te <- as.numeric(pred_te$sigma)

  list(z          = rqr_nbi(y_test, mu = mu_te, sigma = sigma_te, seed = seed),
       mu_test    = mu_te,
       sigma_test = sigma_te,
       success    = TRUE,
       msg        = "",
       n_removed  = n_removed)
}


# ══════════════════════════════════════════════════════════════════
# ZINB (Zero-Inflated NBI) — mixture model
#
# P(Y = 0)   = nu + (1-nu) * f_NBI(0 | mu, sigma)   [zero inflation]
# P(Y = k>0) = (1-nu) * f_NBI(k | mu, sigma)
#
# The conditional distribution given Y>0 is identical to ZANBI:
#   P(Y=k | Y>0) ∝ f_NBI(k) / (1 - f_NBI(0))
# so rqr_count_cond() is reused unchanged.
# ══════════════════════════════════════════════════════════════════

# ── Binary RQR for ZINB ──────────────────────────────────────────
# P(Y=0) = nu + (1-nu)*pNBI(0), so Binary_Z uses this combined zero prob.
rqr_binary_zinb <- function(y, mu, sigma, nu, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p0_nbi  <- pNBI(0L, mu = mu, sigma = sigma)
  p_zero  <- pmax(pmin(nu + (1 - nu) * p0_nbi, 1 - 1e-8), 1e-8)

  detected <- as.integer(round(y)) > 0
  a  <- ifelse(detected, p_zero, 0)
  b  <- ifelse(detected, 1, p_zero)
  lo <- pmax(pmin(a, 1 - 1e-8), 1e-8)
  hi <- pmax(pmin(b, 1 - 1e-8), 1e-8)
  lo <- pmin(lo, hi); hi <- pmax(lo, hi)
  u  <- runif(length(y), min = lo, max = hi)
  qnorm(u)
}

# ── Full ZINB RQR ─────────────────────────────────────────────────
rqr_zinb_full <- function(y, mu, sigma, nu, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y_int <- as.integer(round(y))
  a     <- ifelse(y_int > 0,
                  pZINBI(y_int - 1L, mu = mu, sigma = sigma, nu = nu),
                  0)
  b     <- pZINBI(y_int, mu = mu, sigma = sigma, nu = nu)
  lo    <- pmax(pmin(a, 1 - 1e-8), 1e-8)
  hi    <- pmax(pmin(b, 1 - 1e-8), 1e-8)
  lo    <- pmin(lo, hi); hi <- pmax(lo, hi)
  u     <- runif(length(y), min = lo, max = hi)
  qnorm(u)
}

# ── Count RQR conditional on detection (shared by ZANBI and ZINB) ─
# For Y>0: F_cond(k) = (F_NBI(k) - F_NBI(0)) / (1 - F_NBI(0))
# For Y=0: returns NA
rqr_count_cond <- function(y, mu, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y_int  <- as.integer(round(y))
  p0     <- pNBI(0L, mu = mu, sigma = sigma)
  denom  <- pmax(1 - p0, 1e-8)
  cdf_y  <- (pNBI(y_int,      mu = mu, sigma = sigma) - p0) / denom
  cdf_y1 <- (pNBI(y_int - 1L, mu = mu, sigma = sigma) - p0) / denom
  cdf_y1 <- ifelse(y_int > 0, cdf_y1, 0)
  a  <- pmax(pmin(cdf_y1, 1 - 1e-8), 1e-8)
  b  <- pmax(pmin(cdf_y,  1 - 1e-8), 1e-8)
  lo <- pmin(a, b); hi <- pmax(a, b)
  u  <- runif(length(y), min = lo, max = hi)
  z  <- qnorm(u)
  z[y_int == 0] <- NA_real_
  z
}

# ── Internal: fit ZINB ────────────────────────────────────────────
.fit_zinb <- function(df_tr, mu_fml, sigma_fml, nu_fml, n_cyc) {
  tryCatch(
    gamlss(
      formula       = mu_fml,
      sigma.formula = sigma_fml,
      nu.formula    = nu_fml,
      family        = ZINBI(),
      data          = df_tr,
      control       = gamlss.control(n.cyc = n_cyc, trace = FALSE)
    ),
    error = function(e) e
  )
}

# ── ZINB main fitting function ────────────────────────────────────
# Conservative defaults: max_iter=2, max_remove_frac=0.05.
# Outlier detection uses full ZINB z-score on training data.

fit_zinb_gene <- function(y_train, y_test, X_train, X_test,
                            seed = NULL, n_cyc = 50,
                            outlier_z = 5.0, max_iter = 2L,
                            max_remove_frac = 0.05,
                            nu_type = "intercept") {
  # nu_type = "intercept" : nu ~ 1  (global constant; prevents overfitting)
  # nu_type = "full"      : nu ~ all covariates (same as mu/sigma)
  n_te      <- length(y_test)
  n_tr_orig <- length(y_train)
  safe_names <- sanitize_names(colnames(X_train))
  mu_rhs     <- paste(safe_names, collapse = " + ")

  df_tr <- as.data.frame(X_train)
  colnames(df_tr) <- safe_names
  df_tr$y__ <- as.integer(round(y_train))

  df_te <- as.data.frame(X_test)
  colnames(df_te) <- safe_names

  mu_fml    <- as.formula(paste("y__ ~",  mu_rhs))
  sigma_fml <- as.formula(paste("~",      mu_rhs))
  nu_fml    <- if (nu_type == "full") as.formula(paste("~", mu_rhs)) else as.formula("~ 1")

  na_result <- list(
    z_binary     = rep(NA_real_, n_te),
    z_count_cond = rep(NA_real_, n_te),
    z_full       = rep(NA_real_, n_te),
    mu_test = rep(NA_real_, n_te), sigma_test = rep(NA_real_, n_te),
    nu_test = rep(NA_real_, n_te),
    success = FALSE, msg = "", n_removed = 0L
  )

  keep      <- rep(TRUE, n_tr_orig)
  n_removed <- 0L

  for (iter in seq_len(max_iter)) {
    fit <- .fit_zinb(df_tr[keep, , drop = FALSE], mu_fml, sigma_fml, nu_fml, n_cyc)
    if (inherits(fit, "error")) {
      na_result$msg <- conditionMessage(fit)
      return(na_result)
    }

    pred_tr <- tryCatch(
      predictAll(fit, newdata = df_tr[keep, , drop = FALSE],
                 type = "response", data = df_tr[keep, , drop = FALSE]),
      error = function(e) e
    )
    if (inherits(pred_tr, "error")) break

    z_tr <- rqr_zinb_full(
      df_tr$y__[keep],
      mu    = as.numeric(pred_tr$mu),
      sigma = as.numeric(pred_tr$sigma),
      nu    = as.numeric(pred_tr$nu)
    )
    outlier <- is.finite(z_tr) & (abs(z_tr) > outlier_z)
    if (!any(outlier)) break
    if (sum(outlier) / n_tr_orig > max_remove_frac) break

    idx_keep  <- which(keep)
    keep[idx_keep[outlier]] <- FALSE
    n_removed <- n_removed + sum(outlier)
  }

  pred_te <- tryCatch(
    predictAll(fit, newdata = df_te, type = "response",
               data = df_tr[keep, , drop = FALSE]),
    error = function(e) e
  )
  if (inherits(pred_te, "error")) {
    na_result$msg <- conditionMessage(pred_te)
    return(na_result)
  }

  mu_te    <- as.numeric(pred_te$mu)
  sigma_te <- as.numeric(pred_te$sigma)
  nu_te    <- as.numeric(pred_te$nu)

  list(
    z_binary     = rqr_binary_zinb(y_test, mu = mu_te, sigma = sigma_te,
                                     nu = nu_te, seed = seed),
    z_count_cond = rqr_count_cond(y_test, mu = mu_te, sigma = sigma_te, seed = seed),
    z_full       = rqr_zinb_full(y_test, mu = mu_te, sigma = sigma_te,
                                   nu = nu_te, seed = seed),
    mu_test      = mu_te,
    sigma_test   = sigma_te,
    nu_test      = nu_te,
    success      = TRUE,
    msg          = "",
    n_removed    = n_removed
  )
}


# ══════════════════════════════════════════════════════════════════
# Model Engine helper — full-data NBI training for inference
# Returns coefficient vectors so Python can score without R.
# ══════════════════════════════════════════════════════════════════

train_nbi_coeffs <- function(y_train, X_train,
                              n_cyc = 50,
                              outlier_z = 5.0, max_iter = 2L,
                              max_remove_frac = 0.05) {
  n_tr_orig  <- length(y_train)
  safe_names <- sanitize_names(colnames(X_train))
  mu_rhs     <- paste(safe_names, collapse = " + ")
  mu_fml     <- as.formula(paste("y__ ~", mu_rhs))
  sigma_fml  <- as.formula(paste("~",     mu_rhs))

  df_tr <- as.data.frame(X_train)
  colnames(df_tr) <- safe_names
  df_tr$y__ <- as.integer(round(y_train))

  p          <- ncol(X_train) + 1L
  na_result  <- list(mu_coef = rep(NA_real_, p), sigma_coef = rep(NA_real_, p),
                     success = FALSE, msg = "", n_removed = 0L)

  keep <- rep(TRUE, n_tr_orig); n_removed <- 0L; fit <- NULL

  for (iter in seq_len(max_iter)) {
    fit <- .fit_nbi(df_tr[keep, , drop = FALSE], mu_fml, sigma_fml, n_cyc)
    if (inherits(fit, "error")) { na_result$msg <- conditionMessage(fit); return(na_result) }
    pred_tr <- tryCatch(
      predictAll(fit, newdata = df_tr[keep, , drop = FALSE],
                 type = "response", data = df_tr[keep, , drop = FALSE]),
      error = function(e) e)
    if (inherits(pred_tr, "error")) break
    z_tr    <- rqr_nbi(df_tr$y__[keep],
                        mu = as.numeric(pred_tr$mu), sigma = as.numeric(pred_tr$sigma))
    outlier <- is.finite(z_tr) & (abs(z_tr) > outlier_z)
    if (!any(outlier)) break
    if (sum(outlier) / n_tr_orig > max_remove_frac) break
    keep[which(keep)[outlier]] <- FALSE
    n_removed <- n_removed + sum(outlier)
  }

  if (is.null(fit) || inherits(fit, "error")) return(na_result)
  list(mu_coef = as.numeric(fit$mu.coefficients),
       sigma_coef = as.numeric(fit$sigma.coefficients),
       success = TRUE, msg = "", n_removed = n_removed)
}

# Train ZINBI on full HC data and return coefficient vectors.
# nu_type = "intercept" returns a 1-element nu_coef (the intercept).
train_zinbi_coeffs <- function(y_train, X_train,
                                n_cyc = 50,
                                outlier_z = 5.0, max_iter = 2L,
                                max_remove_frac = 0.05,
                                nu_type = "intercept") {
  n_tr_orig  <- length(y_train)
  safe_names <- sanitize_names(colnames(X_train))
  mu_rhs     <- paste(safe_names, collapse = " + ")

  df_tr <- as.data.frame(X_train)
  colnames(df_tr) <- safe_names
  df_tr$y__ <- as.integer(round(y_train))

  p         <- ncol(X_train) + 1L
  na_result <- list(mu_coef = rep(NA_real_, p), sigma_coef = rep(NA_real_, p),
                    nu_coef = NA_real_,
                    success = FALSE, msg = "", n_removed = 0L)

  mu_fml    <- as.formula(paste("y__ ~",  mu_rhs))
  sigma_fml <- as.formula(paste("~",      mu_rhs))
  nu_fml    <- if (nu_type == "full") as.formula(paste("~", mu_rhs)) else as.formula("~ 1")

  keep <- rep(TRUE, n_tr_orig); n_removed <- 0L; fit <- NULL

  for (iter in seq_len(max_iter)) {
    fit <- .fit_zinb(df_tr[keep, , drop = FALSE], mu_fml, sigma_fml, nu_fml, n_cyc)
    if (inherits(fit, "error")) { na_result$msg <- conditionMessage(fit); return(na_result) }
    pred_tr <- tryCatch(
      predictAll(fit, newdata = df_tr[keep, , drop = FALSE],
                 type = "response", data = df_tr[keep, , drop = FALSE]),
      error = function(e) e)
    if (inherits(pred_tr, "error")) break
    z_tr <- rqr_zinb_full(df_tr$y__[keep],
                            mu = as.numeric(pred_tr$mu),
                            sigma = as.numeric(pred_tr$sigma),
                            nu = as.numeric(pred_tr$nu))
    outlier <- is.finite(z_tr) & (abs(z_tr) > outlier_z)
    if (!any(outlier)) break
    if (sum(outlier) / n_tr_orig > max_remove_frac) break
    keep[which(keep)[outlier]] <- FALSE
    n_removed <- n_removed + sum(outlier)
  }

  if (is.null(fit) || inherits(fit, "error")) return(na_result)
  list(mu_coef    = as.numeric(fit$mu.coefficients),
       sigma_coef = as.numeric(fit$sigma.coefficients),
       nu_coef    = as.numeric(fit$nu.coefficients),
       success    = TRUE, msg = "", n_removed = n_removed)
}


# ══════════════════════════════════════════════════════════════════
# Warm-start NBI fitting
#
# Uses the predictions from a previously fitted model (mu_start,
# sigma_start) as starting values for the GAMLSS inner algorithm,
# instead of the default link-function mean. This reduces the number
# of backfitting cycles needed and preserves learned structure when
# only a small amount of new data is added.
#
# mu_start, sigma_start : numeric vectors (length = nrow(X_train))
#   Predicted fitted values from the old model on the TRAINING data,
#   i.e. exp(X_train_aug %*% mu_coef_old) etc.
# ══════════════════════════════════════════════════════════════════

train_nbi_coeffs_warm <- function(y_train, X_train,
                                   mu_start, sigma_start,
                                   n_cyc = 50,
                                   outlier_z = 5.0, max_iter = 2L,
                                   max_remove_frac = 0.10) {
  n_tr_orig  <- length(y_train)
  safe_names <- sanitize_names(colnames(X_train))
  mu_rhs     <- paste(safe_names, collapse = " + ")
  mu_fml     <- as.formula(paste("y__ ~", mu_rhs))
  sigma_fml  <- as.formula(paste("~",     mu_rhs))

  df_tr <- as.data.frame(X_train)
  colnames(df_tr) <- safe_names
  df_tr$y__ <- as.integer(round(y_train))

  p         <- ncol(X_train) + 1L
  na_result <- list(mu_coef = rep(NA_real_, p), sigma_coef = rep(NA_real_, p),
                    success = FALSE, msg = "", n_removed = 0L)

  keep      <- rep(TRUE, n_tr_orig)
  n_removed <- 0L
  fit       <- NULL
  mu_s      <- as.numeric(mu_start)
  sigma_s   <- as.numeric(sigma_start)

  for (iter in seq_len(max_iter)) {
    df_sub    <- df_tr[keep, , drop = FALSE]
    mu_sub    <- mu_s[keep]
    sigma_sub <- sigma_s[keep]

    fit <- tryCatch(
      gamlss(
        formula       = mu_fml,
        sigma.formula = sigma_fml,
        family        = NBI(),
        data          = df_sub,
        mu.start      = mu_sub,
        sigma.start   = sigma_sub,
        control       = gamlss.control(n.cyc = n_cyc, trace = FALSE)
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      na_result$msg <- conditionMessage(fit); return(na_result)
    }

    pred_tr <- tryCatch(
      predictAll(fit, newdata = df_sub, type = "response", data = df_sub),
      error = function(e) e
    )
    if (inherits(pred_tr, "error")) break

    z_tr    <- rqr_nbi(df_sub$y__,
                        mu = as.numeric(pred_tr$mu),
                        sigma = as.numeric(pred_tr$sigma))
    outlier <- is.finite(z_tr) & (abs(z_tr) > outlier_z)
    if (!any(outlier)) break
    if (sum(outlier) / n_tr_orig > max_remove_frac) break

    keep[which(keep)[outlier]] <- FALSE
    n_removed <- n_removed + sum(outlier)

    # Update starting values for the surviving samples only (next iter)
    mu_s    <- as.numeric(pred_tr$mu)
    sigma_s <- as.numeric(pred_tr$sigma)
  }

  if (is.null(fit) || inherits(fit, "error")) return(na_result)
  list(mu_coef    = as.numeric(fit$mu.coefficients),
       sigma_coef = as.numeric(fit$sigma.coefficients),
       success    = TRUE, msg = "", n_removed = n_removed)
}
