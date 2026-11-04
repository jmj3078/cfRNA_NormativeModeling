suppressPackageStartupMessages({library(glmmTMB); library(parallel); library(jsonlite)})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", grep("--file=", .args, value = TRUE))
source(file.path(dirname(normalizePath(.script_path)), "..", "..", "core", "glmm_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
WORK <- args[1]; CALIB <- args[2]; GENES <- args[3]; PRIOR <- args[4]; FP_PATH <- args[5]; OUT <- args[6]
NREP <- as.integer(args[7])

X <- as.matrix(read.csv(file.path(WORK, "X.csv.gz"), row.names = 1))
batch <- read.csv(file.path(WORK, "batch.csv.gz"), row.names = 1)[[1]]
fp <- fromJSON(FP_PATH)
tau_slope <- as.numeric(fromJSON(PRIOR)$tau_slope)
sn <- sanitize_names(colnames(X)); colnames(X) <- sn
Xa <- cbind(1, X); n <- nrow(Xa); p <- ncol(Xa)
Z <- model.matrix(~ 0 + factor(batch))
lev <- unique(batch); bidx <- match(batch, lev)
calib <- read.csv(CALIB)
calib <- calib[calib$gene %in% read.csv(GENES)$gene, ]
mu_c <- paste0("mu_coef_", 0:(p - 1)); dp_c <- paste0("disp_coef_", 0:(p - 1))
mu_fml <- as.formula(paste("y__ ~", paste(sn, collapse = " + "), "+ (1 | batch__)"))
disp_fml <- as.formula(paste("~", paste(sn, collapse = " + ")))
prior_df <- data.frame(prior = sprintf("normal(0, %.6g)", tau_slope), class = "betad", coef = sn)
fit_d <- function(d) tryCatch(glmmTMB(mu_fml, dispformula = disp_fml, family = nbinom2(), data = d, priors = prior_df),
                              error = function(e) NULL)

scores <- function(fit, y, ta) {
  mu <- as.numeric(predict(fit, type = "response"))
  a_own <- exp(-as.numeric(Xa %*% fixef(fit)$disp))
  tau2 <- as.numeric(VarCorr(fit)$cond[[1]][1, 1])
  mixed_h <- function(alpha) {
    Mw <- cbind(Xa, Z) * sqrt(mu / (1 + alpha * mu))
    A <- solve(crossprod(Mw) + diag(c(rep(0, p), rep(1 / max(tau2, 1e-6), ncol(Z)))))
    pmin(pmax(rowSums((Mw %*% A) * Mw), 0), 1 - 1e-8)
  }
  r_t <- (y - mu) / sqrt(mu + ta * mu^2); h_t <- mixed_h(ta)
  r_o <- (y - mu) / sqrt(mu + a_own * mu^2); h_o <- mixed_h(a_own)
  Xw <- Xa * sqrt(mu / (1 + a_own * mu))
  h_f <- pmin(pmax(rowSums((Xw %*% solve(crossprod(Xw))) * Xw), 0), 1 - 1e-8)
  list(pcis_trend = r_t^2 * h_t / (sum(h_t) * (1 - h_t)^2),
       pcis_own = r_o^2 * h_o / (sum(h_o) * (1 - h_o)^2),
       cooks = r_o^2 * h_f / (p * (1 - h_f)^2),
       alpha_own = exp(mean(log(a_own))))
}

conds <- expand.grid(k = c(0L, 2L, 3L), m = c(5, 10, 20), place = c("random", "same_batch"),
                     stringsAsFactors = FALSE)
conds <- conds[!(conds$k == 0 & (conds$m != 5 | conds$place != "random")), ]
tasks <- expand.grid(g = seq_len(nrow(calib)), c = seq_len(nrow(conds)), rep = seq_len(NREP))

work <- function(t) {
  row <- calib[tasks$g[t], ]; cd <- conds[tasks$c[t], ]; b <- tasks$rep[t]
  set.seed(1000 * tasks$g[t] + 100 * tasks$c[t] + b)
  ta <- row$trend_alpha
  mu0 <- exp(pmin(pmax(as.numeric(Xa %*% as.numeric(row[mu_c])), -30), 30))
  a_true <- exp(-pmin(pmax(as.numeric(Xa %*% as.numeric(row[dp_c])), -30), 30))
  u <- rnorm(length(lev), 0, sqrt(max(row$tau2, 0)))[bidx]
  mu_t <- mu0 * exp(u)
  y <- rnbinom(n, size = 1 / a_true, mu = mu_t)
  inj <- integer(0)
  if (cd$k > 0) {
    pool <- if (cd$place == "same_batch") which(batch == sample(names(which(table(batch) >= 10)), 1)) else seq_len(n)
    inj <- sample(pool, cd$k)
    y[inj] <- rnbinom(cd$k, size = 1 / a_true[inj], mu = cd$m * mu_t[inj])
  }
  d <- as.data.frame(X); d$y__ <- as.integer(y); d$batch__ <- factor(batch)
  fit <- fit_d(d)
  if (is.null(fit)) return(NULL)
  s <- scores(fit, y, ta)
  drop <- pcis_outliers(fit, Xa, y, batch, ta, fp$pcis_cut, fp$max_outlier_frac)
  a_after <- NA_real_
  if (length(drop)) {
    f2 <- fit_d(droplevels(d[-drop, , drop = FALSE]))
    if (!is.null(f2)) a_after <- exp(mean(log(exp(-as.numeric(Xa %*% fixef(f2)$disp)))))
  }
  clean <- setdiff(seq_len(n), inj)
  hit <- function(v, cut) if (length(inj)) mean(v[inj] > cut) else NA_real_
  fp_n <- function(v, cut) sum(v[clean] > cut)
  data.frame(gene = row$gene, k = cd$k, m = cd$m, place = cd$place, rep = b, trend_alpha = ta,
    alpha_true = exp(mean(log(a_true))), alpha_own = s$alpha_own, alpha_after_removal = a_after,
    n_removed = length(drop), removed_injected = sum(inj %in% drop),
    tpr_pcis_trend = hit(s$pcis_trend, fp$pcis_cut), tpr_pcis_own = hit(s$pcis_own, fp$pcis_cut),
    tpr_cooks_1 = hit(s$cooks, 1), tpr_cooks_4n = hit(s$cooks, 4 / n),
    fp_pcis_trend = fp_n(s$pcis_trend, fp$pcis_cut), fp_pcis_own = fp_n(s$pcis_own, fp$pcis_cut),
    fp_cooks_1 = fp_n(s$cooks, 1), fp_cooks_4n = fp_n(s$cooks, 4 / n),
    min_pcis_trend_inj = if (length(inj)) min(s$pcis_trend[inj]) else NA_real_,
    min_pcis_own_inj = if (length(inj)) min(s$pcis_own[inj]) else NA_real_,
    min_cooks_inj = if (length(inj)) min(s$cooks[inj]) else NA_real_)
}

res <- mclapply(seq_len(nrow(tasks)), work, mc.cores = as.integer(fp$cores))
res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
write.csv(res, OUT, row.names = FALSE)
cat(sprintf("sim done: %d rows (%d genes x %d conditions x %d reps)\n", nrow(res), nrow(calib), nrow(conds), NREP))
