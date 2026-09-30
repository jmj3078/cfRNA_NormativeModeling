suppressPackageStartupMessages({library(glmmTMB); library(parallel); library(jsonlite)})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", grep("--file=", .args, value = TRUE))
source(file.path(dirname(normalizePath(.script_path)), "glmm_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
RES <- args[1]; WORK <- args[2]; OUT <- args[3]
NREP <- as.integer(args[4]); CORES <- as.integer(args[5]); PRIOR <- args[6]

X <- as.matrix(read.csv(file.path(WORK, "X.csv.gz"), row.names = 1))
batch <- read.csv(file.path(WORK, "batch.csv.gz"), row.names = 1)[[1]]
sn <- sanitize_names(colnames(X)); colnames(X) <- sn
Xa <- cbind(1, X); n <- nrow(Xa); p <- ncol(Xa)
pri <- as.numeric(fromJSON(PRIOR)$tau_slope)
res <- read.csv(RES); res <- res[toupper(as.character(res$ok)) == "TRUE", ]
mu_c <- paste0("mu_coef_", 0:10); dp_c <- paste0("disp_coef_", 0:10)
lev <- unique(batch); bidx <- match(batch, lev)
TOPK <- 50

# Full PCIS vector + effective df, mirroring pcis_outliers() in glmm_helpers.R.
pcis_vec <- function(fit, y, trend_alpha) {
  mu <- tryCatch(as.numeric(predict(fit, type = "response")), error = function(e) NULL)
  if (is.null(mu) || !all(is.finite(mu))) return(NULL)
  r <- (y - mu) / sqrt(mu + trend_alpha * mu^2); sw <- sqrt(mu / (1 + trend_alpha * mu))
  tau2 <- tryCatch(as.numeric(VarCorr(fit)$cond[[1]][1, 1]), error = function(e) NA_real_)
  Z <- model.matrix(~ 0 + factor(batch)); h <- NULL
  if (isTRUE(is.finite(tau2))) {
    Mw <- cbind(Xa, Z) * sw
    A <- tryCatch(solve(crossprod(Mw) + diag(c(rep(0, p), rep(1 / max(tau2, 1e-6), ncol(Z))))), error = function(e) NULL)
    if (!is.null(A)) h <- rowSums((Mw %*% A) * Mw)
  }
  if (is.null(h)) {
    Xw <- Xa * sw; A <- tryCatch(solve(crossprod(Xw)), error = function(e) NULL)
    if (is.null(A)) return(NULL); h <- rowSums((Xw %*% A) * Xw)
  }
  h <- pmin(pmax(h, 0), 1 - 1e-8); pe <- sum(h)
  if (!is.finite(pe) || pe < 1 || n - pe < 4) return(NULL)
  list(pcis = r^2 * h / (pe * (1 - h)^2), p_eff = pe)
}

work <- function(i) {
  row <- res[i, ]; st <- row$stage; ta <- row$trend_alpha
  if (!is.finite(ta)) return(NULL)
  bc <- as.numeric(row[mu_c]); gc_ <- as.numeric(row[dp_c])
  if (any(!is.finite(bc)) || !is.finite(gc_[1])) return(NULL)
  mu0 <- exp(pmin(pmax(as.numeric(Xa %*% bc), -30), 30))
  alpha <- if (st == "nbi_full_eb") exp(-pmin(pmax(as.numeric(Xa %*% ifelse(is.finite(gc_), gc_, 0)), -30), 30))
           else rep(exp(-gc_[1]), n)
  t2 <- if (is.finite(row$tau2)) max(row$tau2, 0) else 0
  disp_fml <- if (st == "nbi_full_eb") as.formula(paste("~", paste(sn, collapse = "+"))) else as.formula("~ 1")
  prd <- if (st == "nbi_full_eb") data.frame(prior = sprintf("normal(0, %.6g)", pri), class = "betad", coef = sn) else NULL
  out <- list()
  for (b in seq_len(NREP)) {
    u <- if (t2 > 0) rnorm(length(lev), 0, sqrt(t2))[bidx] else rep(0, n)
    ys <- rnbinom(n, size = pmin(pmax(1 / alpha, 1e-6), 1e8), mu = mu0 * exp(u))
    d <- as.data.frame(X); d$y__ <- as.integer(ys); d$batch__ <- factor(batch)
    f <- tryCatch({
      fml <- as.formula(paste("y__ ~", paste(sn, collapse = "+"), "+ (1|batch__)"))
      if (is.null(prd)) glmmTMB(fml, dispformula = disp_fml, family = nbinom2(), data = d)
      else glmmTMB(fml, dispformula = disp_fml, family = nbinom2(), data = d, priors = prd)
    }, error = function(e) NULL)
    if (is.null(f)) next
    pv <- pcis_vec(f, d$y__, ta)
    if (is.null(pv)) next
    top <- sort(pv$pcis, decreasing = TRUE)[seq_len(min(TOPK, n))]
    out[[length(out) + 1]] <- data.frame(gene = row$gene, rep = b, stage = st, n_obs = n,
        log_mu = log(mean(mu0)), trend_alpha = ta, tau2 = t2, p_eff = pv$p_eff,
        rank = seq_along(top), pcis = top)
  }
  if (!length(out)) NULL else do.call(rbind, out)
}

set.seed(11)
r <- mclapply(seq_len(nrow(res)), work, mc.cores = CORES)
df <- do.call(rbind, r[!sapply(r, is.null)])
write.csv(df, OUT, row.names = FALSE)
cat(sprintf("null sim done: genes=%d reps=%d rows=%d -> %s\n",
            length(unique(df$gene)), NREP, nrow(df), OUT))
