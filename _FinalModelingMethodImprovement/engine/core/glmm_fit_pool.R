suppressPackageStartupMessages({
  library(optparse); library(jsonlite)
})

.args <- commandArgs(trailingOnly = FALSE)
.script_path <- sub("--file=", "", grep("--file=", .args, value = TRUE))
source(file.path(dirname(normalizePath(.script_path)), "glmm_helpers.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--x", type = "character"), make_option("--y", type = "character"),
  make_option("--batch", type = "character"), make_option("--genes", type = "character"),
  make_option("--rare-overdisp-thr", type = "double", default = 2.0),
  make_option("--out", type = "character")
)))

X <- as.matrix(read.csv(opt$x, row.names = 1))
Y <- as.matrix(read.csv(opt$y, row.names = 1))
batch <- read.csv(opt$batch, row.names = 1)[[1]]
genes <- read.csv(opt$genes)$gene
n_hc <- nrow(X)
eps <- 1 / (2 * n_hc)
mean_hc <- colMeans(Y)

fit <- fit_pooled_glmm(Y, X, batch, mean_hc, eps, opt$`rare-overdisp-thr`)

out <- list(ok = isTRUE(fit$ok), family = if (isTRUE(fit$ok)) fit$family else NA,
           beta = if (isTRUE(fit$ok)) as.numeric(fit$beta) else numeric(0),
           alpha = if (isTRUE(fit$ok) && !is.na(fit$alpha)) as.numeric(fit$alpha) else NA,
           tau2 = if (isTRUE(fit$ok)) as.numeric(fit$tau2) else NA,
           overdisp_ratio = if (isTRUE(fit$ok)) as.numeric(fit$overdisp_ratio) else NA,
           mult_lo = if (isTRUE(fit$ok)) as.numeric(fit$mult_lo) else NA,
           mult_hi = if (isTRUE(fit$ok)) as.numeric(fit$mult_hi) else NA,
           eps = eps, gene = genes, mean_hc = as.numeric(mean_hc))
write(toJSON(out, auto_unbox = TRUE, na = "null"), opt$out)
cat("DONE\n")
