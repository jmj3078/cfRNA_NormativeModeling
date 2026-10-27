# PEER (Stegle et al. 2012) factor inference for one fold's TRAIN matrix.
#
# Runs in peer_env (R 3.4.1 + r-peer 1.3 + libgfortran3) -- r-peer is pinned to R 3.4 and
# cannot coexist with OUTRIDER's R 4.3, so this is a file-mediated step: Python writes the
# standardized train matrix, this script writes back the frozen factor weights W and the
# per-gene residual variances, and Python does the held-out projection.
#
# Only W and the residual variances cross the boundary, because those are the frozen
# quantities; the held-out factor values are a posterior mean computed against them, which
# is the inductive step PEER's own API does not expose.
#
# Usage: Rscript fit_peer.R <train_csv> <out_prefix> <K> <max_iter>
suppressMessages(library(peer))

args <- commandArgs(trailingOnly = TRUE)
train_csv <- args[1]; out_prefix <- args[2]
K <- as.integer(args[3]); max_iter <- as.integer(args[4])

Y <- as.matrix(read.csv(train_csv, row.names = 1, check.names = FALSE))  # samples x genes
cat("train matrix:", dim(Y), "\n")

m <- PEER()
PEER_setPhenoMean(m, Y)
PEER_setNk(m, K)
PEER_setNmax_iterations(m, max_iter)
PEER_setAdd_mean(m, FALSE)   # Python already centered each gene
PEER_update(m)

W <- PEER_getW(m)                  # genes x K -- the frozen basis
eps <- PEER_getEps(m)              # per-gene noise PRECISION (1/variance)
alpha <- PEER_getAlpha(m)          # ARD precision per factor; large => pruned
# PEER_getResidualVars is not per-gene (length 3K), so Eps is the right frozen quantity.
stopifnot(length(as.vector(eps)) == ncol(Y))

write.csv(W, paste0(out_prefix, "_W.csv"), row.names = FALSE)
write.csv(data.frame(eps = as.vector(eps)), paste0(out_prefix, "_eps.csv"), row.names = FALSE)
write.csv(data.frame(alpha = as.vector(alpha)), paste0(out_prefix, "_alpha.csv"), row.names = FALSE)
cat("W:", dim(W), "| effective factors (alpha < 10x min):",
    sum(as.vector(alpha) < 10 * min(as.vector(alpha))), "of", K, "\n")
cat("PEER_FIT_DONE\n")
