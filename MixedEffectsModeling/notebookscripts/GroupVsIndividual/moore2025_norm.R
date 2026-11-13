args = commandArgs(trailingOnly = TRUE)
repo = normalizePath(args[1])
out = args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out = normalizePath(out)
setwd(repo)
source("normalization.R")

rna_counts = read.csv("pdac_genecount.csv.gz", row.names = 1)
intrinsic_genes = read.csv("updated_intrinsic.csv")
intrinsic_genes = intrinsic_genes[!is.na(intrinsic_genes$Intrinsic), ]
write.csv(tmm(cf_norm(rna_counts, intrinsic_genes)), gzfile(file.path(out, "cfnorm.csv.gz")))
write.csv(cpm(rna_counts), gzfile(file.path(out, "cpm.csv.gz")))
write.csv(tpm(rna_counts), gzfile(file.path(out, "tpm.csv.gz")))
