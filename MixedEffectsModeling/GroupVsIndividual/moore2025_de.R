args = commandArgs(trailingOnly = TRUE)
repo = normalizePath(args[1])
out = args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
out = normalizePath(out)
setwd(repo)
source("normalization.R")
source("DE_wilcox.R")

rna_counts = read.csv("pdac_genecount.csv.gz", row.names = 1)
meta = read.csv("pdac_meta.csv")
meta$Group[meta$Group == "Islet Cell Tumor"] = "Other Cancer"
meta$Group[meta$Group == "Cancer_other"] = "Other Cancer"
meta$Group[meta$Group == "Benign Pancreas"] = "Benign"
ids_train = meta$Cohort == "CEDAR"
ids_val = meta$Cohort == "BCC"

intrinsic_genes = read.csv("updated_intrinsic.csv")
intrinsic_genes = intrinsic_genes[!is.na(intrinsic_genes$Intrinsic), ]
rna_counts_norm = tmm(cf_norm(rna_counts, intrinsic_genes))

train_counts_norm = rna_counts_norm[, ids_train]
val_counts_norm = rna_counts_norm[, ids_val]
V = apply(train_counts_norm, 1, var)
train_counts_norm = train_counts_norm[V > 1, ]
val_counts_norm = val_counts_norm[V > 1, ]

DE = DE_wilcox(train_counts_norm, meta$Group[ids_train] == "PDAC")
DE_val = DE_wilcox(val_counts_norm, meta$Group[ids_val] == "PDAC")
cat("CEDAR genes tested:", nrow(DE), "| FDR<0.05:", sum(as.numeric(DE$FDR) < 0.05), "\n")
cat(sort(DE$genes[as.numeric(DE$FDR) < 0.05]), "\n")
write.csv(DE, file.path(out, "moore2025_de_cedar.csv"), row.names = FALSE)
write.csv(DE_val, file.path(out, "moore2025_de_bcc.csv"), row.names = FALSE)
