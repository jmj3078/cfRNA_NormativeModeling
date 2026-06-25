# 1. 라이브러리 로드
library(RUVSeq)
library(EDASeq)
library(edgeR)
library(dplyr)
library(tidyr)
library(limma)
library(readxl)
library(tibble)
# -----------------------------------------------------------------------------
# STEP 1: 데이터 로딩 및 병합
# -----------------------------------------------------------------------------
meta     <- read_excel("/project/cfRNA_Disentaglement/OpenAccess_nfcore/Meta/Meta_Processed.xlsx") %>% column_to_rownames(var = colnames(.)[1])
annot    <- read.table("/project/cfRNA_Disentaglement/Data/GECODEv49_Annot.tsv", header=TRUE, row.names=1, sep='\t')
palangodb <- read.table("/project/cfRNA_Disentaglement/Data/PalangoDB_CellTypeMarkers.tsv", header=TRUE, sep='\t')
counts <- read.table("/project/cfRNA_Disentaglement/OpenAccess_nfcore/Merged_filtered_salmon_quant_gene.tsv", header=TRUE, row.names=1, sep='\t')
counts[is.na(counts)] <- 0
lib_sizes <- colSums(counts)
zero_samples <- names(lib_sizes[lib_sizes == 0])
counts <- counts[, order(colnames(counts))]

if (length(zero_samples) > 0) {
  print(zero_samples)
  counts <- counts[, lib_sizes > 0]
  meta <- meta[colnames(counts), ]
  cat("제거 완료. 남은 샘플 수:", ncol(counts), "\n")
} else {
  cat("모든 샘플의 Library Size가 정상(>0)입니다.\n")
}

print(dim(counts))
print(dim(meta))
# -----------------------------------------------------------------------------
# STEP 2: 전처리 (Metadata & Gene Filtering)
# -----------------------------------------------------------------------------
# 2.2 샘플 동기화
common_samples <- intersect(colnames(counts), rownames(meta))
counts_common <- counts[, common_samples]
meta_common   <- meta[common_samples, ]
dim(counts_common)
dim(meta_common)
stopifnot(all(colnames(counts_common) == rownames(meta_common)))

# 2.3 Protein Coding 유전자 필터링
pc_genes  <- rownames(annot)[annot$GeneType == "protein_coding"]
counts_pc <- counts_common[rownames(counts_common) %in% pc_genes, ]
annot_pc  <- annot[rownames(counts_pc), ] # 순서 동기화
valid_len_gc <- !is.na(annot_pc$Length) & annot_pc$Length > 0 & !is.na(annot_pc$GC_Percent)
counts_pc    <- counts_pc[valid_len_gc, ]
annot_pc     <- annot_pc[valid_len_gc, ]
dim(counts_pc)

# 2.4 Platelet Control Genes 추출
platelet_syms <- palangodb[palangodb$cell.type == "Platelets", "official.gene.symbol"]
platelet_ids  <- rownames(annot_pc)[annot_pc$GeneName %in% platelet_syms]
# -----------------------------------------------------------------------------
# STEP 3: 정규화 수행
# -----------------------------------------------------------------------------
out_base_dir <- "/project/cfRNA_Disentaglement/OpenAccess_nfcore/Processed/"
if(!dir.exists(out_base_dir)) dir.create(out_base_dir, recursive = TRUE)
meta$Group <- paste(meta$Author, meta$tissue, sep = "_")
group_list <- unique(meta$Group)
print(group_list)
# Group을 기준으로 루프 실행
clean_matrix <- function(mat) {
    mat[is.na(mat)] <- 0
    mat[is.infinite(mat)] <- 0
    return(mat)
}

group_list <- unique(meta$Group)
message(paste("Start processing for", length(group_list), "groups:", paste(group_list, collapse=", ")))

for (group in group_list) {
    message(paste0("\n>>> Processing Group: ", group, " <<<"))
    
    # 수정 1: meta_common -> meta 통일
    target_samples <- rownames(meta)[meta$Group == group]
    group_samples <- intersect(target_samples, colnames(counts_pc))
    if (length(group_samples) == 0) {
        message(paste("  [SKIP] No matching samples found in counts_pc for group:", group))
        next
    }
    
    sub_counts <- counts_pc[, group_samples, drop=FALSE]
    # 0 라이브러리 샘플 필터링
    valid_samples <- colSums(sub_counts) > 0
    if (!all(valid_samples)) {
        sub_counts <- sub_counts[, valid_samples, drop=FALSE]
        if(ncol(sub_counts) == 0) {
            message("  [SKIP] No valid samples left in this group.")
            next
        }
    }
    
    # Low Expression Gene Filtering
    keep_genes <- rowSums(sub_counts) > 0
    sub_counts <- sub_counts[keep_genes, ]
    sub_annot  <- annot_pc[keep_genes, ]
    
    control_genes <- intersect(platelet_ids, rownames(sub_counts))
    message(paste("   Samples:", ncol(sub_counts), "| Genes:", nrow(sub_counts), "| Control Genes:", length(control_genes)))
    
    results_list <- list()
    results_list[["Raw"]] <- sub_counts
    
    # ---------------------------------------------------------
    # 3.1 Basic Normalization
    # ---------------------------------------------------------
    dge <- DGEList(counts = as.matrix(sub_counts))
    dge <- calcNormFactors(dge, method = "TMM")
    results_list[["TMM_log2"]] <- clean_matrix(cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))
    
    gene_len_kb <- sub_annot$Length / 1000
    rpk <- sub_counts / gene_len_kb
    scale_factor <- colSums(rpk)
    scale_factor[scale_factor == 0] <- 1
    tpm_val <- t(t(rpk) * 1e6 / scale_factor)
    results_list[["TPM_log2"]] <- clean_matrix(log2(tpm_val + 1))
    
    results_list[["FPKM_log2"]] <- clean_matrix(rpkm(as.matrix(sub_counts), 
                                                     gene.length = sub_annot$Length, 
                                                     log = TRUE, 
                                                     prior.count = 1))
    
    # ---------------------------------------------------------
    # 3.2 EDASeq 
    # ---------------------------------------------------------
    eda_set <- newSeqExpressionSet(as.matrix(sub_counts),
                                   featureData = data.frame(
                                       gc = sub_annot$GC_Percent, 
                                       length = sub_annot$Length,
                                       row.names = rownames(sub_counts)
                                   ))
    
    set_within <- withinLaneNormalization(withinLaneNormalization(eda_set, "gc", which="full"), "length", which="full")
    norm_eda   <- betweenLaneNormalization(set_within, which="upper")
    results_list[["EDA_Full_All"]] <- clean_matrix(log2(normCounts(norm_eda) + 1))
    
    # ---------------------------------------------------------
    # 3.3 RUVg
    # ---------------------------------------------------------
    input_ruv <- results_list[["TMM_log2"]]
    if (length(control_genes) < 5) {
        message("   [WARNING] Too few control genes found. Skipping RUV.")
    } else {
        for (k in c(1, 2, 3)) {
            set_ruvg <- RUVg(input_ruv, control_genes, k = k, isLog = TRUE, center = TRUE)
            W_platelet <- set_ruvg$W
            
            # 수정 2: cohort -> group 변수명 동기화
            write.csv(W_platelet, file.path(out_base_dir, paste0(group, "_W_factor_k", k, ".csv")))
            
            key_ruvg <- paste0("RUVg_Platelet_k", k)
            results_list[[key_ruvg]] <- clean_matrix(set_ruvg$normalizedCounts)
            
            key_full <- paste0("Proposed_Full_k", k)
            results_list[[key_full]] <- clean_matrix(
                removeBatchEffect(results_list[["EDA_Full_All"]], covariates = W_platelet)
            )
        }
    }
    
    for(name in names(results_list)) {
        file_path <- paste0(out_base_dir, "Norm_", group, "_", name, ".csv")
        write.csv(as.data.frame(results_list[[name]]), file = file_path, row.names = TRUE)
    }
    message(paste("   Saved files for:", group))
}

# -----------------------------------------------------------------------------
# STEP 4: Summary 병합 및 Metadata 업데이트
# -----------------------------------------------------------------------------
# summary_files <- list.files(path = base_path, 
#                             pattern = "_featureCounts_summary.tsv$", 
#                             recursive = TRUE, 
#                             full.names = TRUE)

# print(paste("Found summary files:", length(summary_files)))
# summary_combined <- NULL

# for (i in 1:length(summary_files)) {
#     file_path <- summary_files[i]
#     current_sum <- read.delim(file_path, header=TRUE, sep="\t", 
#                               check.names=FALSE, stringsAsFactors=FALSE)
#     if (is.null(summary_combined)) {
#         summary_combined <- current_sum
#     } else {
#         summary_combined <- merge(summary_combined, current_sum, by = 1, all = TRUE)
#     }
# }

# summary_combined[is.na(summary_combined)] <- 0
# rownames(summary_combined) <- summary_combined[, 1]
# summary_combined <- summary_combined[, -1]
# summary_t <- as.data.frame(t(summary_combined))

# # 수정 4: row.names 소실 방지 파이프라인 적용
# library(tibble)
# meta <- merge(meta, summary_t, by = "row.names", all.x = TRUE) %>%
#     column_to_rownames("Row.names")

head(meta)
dim(meta)
write.csv(meta, paste0(out_base_dir, "Meta_Processed.csv"))
write.csv(annot_pc, paste0(out_base_dir, "Annot_PC.csv")) # 저장 객체를 전체 annot에서 실제로 사용된 annot_pc로 변경 권장

