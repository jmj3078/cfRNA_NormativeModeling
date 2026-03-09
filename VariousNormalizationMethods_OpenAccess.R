# 1. 라이브러리 로드
library(RUVSeq)
library(EDASeq)
library(edgeR)
library(dplyr)
library(tidyr)
library(limma)
# -----------------------------------------------------------------------------
# STEP 1: 데이터 로딩 및 병합
# -----------------------------------------------------------------------------
meta     <- read.table("/project/cfRNA_Disentaglement/OpenAccess_nfcore/Meta/OpenAccessMeta_Processed.tsv", header=TRUE, row.names=1, sep='\t')
annot    <- read.table("/project/cfRNA_Disentaglement/Data/GECODEv49_Annot.tsv", header=TRUE, row.names=1, sep='\t')
palangodb <- read.table("/project/cfRNA_Disentaglement/Data/PalangoDB_CellTypeMarkers.tsv", header=TRUE, sep='\t')

library(fs) 
base_path <- "/project/cfRNA_Disentaglement/Data/OpenAccess/"
target_files <- list.files(path = base_path, 
                           pattern = "_featureCounts_matrix.tsv$", 
                           recursive = TRUE, 
                           full.names = TRUE)
counts <- NULL

for (i in 1:length(target_files)) {
  file_path <- target_files[i]
  cat(sprintf("[%d/%d] Reading : %s\n", i, length(target_files), basename(file_path)))
  current_data <- read.delim(file_path, 
                             header = TRUE, 
                             sep = "\t", 
                             check.names = FALSE, 
                             stringsAsFactors = FALSE)
  
  if (is.null(counts)) {
    counts <- current_data
  } else {
    counts <- merge(counts, current_data, by = 1, all = TRUE)
  }
}

rownames(counts) <- counts[, 1]
counts <- counts[, -1]
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

target_indices <- meta$tissue %in% c("Plasma", "Serum")
meta_filtered <- meta[target_indices, ]
counts_filtered <- counts[, rownames(meta_filtered)]
dim(meta_filtered)
dim(counts_filtered)
write.csv(counts_filtered, "/project/cfRNA_Disentaglement/Data/OpenAccess/Processed/Merged_featureCounts.csv")


# -----------------------------------------------------------------------------
# STEP 2: 전처리 (Metadata & Gene Filtering)
# -----------------------------------------------------------------------------
# 2.2 샘플 동기화
common_samples <- intersect(colnames(counts), rownames(meta))
dim(counts)
dim(meta)

counts <- counts[, common_samples]
meta   <- meta[common_samples, ]
stopifnot(all(colnames(counts) == rownames(meta)))
# 2.3 Protein Coding 유전자 필터링
pc_genes  <- rownames(annot)[annot$GeneType == "protein_coding"]

counts_pc <- counts[rownames(counts) %in% pc_genes, ]
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
out_base_dir <- "/project/cfRNA_Disentaglement/Data/OpenAccess/Processed/"
if(!dir.exists(out_base_dir)) dir.create(out_base_dir, recursive = TRUE)
cohort_list <- unique(meta$BioProject)
cohort_list
clean_matrix <- function(mat) {
    mat[is.na(mat)] <- 0
    mat[is.infinite(mat)] <- 0
    return(mat)
}

message(paste("Start processing for", length(cohort_list), "cohorts:", paste(cohort_list, collapse=", ")))

for (cohort in cohort_list) {

    message(paste0("\n>>> Processing Cohort: ", cohort, " <<<"))
    
    # 2.1 Subset Data (해당 코호트 샘플만 추출)
    cohort_samples <- rownames(meta)[meta$BioProject == cohort]
    sub_counts <- counts_pc[, cohort_samples, drop=FALSE]
    
    # 2.2 Low Expression Gene Filtering (이 코호트 내에서 발현 안되는 유전자 제거)
    # 중요: TMM/EDASeq은 0이 너무 많으면 에러가 나거나 왜곡됨
    keep_genes <- rowSums(sub_counts) > 0
    sub_counts <- sub_counts[keep_genes, ]
    sub_annot  <- annot_pc[keep_genes, ]
    
    # Platelet Gene도 현재 살아남은 유전자와 교집합
    control_genes <- intersect(platelet_ids, rownames(sub_counts))
    
    message(paste("   Samples:", ncol(sub_counts), "| Genes:", nrow(sub_counts), "| Control Genes:", length(control_genes)))
    
    # 저장용 리스트 초기화
    results_list <- list()
    results_list[["Raw"]] <- sub_counts
    
    # ---------------------------------------------------------
    # 3.1 Basic Normalization (TMM, TPM, FPKM)
    # ---------------------------------------------------------
    # TMM
    dge <- DGEList(counts = as.matrix(sub_counts))
    dge <- calcNormFactors(dge, method = "TMM")
    results_list[["TMM_log2"]] <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
    
    # TPM
    gene_len_kb <- sub_annot$Length / 1000
    rpk <- sub_counts / gene_len_kb
    scale_factor <- colSums(rpk)
    scale_factor[scale_factor == 0] <- 1
    tpm_val <- t(t(rpk) * 1e6 / scale_factor)
    results_list[["TPM_log2"]] <- clean_matrix(log2(tpm_val + 1))
    
    # FPKM
    results_list[["FPKM_log2"]] <- clean_matrix(rpkm(as.matrix(sub_counts), 
                                                     gene.length = sub_annot$Length, 
                                                     log = TRUE, 
                                                     prior.count = 1))
    # ---------------------------------------------------------
    # 3.2 EDASeq (GC/Length) - 코호트 내부 Bias 보정
    # ---------------------------------------------------------
    eda_set <- newSeqExpressionSet(as.matrix(sub_counts),
                                   featureData = data.frame(
                                       gc = sub_annot$GC_Percent, 
                                       length = sub_annot$Length,
                                       row.names = rownames(sub_counts)
                                   ))
    
    # GC -> Length -> Sequencing Depth 순차 보정
    set_within <- withinLaneNormalization(withinLaneNormalization(eda_set, "gc", which="full"), "length", which="full")
    norm_eda   <- betweenLaneNormalization(set_within, which="upper")
    
    results_list[["EDA_Full_All"]] <- clean_matrix(log2(normCounts(norm_eda) + 1))
    
    # ---------------------------------------------------------
    # 3.3 RUVg & Proposed Full (EDA + RUV)
    # ---------------------------------------------------------
    # RUV는 TMM Log값을 인풋으로 사용 (논문 권장)
    input_ruv <- results_list[["TMM_log2"]]
    if (length(control_genes) < 5) {
        message("   [WARNING] Too few control genes found. Skipping RUV.")
    } else {
        for (k in c(1, 2, 3)) {
            # RUV 계산
            set_ruvg <- RUVg(input_ruv, control_genes, k = k, isLog = TRUE, center = TRUE)
            # W (Noise Factor) 저장
            W_platelet <- set_ruvg$W
            write.csv(W_platelet, file.path(out_base_dir, paste0(cohort, "_W_factor_k", k, ".csv")))
            # 1) 순수 RUV 결과
            key_ruvg <- paste0("RUVg_Platelet_k", k)
            results_list[[key_ruvg]] <- clean_matrix(set_ruvg$normalizedCounts)
            # 2) Proposed: EDA로 GC 잡고 + RUV로 Platelet 잡기 (Within Cohort)
            key_full <- paste0("Proposed_Full_k", k)
            results_list[[key_full]] <- clean_matrix(
                removeBatchEffect(results_list[["EDA_Full_All"]], covariates = W_platelet)
            )
        }
    }
    # ---------------------------------------------------------
    # 4. 저장 (Cohort 접두사 붙여서 저장)
    # ---------------------------------------------------------
    for(name in names(results_list)) {
        # 파일명 예시: Norm_Seq1_TMM_log2.csv
        file_path <- paste0(out_base_dir, "Norm_", cohort, "_", name, ".csv")
        write.csv(as.data.frame(results_list[[name]]), file = file_path, row.names = TRUE)
    }
    message(paste("   Saved files for:", cohort))
}

summary_files <- list.files(path = base_path, 
                            pattern = "_featureCounts_summary.tsv$", 
                            recursive = TRUE, 
                            full.names = TRUE)

print(paste("Found summary files:", length(summary_files)))
summary_combined <- NULL

for (i in 1:length(summary_files)) {
  file_path <- summary_files[i]
    current_sum <- read.delim(file_path, header=TRUE, sep="\t", 
                            check.names=FALSE, stringsAsFactors=FALSE)
  if (is.null(summary_combined)) {
    summary_combined <- current_sum
  } else {
    summary_combined <- merge(summary_combined, current_sum, by = 1, all = TRUE)
  }
}

summary_combined[is.na(summary_combined)] <- 0
rownames(summary_combined) <- summary_combined[, 1]
summary_combined <- summary_combined[, -1]
summary_t <- as.data.frame(t(summary_combined))

meta_new <- merge(meta_filtered, summary_t, 
                  by = "row.names", all.x = TRUE)
head(meta_new)
dim(meta_new)
write.csv(meta_new, paste0(out_base_dir, "Meta_Processed.csv"))
write.csv(annot, paste0(out_base_dir, "Annot.csv"))

