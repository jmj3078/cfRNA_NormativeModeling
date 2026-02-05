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
counts_a <- read.table("/project/cfRNA_Disentaglement/Data/RPM/RPM_Lab/RPM10_12_13_14_16_raw.csv", header=TRUE, row.names=1, sep=',')
counts_b <- read.table("/project/cfRNA_Disentaglement/Data/RPM/Compgen_Lab/RPM15_seq.csv", header=TRUE, row.names=1, sep=',')
meta     <- read.table("/project/cfRNA_Disentaglement/Data/RPM/RPM_Lab/Meta.csv", header=TRUE, row.names=1, sep=',')
annot    <- read.table("/project/cfRNA_Disentaglement/Data/GECODEv49_Annot.tsv", header=TRUE, row.names=1, sep='\t')
palangodb <- read.table("/project/cfRNA_Disentaglement/Data/PalangoDB_CellTypeMarkers.tsv", header=TRUE, sep='\t')

out_base_dir <- "/project/cfRNA_Disentaglement/Data/RPM/Processed/"
if(!dir.exists(out_base_dir)) dir.create(out_base_dir, recursive = TRUE)

# Count Matrix 병합 및 정렬
counts <- merge(counts_a, counts_b, by = 0)
rownames(counts) <- counts$Row.names
counts$Row.names <- NULL
counts <- counts[, order(colnames(counts))]
sample_sums <- colSums(counts)

# 총합이 0인 샘플 찾기
zero_samples <- names(sample_sums)[sample_sums == 0]

if (length(zero_samples) > 0) {
    message(paste("!!! 경고: 총 Read Count가 0인 샘플이 발견되어 제거합니다:", length(zero_samples), "개"))
    print(zero_samples) # 어떤 샘플이 지워지는지 확인
    
    # 0보다 큰 샘플만 남김
    counts <- counts[, sample_sums > 0]
} else {
    message(">> 모든 샘플이 정상입니다 (Read Count > 0).")
}

# -----------------------------------------------------------------------------
message("Data Loaded:")
message(paste("Total Samples:", ncol(counts)))
message(paste("Total Genes:", nrow(counts)))
write.csv(counts, paste0(out_base_dir, "featureCounts_merged.csv"))
# -----------------------------------------------------------------------------
# STEP 2: 전처리 (Metadata & Gene Filtering)
# -----------------------------------------------------------------------------
# 2.1 Metadata 정제 (Python 로직 반영)
meta <- meta %>%
  mutate(
    Batch_Granular = paste(Sample.source, Seq_ID, Batch_tube, Batch_centrifuge, Batch_RNAext, sep = "-"),
    Subtype = ifelse(Subtype == Type | Subtype == "(NA)", NA, Subtype),
    Responder = ifelse(Responder == 1, "ICI-Responder", "ICI-Nonresponder")
  ) %>%
  unite("Type_Granular", Type, Subtype, Responder, sep = "_", remove = FALSE, na.rm = TRUE)

# 2.2 샘플 동기화 (Intersection)
common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta   <- meta[common_samples, ]

# 데이터 정합성 체크
stopifnot(all(colnames(counts) == rownames(meta)))

# 2.3 Protein Coding & Valid Annotation 필터링 (Reference와 동일 로직)
pc_genes  <- rownames(annot)[annot$GeneType == "protein_coding"]
counts_pc <- counts[rownames(counts) %in% pc_genes, ]
annot_pc  <- annot[rownames(counts_pc), ] # 순서 동기화

# 길이(Length)와 GC 함량이 유효한 유전자만 남김 (EDASeq 필수 조건)
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
cohort_list <- unique(meta$Seq_ID)
cohort_list

clean_matrix <- function(mat) {
    mat[is.na(mat)] <- 0
    mat[is.infinite(mat)] <- 0
    return(mat)
}

message(paste("Start processing for", length(cohort_list), "cohorts:", paste(cohort_list, collapse=", ")))

for (cohort in cohort_list) {
    
    message(paste0("\n>>> Processing Cohort: ", cohort, " <<<"))
    
    # 3.1 Subset Data (해당 코호트 샘플만 추출)
    cohort_samples <- rownames(meta)[meta$Seq_ID == cohort]
    sub_counts <- counts_pc[, cohort_samples, drop=FALSE]
    
    # 3.2 Low Expression Gene Filtering (이 코호트 내에서 발현 안되는 유전자 제거)
    # Reference 로직: 해당 배치 내에서 RowSum > 0 인 유전자만 유지
    keep_genes <- rowSums(sub_counts) > 0
    sub_counts <- sub_counts[keep_genes, ]
    sub_annot  <- annot_pc[keep_genes, ]
    
    # Platelet Gene도 현재 살아남은 유전자와 교집합 업데이트
    control_genes <- intersect(platelet_ids, rownames(sub_counts))
    
    message(paste("   Samples:", ncol(sub_counts), "| Genes:", nrow(sub_counts), "| Control Genes:", length(control_genes)))
    
    # 결과 저장용 리스트
    results_list <- list()
    results_list[["Raw"]] <- sub_counts
    
    # ---------------------------------------------------------
    # 3.3 Basic Normalization (TMM, TPM, FPKM)
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
    # 3.4 EDASeq (GC/Length) - 코호트 내부 Bias 보정
    # ---------------------------------------------------------
    eda_set <- newSeqExpressionSet(as.matrix(sub_counts),
                                   featureData = data.frame(
                                       gc = sub_annot$GC_Percent, 
                                       length = sub_annot$Length,
                                       row.names = rownames(sub_counts)
                                   ))
    
    # GC -> Length -> Sequencing Depth(Upper Quartile) 순차 보정
    set_within <- withinLaneNormalization(withinLaneNormalization(eda_set, "gc", which="full"), "length", which="full")
    norm_eda   <- betweenLaneNormalization(set_within, which="upper")
    
    results_list[["EDA_Full_All"]] <- clean_matrix(log2(normCounts(norm_eda) + 1))
    
    # ---------------------------------------------------------
    # 3.5 RUVg & Proposed Full (EDA + RUV)
    # ---------------------------------------------------------
    # RUV 입력: TMM log값 사용
    input_ruv <- results_list[["TMM_log2"]]
    
    if (length(control_genes) < 5) {
        message("   [WARNING] Too few control genes found. Skipping RUV.")
    } else {
        for (k in c(1, 2, 3)) {
            # 1) RUVg 계산
            set_ruvg <- RUVg(input_ruv, control_genes, k = k, isLog = TRUE, center = TRUE)
            
            # W (Noise Factor) 저장 -> 추후 공변량 분석용
            W_platelet <- set_ruvg$W
            write.csv(W_platelet, file.path(out_base_dir, paste0(cohort, "_W_factor_k", k, ".csv")))
            
            # RUVg 결과 저장
            key_ruvg <- paste0("RUVg_Platelet_k", k)
            results_list[[key_ruvg]] <- clean_matrix(set_ruvg$normalizedCounts)
            
            # 2) Proposed: EDA(GC/Len)로 1차 보정 + RUV(Platelet)로 2차 보정
            # limma::removeBatchEffect를 사용하여 EDA 결과에서 W 요인 제거
            key_full <- paste0("Proposed_Full_k", k)
            results_list[[key_full]] <- clean_matrix(
                removeBatchEffect(results_list[["EDA_Full_All"]], covariates = W_platelet)
            )
        }
    }
    
    # ---------------------------------------------------------
    # 4. 개별 코호트 결과 저장
    # ---------------------------------------------------------
    for(name in names(results_list)) {
        # 파일명 형식: Norm_{SeqID}_{Method}.csv
        file_path <- paste0(out_base_dir, "Norm_", cohort, "_", name, ".csv")
        write.csv(as.data.frame(results_list[[name]]), file = file_path, row.names = TRUE)
    }
    message(paste("   Saved files for:", cohort))
}

write.csv(annot, paste0(out_base_dir, "Annot_Processed.csv"))
write.csv(meta, paste0(out_base_dir, "Meta_Processed.csv"))

