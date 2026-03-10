#!/bin/bash
TARGET_DIR="/project/cfRNA_Disentaglement/Data/OpenAccess/Processed"

# 1. 삭제 대상 파일 목록 사전 확인 (Dry-run)
echo "=== 삭제 대상 파일 목록 ==="
find "$TARGET_DIR" -type f -name "*.csv" ! -name "*Plasma*" ! -name "*Serum*"

# 2. 사용자 확인
find "$TARGET_DIR" -type f -name "*.csv" ! -name "*Plasma*" ! -name "*Serum*" -exec rm -v {} +