#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")"

LOG_ENGINE=run_model_engine_v2.log
LOG_CV=cv_model_engine_v2.log

echo "[$(date)] START run_model_engine_v2.py" | tee "$LOG_ENGINE"
python run_model_engine_v2.py >> "$LOG_ENGINE" 2>&1
ENGINE_RC=$?
echo "[$(date)] END run_model_engine_v2.py rc=$ENGINE_RC" | tee -a "$LOG_ENGINE"

if [ $ENGINE_RC -ne 0 ]; then
    echo "[$(date)] engine training failed, skipping CV" | tee -a "$LOG_ENGINE"
    exit $ENGINE_RC
fi

echo "[$(date)] START cv_model_engine_v2.py" | tee "$LOG_CV"
python cv_model_engine_v2.py >> "$LOG_CV" 2>&1
CV_RC=$?
echo "[$(date)] END cv_model_engine_v2.py rc=$CV_RC" | tee -a "$LOG_CV"

exit $CV_RC
