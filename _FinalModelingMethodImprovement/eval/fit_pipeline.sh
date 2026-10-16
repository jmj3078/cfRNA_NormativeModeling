#!/bin/bash
# Fit one pipeline from scratch, in two steps:
#   bash eval/fit_pipeline.sh <pipeline> calib   EB calibration fits, PCIS null, rate table; put the 1e-4 cut
#                                                into config.PIPELINES[<pipeline>]["pcis_cut"]
#   bash eval/fit_pipeline.sh <pipeline> fit     full engine, CV, LOBO, then eval/build_results.sh
# Logs: runs/<pipeline>/log_fit_*.txt
set -e
export FMMI_PIPELINE=$1
PY=${PY:-~/miniconda3/envs/scRNA/bin/python}
RSCRIPT=$(dirname $PY)/Rscript
cd "$(dirname "$0")/../.."
F=$PWD/_FinalModelingMethodImprovement
RUN=$F/runs/$1
mkdir -p $RUN/pcis_calibration
if [ "$2" = calib ]; then
  $PY $F/engine/core/run_engine.py --calib-only > $RUN/log_fit_calib.txt 2>&1
  $RSCRIPT $F/engine/core/pcis_null.R $RUN/full/calib_fits.csv $RUN/tmp/glmm_train $RUN/pcis_calibration/null_pcis.csv \
    1 12 $RUN/full/disp_prior.json >> $RUN/log_fit_calib.txt 2>&1
  $PY $F/engine/core/pcis_calibration.py $RUN/pcis_calibration/null_pcis.csv | tee -a $RUN/log_fit_calib.txt
elif [ "$2" = fit ]; then
  $PY $F/engine/core/run_engine.py > $RUN/log_fit_full.txt 2>&1
  $PY $F/engine/validation/cv_engine.py > $RUN/log_fit_cv.txt 2>&1
  $PY $F/engine/validation/lobo_engine.py > $RUN/log_fit_lobo.txt 2>&1
  bash $F/eval/build_results.sh $1
else
  echo "usage: bash eval/fit_pipeline.sh <pipeline> calib|fit"; exit 1
fi
