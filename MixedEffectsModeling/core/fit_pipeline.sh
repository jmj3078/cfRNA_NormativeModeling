#!/bin/bash
# Fit one pipeline of config.PIPELINES from scratch, in two steps:
#   bash MixedEffectsModeling/core/fit_pipeline.sh <pipeline> calib   EB calibration fits, PCIS null, rate table; put the
#                                                                     printed 1e-4 cut into config.PIPELINES[<pipeline>]["pcis_cut"]
#   bash MixedEffectsModeling/core/fit_pipeline.sh <pipeline> fit     full engine, CV, LOBO
# Logs: Notebook_<pipeline>/logs/
set -e
export NM_PIPELINE=$1
PY=${PY:-~/miniconda3/envs/scRNA/bin/python}
cd "$(dirname "$0")/../.."
M=$PWD/MixedEffectsModeling
cfg() { $PY -c "import sys; sys.path.insert(0, '.'); import config; print(config.$1)"; }
LOG=$(cfg PIPELINE_DIR)/logs
mkdir -p $LOG
if [ "$2" = calib ]; then
  ENGINE=$(cfg ENGINE_MIXED_DIR); CAL=$(cfg PCIS_CAL_DIR)
  mkdir -p $CAL
  $PY $M/core/run_engine.py --calib-only > $LOG/fit_calib.txt 2>&1
  $(cfg RSCRIPT) $M/core/pcis_null.R $ENGINE/calib_fits.csv $(cfg TMP_ROOT)/glmm_train $CAL/null_pcis.csv \
    5 12 $ENGINE/disp_prior.json >> $LOG/fit_calib.txt 2>&1
  $PY $M/core/pcis_calibration.py $CAL/null_pcis.csv | tee -a $LOG/fit_calib.txt
  $PY -c "import pandas as pd; t = pd.read_csv('$CAL/pcis_rate_table.csv'); print('pcis_cut (1e-4):', t.loc[t.target_per_obs_rate == 1e-4, 'cut'].item())" | tee -a $LOG/fit_calib.txt
elif [ "$2" = fit ]; then
  $PY $M/core/run_engine.py > $LOG/fit_full.txt 2>&1
  $PY $M/validation/cv_engine.py > $LOG/fit_cv.txt 2>&1
  $PY $M/validation/lobo_engine.py > $LOG/fit_lobo.txt 2>&1
else
  echo "usage: bash MixedEffectsModeling/core/fit_pipeline.sh <pipeline> calib|fit"; exit 1
fi
