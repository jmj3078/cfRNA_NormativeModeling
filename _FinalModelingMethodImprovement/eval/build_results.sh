#!/bin/bash
# Rebuild every result row for one pipeline: bash eval/build_results.sh <pipeline>
# Needs runs/<pipeline>/{full,cv,lobo} and the frozen splits / in-silico cells.
set -e
export FMMI_PIPELINE=$1
PY=${PY:-~/miniconda3/envs/scRNA/bin/python}
cd "$(dirname "$0")/../.."
E=_FinalModelingMethodImprovement/eval
LOG=_FinalModelingMethodImprovement/runs/$1
[ -f _FinalModelingMethodImprovement/results/insilico_sweep.csv ] || $PY $E/insilico_sweep.py import
[ -f $LOG/full/scores/Z_hc_shash.npy ] || $PY $E/score_full.py > $LOG/log_score_full.txt 2>&1
$PY $E/insilico_sweep.py score > $LOG/log_insilico.txt 2>&1
$PY $E/adaptation.py insilico > $LOG/log_adapt.txt 2>&1
$PY $E/detection_limit.py > $LOG/log_dl.txt 2>&1
$PY $E/detection_limit.py --adapt > $LOG/log_dl_adapt.txt 2>&1
echo "results rebuilt for $1"
