#!/bin/bash
# Layer C under the local reading (docs/revision_addendum.md, Section 2): EEG C1 (false comfort, R = 8, with
# the guard) under T1 and T3; mouse V1 C3 (mislabel rate, R = 8 and 16, no guard) under T1 and T3.
# This is the chain that produced results/revision/layer_c/ (run 2026-10-04 to 05; about 19 h on 36 cores).
cd "$(dirname "$0")/../.." || exit 1
L=results/revision/layer_c
mkdir -p $L
run() {  # dataset scenario R guard-flag
  local tag="$1_$2_tsne_R$3"
  python3 analyses/calibration/run_study_calibration.py --dataset $1 --scenario $2 --R $3 --oracle > $L/oracle_$tag.log 2>&1 || return 1
  python3 analyses/calibration/run_study_calibration.py --dataset $1 --scenario $2 --R $3 $4 --studies 1000 > $L/studies_$tag.log 2>&1 || return 1
  python3 analyses/calibration/analyze_study_calibration.py --dataset $1 --scenario $2 --R $3 > /dev/null 2>&1 || return 1
}
run EEG T1 8 "" && run EEG T3 8 "" && run mouse T1 8 --no-guard && run mouse T1 16 --no-guard \
  && run mouse T3 8 --no-guard && run mouse T3 16 --no-guard
echo "layer C exit $?" > $L/chain_done.txt
