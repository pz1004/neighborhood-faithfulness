#!/usr/bin/env bash
# Stage 3 -- the exploratory THINGS-EEG2 pilot (plan Amendment 2; docs/gate3_pilot_things.md).
# About 3 minutes.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_things.py > results/things_pilot/test_things.log
python3 analyses/eeg_pilot/run_things_need_pilot.py > results/things_pilot/run.log
python3 analyses/eeg_pilot/analyze_things_need_pilot.py    # -> results/things_pilot/pilot_analysis.txt, pilot_summary.json
python3 analyses/eeg_pilot/explore_things_snr.py > results/things_pilot/explore_snr.txt
