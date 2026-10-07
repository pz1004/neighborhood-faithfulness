#!/usr/bin/env bash
# Stage 4 -- the confirmatory EEG analysis (docs/note_protocol.md, plan Amendment 3).
# About 8 minutes.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_things_confirmatory.py > results/things_confirmatory/test_things_confirmatory.log
python3 analyses/eeg_confirmatory/run_things_confirmatory.py > results/things_confirmatory/run.log
python3 analyses/eeg_confirmatory/analyze_things_confirmatory.py    # -> confirmatory_analysis.txt, confirmatory_summary.json
