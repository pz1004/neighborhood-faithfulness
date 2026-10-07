#!/usr/bin/env bash
# Stage 5 -- mouse V1 on the Allen two-photon recordings (docs/note_protocol_allen.md,
# plan Amendment 4). About 3 minutes.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_allen_confirmatory.py > results/allen_confirmatory/test_allen_confirmatory.log
python3 analyses/mouse_v1/run_allen_confirmatory.py > results/allen_confirmatory/run.log
python3 analyses/mouse_v1/analyze_allen_confirmatory.py    # -> allen_analysis.txt, allen_summary.json
python3 analyses/mouse_v1/explore_allen_snr.py > results/allen_confirmatory/explore_snr.txt
python3 analyses/mouse_v1/describe_allen_subjects.py       # -> subjects_described.json
