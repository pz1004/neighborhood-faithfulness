#!/usr/bin/env bash
# Stage 6 -- six common checks against the truth, in both datasets (docs/checks_protocol.md,
# plan Amendment 5). About 75 minutes.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_checks.py > results/checks/test_checks.log
python3 analyses/checks/run_checks.py --dataset things > results/checks/things/run.log
python3 analyses/checks/run_checks.py --dataset allen > results/checks/allen/run.log
python3 analyses/checks/analyze_checks.py    # -> results/checks/checks_analysis.txt, checks_summary.json
