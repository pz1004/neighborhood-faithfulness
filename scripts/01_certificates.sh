#!/usr/bin/env bash
# Stage 1 -- per-point certificates on synthetic data (docs/research_plan.md, Gates 1 and 2).
# Only scikit-learn's bundled digits; about 15 minutes on 36 threads.
# run_gate2_validity.py resumes from results/gate2/validity.jsonl when that file exists;
# delete it to recompute every draw.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_reach.py > results/gate2/test_reach.log
python3 analyses/certificates/run_gate2_validity.py --workers 30 > results/gate2/validity.log
python3 analyses/certificates/run_gate2_power.py --workers 32 > results/gate2/power.log
python3 analyses/certificates/analyze_gate2.py     # -> results/gate2/{validity,gate2}_analysis.txt, gate2_summary.json
python3 analyses/certificates/run_gate2_synthetic_ceiling.py > results/gate2/synthetic_ceiling.txt
python3 analyses/certificates/run_gate2_need_crosstab.py > results/gate2/need_crosstab.txt
