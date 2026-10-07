#!/usr/bin/env bash
# Stage 8 -- the trial manifest, the session-disjoint sensitivity analysis, the per-point
# rerun and the benchmark of the checks (docs/revision_protocol.md, Section 6). About 1 hour.
# run_points.py re-creates the primary and cv-MDS maps and every check, and reports how many
# committed counts it did not reproduce (0 in the committed logs).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 tests/test_stage2.py > results/revision/test_stage2.log
python3 analyses/benchmark/make_trial_manifest.py    # -> results/manifest/
python3 analyses/benchmark/run_session_disjoint.py > results/session_disjoint/run.log
python3 analyses/benchmark/analyze_session_disjoint.py    # -> results/session_disjoint/analysis.txt, summary.json
python3 analyses/benchmark/run_points.py --dataset things > results/points/things.log
python3 analyses/benchmark/run_points.py --dataset allen > results/points/allen.log
python3 analyses/benchmark/analyze_benchmark.py    # -> results/benchmark/benchmark.txt, summary.json
