#!/usr/bin/env bash
# Stage 7 -- the review audit's evidence, then the simulation study that validates the
# estimator (docs/revision_protocol.md, plan Amendment 6). About 2.5 hours.
# The development runs fixed the noise grid, the estimator and the thresholds in
# docs/revision_addendum.md before any validation seed was used; run_validation.py refuses
# the validation phase until the addendum and the development decisions are committed.
# Re-run, analyze_validation.py --phase dev writes dev/decisions_rerun.json and leaves the
# committed decisions alone.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 analyses/validation/analyze_review_audit.py    # -> results/review_audit/audit.txt, audit_summary.json
python3 tests/test_estimand.py > results/revision/test_estimand.log
python3 tests/test_validation.py > results/revision/test_validation.log
python3 analyses/validation/run_validation.py --phase dev --alpha    # T3's alpha from the bases; no simulation
python3 analyses/validation/run_validation.py --phase dev > results/revision/dev/run.log
python3 analyses/validation/run_validation.py --phase dev --lams 0.25,0.35 --append \
  --cells T1_l0.25,T2_l0.25,T3_l0.25,T1_l0.35,T2_l0.35,T3_l0.35 > results/revision/dev/run_extension.log
python3 analyses/validation/analyze_validation.py --phase dev --lams 0.1,0.25,0.35,0.5,1,2 > results/revision/dev/analysis.txt
python3 analyses/validation/run_validation.py --phase validation > results/revision/validation/run.log
# The committed analysis.txt also holds the script's warnings, hence 2>&1.
python3 analyses/validation/analyze_validation.py --phase validation > results/revision/validation/analysis.txt 2>&1
