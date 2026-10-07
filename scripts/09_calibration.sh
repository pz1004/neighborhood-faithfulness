#!/usr/bin/env bash
# Stage 9 -- every real result's observable status and reporting class (Gate A), then the
# study-level calibration, Layer C: 1,000 simulated studies in each of six scenarios
# (docs/revision_protocol.md, Sections 3.4 to 5). About 19 hours on 36 threads.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
python3 analyses/calibration/analyze_gate_a.py    # -> results/gate_a/classification.json, classification.txt
python3 tests/test_layer_c.py > results/revision/test_layer_c.log
bash analyses/calibration/run_layer_c.sh          # the six runs -> results/revision/layer_c/
grep -qx "layer C exit 0" results/revision/layer_c/chain_done.txt \
  || { echo "Layer C failed; see results/revision/layer_c/*.log" >&2; exit 1; }
python3 analyses/calibration/analyze_layer_c_rule.py    # -> results/revision/layer_c/rule.txt, rule.json
