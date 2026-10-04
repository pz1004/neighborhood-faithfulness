#!/usr/bin/env python3
"""Checks of the Layer C code (run_study_calibration.py, analyze_study_calibration.py) on a synthetic base.

No data file is opened: the base subject is replaced by a synthetic one in the THINGS layout. It
checks the following:
- a study unit returns the published procedure's pooled counts (10 splits x 3 seeds), deterministically;
- its plasmode carries the guard's counts;
- an oracle unit returns exact counts from three maps;
- the analysis runs end to end on a few studies, with the bounds ordered and the truth recovered.
"""
import os, sys, gzip, json, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

from reach import sim
import run_study_calibration as SC
import analyze_study_calibration as AS

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_090_006)
Z = rng.normal(size=(200, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
mu = 1.2 * Z @ rng.normal(size=(4, 60)) / 2.0
feat = mu[:, None, :] + rng.normal(size=(200, 80, 60))
synthetic = sim.Base("EEG", 1, feat=feat)
SC.base_of = lambda ds, s: synthetic                      # every virtual subject uses the synthetic base

u1 = SC.study_unit("EEG", "T1", "tsne", 0, 0, 0, guard=True)
u2 = SC.study_unit("EEG", "T1", "tsne", 0, 0, 0, guard=False)
check("a study unit pools 10 splits x 3 seeds against every published reference size and the truth",
      sorted(u1["counts"]) == sorted(["plug|8", "plug|20", "plug|40", "exact|0"])
      and u1["counts"]["exact|0"][4] == 30 * 200, f"{u1['counts']['exact|0'][5]:.0f} faithful calls pooled")
check("a study unit is deterministic in its seed", u1["counts"] == u2["counts"])
check("the plasmode carries the guard's pooled counts over 20 draws",
      u1["plasmode"] is not None and u1["plasmode"]["exact|0"][4] == 20 * 200 and u2["plasmode"] is None)
o = SC.oracle_unit("EEG", "T1", "tsne", 0, 0)
check("an oracle unit counts three maps against the truth", o["n1"] > 0 and 0 <= o["n10"] <= o["n1"], str(o))

# the analysis on three studies of two subjects each, and a five-draw oracle
tmp = tempfile.mkdtemp()
SC.OUT = tmp
SC.SUBJECTS = 2
with gzip.open(f"{tmp}/oracle_EEG_T1_tsne.jsonl.gz", "wt") as f:
    for v in range(3):
        for d in range(5):
            f.write(json.dumps(SC.oracle_unit("EEG", "T1", "tsne", v, d)) + "\n")
with gzip.open(f"{tmp}/studies_EEG_T1_tsne.jsonl.gz", "wt") as f:
    for st in range(3):
        for j in range(2):
            f.write(json.dumps(SC.study_unit("EEG", "T1", "tsne", st, j, (st + j) % 3, guard=st < 2)) + "\n")
sys.argv = ["analyze_study_calibration.py", "--scenario", "T1"]
AS.main()
res = json.load(open(f"{tmp}/calibration_EEG_T1_tsne.json"))
check("the analysis runs end to end, with the guard computed where the plasmode ran",
      res["studies"] == 3 and res["guard"]["studies"] == 2 and 0 <= res["Theta"] <= 1,
      f"Theta {res['Theta']:.3f}, study medians sqrt {res['medians']['sqrt']:.3f}")
check("the estimates land near the truth in this easy regime", abs(res["medians"]["sqrt"] - res["Theta"]) < 0.1)

print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
