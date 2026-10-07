#!/usr/bin/env python3
"""Checks of the revision's real-data code (docs/revision_protocol.md, Section 6) on synthetic data.

No data file is opened. It checks the following:
- run_points.scores_for_budget gives exactly the labels of reach/checks.labels_for_budget, with
  the same generators, and its scores reproduce them by the published rules;
- the continuous overlap reaches m exactly where the reach is at most k';
- Kruskal stress-1 is 0 for an exact fit;
- the session-disjoint blocks and their control have the stated layouts;
- the EEG manifest's chronological ranks are a permutation within every session.
The bootstrap uses 5 maps instead of 50 here, to keep the run short.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

from reach import checks, things, target as tg
from analyses.benchmark import run_points as RP
from analyses.benchmark import run_session_disjoint as SD
from analyses.benchmark import make_trial_manifest as TM

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_090_005)
checks.B_BOOT = 5

# scores_for_budget against labels_for_budget: a 4-dimensional signal in 30 raw features, 8 trials
n, R, p = 80, 8, 30
Z = rng.normal(size=(n, 4))
mu = 3.0 * Z @ rng.normal(size=(4, p)) / 2.0
raw = mu[:, None, :] + rng.normal(size=(n, R, p))
W, m0 = rng.normal(size=(p, 20)) / np.sqrt(p), raw.mean(axis=(0, 1))
project = lambda A: (A - m0) @ W
Xu, c_raw = project(raw), raw.mean(axis=1)
seeds = (11, 12, 13)
gens = lambda base: [np.random.default_rng((base, j)) for j in range(3)]
ref = checks.labels_for_budget(Xu, c_raw, project, things.halves(R), seeds, gens(1), gens(2), things.cv_sq_dist)
new = RP.scores_for_budget(Xu, c_raw, project, things.halves(R), seeds, gens(1), gens(2), gens(3))
same = all((E1 == E2).all() and all((l1[X] == l2[X]).all() for X in checks.NAMES)
           for (E1, l1), (E2, l2, _) in zip(ref, new))
check("per-point rerun: every check's labels equal labels_for_budget's, with the same generators", same,
      "passes per check: " + ", ".join(f"{X} {int(new[0][1][X].sum())}" for X in checks.NAMES))
sc = new[0][2]
rules = ((sc["naive"] <= RP.KP) == new[0][1]["N"]).all() and ((sc["scd"] >= sc["scd_cut"]) == new[0][1]["SCD"]).all() \
    and ((sc["cv"] <= RP.KP) == new[0][1]["CV"]).all() and ((sc["seed_count"] == 2) == new[0][1]["SEED"]).all() \
    and ((sc["shA"] & sc["shB"]) == new[0][1]["SH"]).all()
check("per-point rerun: the stored scores give the labels by the published rules", rules)

# continuous overlap against the reach
P5 = rng.normal(size=(n, 5))
E_idx = tg.displayed(P5 + rng.normal(scale=0.5, size=P5.shape), 10)[0]        # some faithful, some not
rk = tg.rank_matrix(tg.pairwise(P5))
reach = tg.reach_from_ranks(rk, E_idx, RP.M)
overlap = (np.take_along_axis(rk, E_idx, axis=1) <= RP.KP).sum(axis=1)
check("the continuous overlap reaches m exactly where the reach is at most k'",
      ((overlap >= RP.M) == (reach <= RP.KP)).all(), f"{int((reach <= RP.KP).sum())} faithful of {n}")

# Kruskal stress-1
pts = rng.normal(size=(30, 2))
D = tg.pairwise(pts)
np.fill_diagonal(D, 0.0)
check("Kruskal stress-1 is 0 for an exact fit and positive otherwise",
      RP.stress1(D, pts) < 1e-12 and RP.stress1(D, pts + rng.normal(scale=0.1, size=pts.shape)) > 0)

# session-disjoint blocks and their control
ses = things.sessions()
b = SD.disjoint_blocks(ses, 0, 2, np.random.default_rng(0))
ok = (ses[b["F"]] == 0).all() and len(b["F"]) == 8 and (ses[b["U"]] == 2).all() and len(b["U"]) == 16 \
    and len(b["H"]) == 40 and set(ses[b["H"]]) == {1, 3} \
    and all((np.bincount(ses[b["H"][j * Rh:(j + 1) * Rh]], minlength=4)[[1, 3]] == Rh // 2).all()
            for Rh in (8, 20, 40) for j in range(40 // Rh)) \
    and len(set(b["F"]) | set(b["U"]) | set(b["H"])) == 64
check("session-disjoint blocks: F from one session, U from another, H balanced over the other two", ok)
c = SD.control_blocks(ses, np.random.default_rng(0))
ok = all(len(c[k]) == m and (np.bincount(ses[c[k]], minlength=4) == m // 4).all() for k, m in (("F", 8), ("U", 16), ("H", 40))) \
    and len(set(c["F"]) | set(c["U"]) | set(c["H"])) == 64 \
    and all((np.bincount(ses[c["U"][:Rr]], minlength=4) == Rr // 4).all() for Rr in (8, 16))
check("control blocks: 2, 4 and 10 repetitions per session, every prefix of U balanced", ok)

# the EEG manifest
rows, perm20 = TM.eeg_manifest()
ok = sorted(perm20.tolist()) == list(range(20)) and all(
    sorted(r["chronological_rank_in_session"] for r in rows if r["subject"] == 1 and r["session"] == s_) == list(range(20))
    for s_ in range(1, 5))
check("EEG manifest: chronological ranks are a permutation of 0-19 within every session", ok, f"{len(rows)} rows")
check("EEG manifest: every repetition has a block in each confirmatory split",
      all(r[f"split{j}"] in ("F", "U", "H") for r in rows for j in range(10)))

print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
