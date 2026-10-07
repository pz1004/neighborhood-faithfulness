#!/usr/bin/env python3
"""Checks of the confirmatory code (docs/note_protocol.md) on synthetic data with a known truth.

No THINGS-EEG2 file is opened. The synthetic set has the THINGS layout (200 images, 80
repetitions in 4 sessions), Gaussian noise, and a low-dimensional signal. Because the truth
is known, the extrapolated estimates can be compared with the exact values that the same
code computes.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from reach import things
from analyses.eeg_confirmatory import run_things_confirmatory as run
from analyses.eeg_confirmatory import analyze_things_confirmatory as ana

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_090_002)

# cross-validated distances are unbiased; halves are session-stratified
mu = rng.normal(size=(30, 6))
d_true = ((mu[:, None] - mu[None]) ** 2).sum(-1)
acc = np.zeros_like(d_true)
for _ in range(400):
    Xb = mu[:, None, :] + 2.0 * rng.normal(size=(30, 16, 6))
    d2 = things.cv_sq_dist(Xb, *things.halves(16))
    np.fill_diagonal(d2, 0.0)
    acc += d2
iu = np.triu_indices(30, 1)
err = np.abs(acc / 400 - d_true)[iu].max()
check("cross-validated squared distances are unbiased", err < 1.5, f"max |mean - truth| = {err:.2f} for noise sd 2")
ses = things.sessions()
b = things.split_blocks(ses, np.random.default_rng(0))
ok = all((np.bincount(ses[b["U"][:R][h]], minlength=4) == R // 8).all()
         for R in (8, 16, 32) for h in things.halves(R))
check("halves of every user prefix hold R/8 repetitions per session", ok)
t = np.arange(-20, 80) / 100
check("time windows select 80 and 46 samples",
      things.features_window(np.zeros((1, 1, 17, 100)), t, 0.0, 0.79).shape[-1] == 80 * 17
      and things.features_window(np.zeros((1, 1, 17, 100)), t, 0.05, 0.50).shape[-1] == 46 * 17)

# end to end with a known truth: a 4-dimensional signal in 60 features, unit noise per trial
n, p = 200, 60
Z = rng.normal(size=(n, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
mu = 2.5 * Z @ rng.normal(size=(4, p)) / np.sqrt(4)
X = mu[:, None, :] + rng.normal(size=(n, 80, p))
feats, truth = {"full": X, "early": X[:, :, :40]}, {"full": mu, "early": mu[:, :40]}
cells = [(nm, c) for nm, c in run.CELLS if nm in ("primary", "rsa_cv_mds", "pca2", "unwhitened", "early_window")]
rows = []
for split in range(4):
    r, _ = run.split_rows(1, split, feats, ses, truth=truth, cells=cells, R_list=(8, 16), n_seeds=1)
    rows += r
ref = [r for r in rows if r["Rh"] != "exact"]
check("counts are internally consistent",
      all(r["n10"] <= r["n1"] * r["nb"] and r["dis"] <= r["n"] * r["nb"] and r["t1"] <= r["n"] * r["nb"] for r in rows))
dis = ana.disagreement(ref)
est = ana.estimates(dis)
for nm, _ in cells:
    for R in (8, 16):
        ex = [r for r in rows if r["Rh"] == "exact" and r["cell"] == nm and r["R"] == R]
        n1, n10 = sum(r["n1"] for r in ex), sum(r["n10"] for r in ex)
        fc_exact = n10 / n1 if n1 else float("nan")
        e1, e2 = est[(1, nm, R)]["FC"]
        d8, d40 = dis[(1, nm, R)][8][1], dis[(1, nm, R)][40][1]
        print(f"   {nm:13s} R={R:2d}: NF {est[(1, nm, R)]['NF']:.3f}  exact FC {fc_exact:.3f}  "
              f"d_FC(8) {d8:.3f} d_FC(40) {d40:.3f}  estimates {e1:.3f} / {e2:.3f}")
        if nm == "primary":
            check(f"primary R={R}: the reference moves toward the truth as R_h grows",
                  abs(d40 - fc_exact) <= abs(d8 - fc_exact) + 0.02)
            check(f"primary R={R}: an extrapolant lands within 0.10 of the exact false comfort",
                  min(abs(e1 - fc_exact), abs(e2 - fc_exact)) < 0.10)

# the semi-synthetic check and its bias summary run on the same layout
run.N_SIM = 3
_, b0 = run.split_rows(1, 0, feats, ses, cells=cells[:1], R_list=(8,), n_seeds=1)
sim = run.semi_synthetic(1, feats, b0)
bias = ana.semi_synthetic_bias(sim)
check("semi-synthetic check and bias summary run", np.isfinite(bias["q"]["median_bias"]),
      f"bias FC {bias['FC']['median_bias']:+.3f}, q {bias['q']['median_bias']:+.3f} (3 draws)")
# every registered cell runs, and the analysis summary runs on its output
rows_all, _ = run.split_rows(1, 0, feats, ses, truth=truth, R_list=(8,), n_seeds=1)
names = {r["cell"] for r in rows_all}
check("all 12 registered cells run", names == {nm for nm, _ in run.CELLS}, f"{len(names)} cells")
est_all = ana.estimates(ana.disagreement([r for r in rows_all if r["Rh"] != "exact"]))
summ = ana.summarize(est_all, "primary", 8, "FC", np.random.default_rng(0))
check("the analysis summary runs on every cell's output",
      all(np.isfinite(est_all[(1, nm, 8)]["NF"]) for nm, _ in run.CELLS) and np.isfinite(summ["L"]))
print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
