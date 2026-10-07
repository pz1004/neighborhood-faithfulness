#!/usr/bin/env python3
"""Checks of the second-dataset code (docs/note_protocol_allen.md) on synthetic files.

No Allen recording is opened. A synthetic HDF5 file with the recordings' NWB 1.0 layout is
written with known per-presentation responses. It checks the following:
- the loader recovers those responses exactly;
- the blocks and halves are balanced over the strata;
- every registered cell runs;
- the extrapolated false comfort lands near the exact value computed from the known truth.
"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import h5py
from reach import allen
from analyses.mouse_v1 import run_allen_confirmatory as run
from analyses.mouse_v1 import analyze_allen_confirmatory as ana

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_090_003)
n_cells, fps, gap = 120, 30.0, 20              # presentations 20 frames apart: windows never overlap
Z = rng.normal(size=(118, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
mu = 0.6 * Z @ rng.normal(size=(4, n_cells)) / 2.0          # true mean response per image
order = rng.permutation(np.repeat(np.arange(-1, 118), 50))  # images and blanks, random order
starts = 30 + gap * np.arange(len(order))
n_frames = int(starts[-1] + 60)
trial = rng.normal(size=(len(order), n_cells))              # trial noise
dff = np.zeros((n_cells, n_frames), np.float32)
for p, (i, s) in enumerate(zip(order, starts)):
    if i >= 0:
        dff[:, s:s + 15] = (mu[i] + trial[p])[:, None]
with tempfile.TemporaryDirectory() as d:
    path = os.path.join(d, "synthetic.nwb")
    with h5py.File(path, "w") as f:
        f[allen.DFF + "/data"] = dff
        f[allen.DFF + "/timestamps"] = np.arange(n_frames) / fps
        f[allen.NS + "/data"] = order.astype(np.int64)
        f[allen.NS + "/frame_duration"] = np.column_stack([starts, starts + 7]).astype(np.int64)
    feats, info = allen.load_responses(path)
expect = np.empty((118, 50, n_cells))
for i in range(118):
    idx = np.flatnonzero(order == i)                         # already in time order
    expect[i] = mu[i] + trial[idx]
check("loader recovers per-presentation responses (0.5 s window)", np.allclose(feats["full"], expect, atol=1e-5),
      f"fps {info['fps']:.1f}, {info['n_cells']} cells")
check("0.25 s window uses the first 8 frames",
      np.allclose(feats["early"], expect, atol=1e-5))   # the planted response is constant over 15 frames

st = allen.strata()
b = allen.split_blocks(st, np.random.default_rng(0))
check("blocks partition the 50 repetitions",
      sorted(np.concatenate([b["F"], b["U"], b["H"]]).tolist()) == list(range(50)))
bal = all((np.bincount(st[b["U"][:R]], minlength=2) == R // 2).all() for R in (8, 16))
bal &= all((np.bincount(st[b["H"][j:j + Rh]], minlength=2) == Rh // 2).all()
           for Rh in (8, 12, 24) for j in range(0, 24, Rh))
bal &= all((np.bincount(st[b["U"][:R][h]], minlength=2) == R // 4).all() for R in (8, 16) for h in allen.halves(R))
check("user prefixes, reference sub-blocks and halves are stratum-balanced", bal)

truth = {"full": mu, "early": mu}
rows = []
for split in range(4):
    r, _ = run.split_rows(1, split, feats, st, truth=truth,
                          cells=[(nm, c) for nm, c in run.CELLS if nm == "primary"], n_seeds=1)
    rows += r
est = ana.estimates(ana.disagreement([r for r in rows if r["Rh"] != "exact"]))
for R in (8, 16):
    ex = [r for r in rows if r["Rh"] == "exact" and r["R"] == R]
    fc = sum(r["n10"] for r in ex) / max(1, sum(r["n1"] for r in ex))
    e1, e2 = est[(1, "primary", R)]["FC"]
    print(f"   primary R={R}: NF {est[(1, 'primary', R)]['NF']:.3f}  exact FC {fc:.3f}  estimates {e1:.3f} / {e2:.3f}")
    check(f"primary R={R}: an extrapolant lands within 0.10 of the exact false comfort",
          min(abs(e1 - fc), abs(e2 - fc)) < 0.10)

rows_all, b0 = run.split_rows(1, 0, feats, st, truth=truth, R_list=(8,), n_seeds=1)
check("all 12 registered cells run", {r["cell"] for r in rows_all} == {nm for nm, _ in run.CELLS})
est_all = ana.estimates(ana.disagreement([r for r in rows_all if r["Rh"] != "exact"]))
check("the analysis summary runs on every cell's output",
      np.isfinite(ana.summarize(est_all, "primary", 8, "FC", np.random.default_rng(0))["L"]))
run.N_SIM = 2
bias = ana.semi_synthetic_bias(run.semi_synthetic(1, feats, b0))
check("semi-synthetic check and bias summary run", np.isfinite(bias["q"]["median_bias"]),
      f"bias FC {bias['FC']['median_bias']:+.3f}, q {bias['q']['median_bias']:+.3f} (2 draws)")
print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
