#!/usr/bin/env python3
"""Checks of the common-checks code (docs/checks_protocol.md) on synthetic data.

No THINGS-EEG2 or Allen file is opened. It checks the following:
- per-point trustworthiness averages to scikit-learn's trustworthiness;
- scDEED's reliability score matches a direct transcription of its definition;
- the recurrence and bootstrap checks behave as defined in limiting cases;
- the re-created maps reproduce the confirmatory code's plug-in counts exactly;
- every check runs end to end on both layouts, and the estimator lands near the exact
  false comfort computed from a known truth.
The bootstrap uses 5 maps instead of 50 here, to keep the run short.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from sklearn.manifold import trustworthiness

from reach import allen, things, checks, target as tg
from analyses.eeg_confirmatory import run_things_confirmatory as RT
from analyses.checks import run_checks as run
from analyses.checks import analyze_checks as ana

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_100_001)

# per-point trustworthiness
X = rng.normal(size=(120, 8))
emb = checks.tsne(X, 1)
E, _ = tg.displayed(emb, 10)
tw = checks.trustworthiness_points(tg.rank_matrix(tg.pairwise(X)), E, 10)
ref = trustworthiness(X, emb, n_neighbors=10)
check("per-point trustworthiness averages to scikit-learn's", abs(tw.mean() - ref) < 1e-12,
      f"{tw.mean():.6f} vs {ref:.6f}")
check("rank matching passes the requested number without ties",
      checks.top_matched(np.arange(10.0), 3).sum() == 3 and checks.top_matched(np.arange(10.0), 0).sum() == 0)
check("rank matching keeps ties at the cutoff", checks.top_matched(np.array([3, 2, 2, 2, 1.0]), 2).sum() == 4)

# scDEED reliability, against a direct transcription of Methods step 4
pre, e2 = rng.normal(size=(60, 5)), rng.normal(size=(60, 2))
got = checks.scdeed_scores(pre, e2)
want = []
for i in range(60):
    dp = np.linalg.norm(pre - pre[i], axis=1)
    de = np.linalg.norm(e2 - e2[i], axis=1)
    dp[i] = de[i] = np.inf
    m = 30
    v1 = de[np.argsort(dp, kind="stable")[:m]]
    v2 = np.sort(de)[:m]
    want.append(np.corrcoef(v1, v2)[0, 1])
check("scDEED score matches its definition", np.allclose(got, want), f"max diff {np.abs(got - want).max():.1e}")
check("scDEED score is 1 when the map equals the data", np.allclose(checks.scdeed_scores(e2, e2), 1.0))
share = checks.scdeed_trustworthy(X, emb, 1, np.random.default_rng(0), X, lambda A: A).mean()
check("scDEED labels few points trustworthy when the data have no structure", share < 0.25, f"{share:.3f}")

# recurrence and bootstrap in limiting cases
check("a neighborhood recurs in its own map", checks.recurs(E, emb, 8, 20).all())
Xu = np.repeat(X[:, None, :], 8, axis=1)                  # identical trials: every resample is c
sup = checks.bootstrap_support(Xu, E, 1, 8, 20, np.random.default_rng(0), B=2)
check("bootstrap support is 1 when every resample equals the averages", (sup == 1).all())

# end to end on the THINGS layout with a known truth
checks.B_BOOT = 5
n, p = 200, 60
Z = rng.normal(size=(n, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
mu = 2.5 * Z @ rng.normal(size=(4, p)) / np.sqrt(4)
Xs = mu[:, None, :] + rng.normal(size=(n, 80, p))
rows = []
for split in range(2):
    rows += run.split_task("things", 1, split, feat=Xs, truth=mu, R_list=(8, 16))
ref_rows = []
for split in range(2):
    r, _ = RT.split_rows(1, split, {"full": Xs}, things.sessions(), truth={"full": mu},
                         cells=[c for c in RT.CELLS if c[0] == "primary"], R_list=(8, 16))
    ref_rows += r
key = lambda r: (r["split"], r["R"], r["seed"], r["Rh"])
old = {key(r): (r["n1"], r["n10"], r["dis"], r["t1"]) for r in ref_rows}
new = {key(r): (r["n1"], r["n10"], r["dis"], r["t1"]) for r in rows}
check("re-created maps reproduce the confirmatory plug-in counts exactly", old == new, f"{len(old)} rows")
check("counts are internally consistent",
      all(r[f"{X}_a"] <= r[f"{X}_pass"] * r["nb"] and r[f"{X}_b"] <= r["n10"] for r in rows for X in checks.NAMES)
      and all(r["N_a"] == r["n10"] and r["N_b"] == 0 for r in rows))
check("TW passes at least as many points as N, per map",
      all(r["TW_pass"] >= r["n1"] for r in rows))
check("joint passes with N are consistent",
      all(r[f"{X}_nx"] <= min(r[f"{X}_pass"], r["n1"]) for r in rows for X in checks.NAMES)
      and all(r["N_nx"] == r["n1"] for r in rows))
est = ana.estimates(rows, RT.RH)
ex = {}
for R in (8, 16):
    sub = [r for r in rows if r["Rh"] == "exact" and r["R"] == R]
    for X in checks.NAMES:
        a, ps = sum(r[f"{X}_a"] for r in sub), sum(r[f"{X}_pass"] for r in sub)
        ex[(X, R)] = a / ps if ps else float("nan")
        e = est[(1, R)][X]
        print(f"   R={R:2d} {X:5s} pass {e['PASS']:.3f}  exact FC {ex[(X, R)]:.3f}  "
              f"estimates {e['FC'][0]:.3f} / {e['FC'][1]:.3f}  catch {e['CATCH'][0]:.3f} / {e['CATCH'][1]:.3f}"
              f"  with N: pass {e['PASS_NX']:.3f} FC {e['FC_NX'][0]:.3f} / {e['FC_NX'][1]:.3f}")
for R in (8, 16):
    e = est[(1, R)]["N"]["FC"]
    check(f"R={R}: N's estimate lands within 0.10 of its exact false comfort",
          min(abs(e[0] - ex[("N", R)]), abs(e[1] - ex[("N", R)])) < 0.10)
check("scDEED labels most points trustworthy on strongly structured data", est[(1, 8)]["SCD"]["PASS"] > 0.5,
      f"{est[(1, 8)]['SCD']['PASS']:.3f}")
ok = [X for X in checks.NAMES if est[(1, 8)][X]["PASS"] >= 0.02]
close = [X for X in ok if min(abs(v - ex[(X, 8)]) for v in est[(1, 8)][X]["FC"]) < 0.15]
check("R=8: every check passing at least 2% of points lands within 0.15 of its exact FC",
      close == ok, f"{len(ok)} checks: " + ", ".join(ok))

# the Allen layout and the semi-synthetic check
Za = rng.normal(size=(118, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
mua = 0.8 * Za @ rng.normal(size=(4, 150)) / 2.0
Xa = mua[:, None, :] + rng.normal(size=(118, 50, 150))
ra = run.split_task("allen", 1, 0, feat=Xa, truth=mua, R_list=(8,))
check("the Allen layout runs every check", {r["Rh"] for r in ra} == {8, 12, 24, "exact"}
      and all(f"{X}_pass" in ra[0] for X in checks.NAMES), f"{len(ra)} rows")
check("Allen halves hold R/4 repetitions per stratum",
      all((np.bincount(allen.strata()[allen.split_blocks(allen.strata(), np.random.default_rng(0))["U"][:8][h]],
                       minlength=2) == 2).all() for h in allen.halves(8)))
sim = run.sim_task("things", 1, feat=Xs, n_sim=1)
bias = ana.semi_synthetic_bias(sim, RT.RH)
check("semi-synthetic check and bias summary run",
      {r["Rh"] for r in sim} == {8, 20, 40, "exact"} and np.isfinite(bias["N"]["bias_lo"]),
      "FC bias of N " + f"{bias['N']['bias_lo']:+.3f} / {bias['N']['bias_hi']:+.3f} (1 draw)")
S = dict(FC={X: ana.summarize(est, 8, X, "FC", np.random.default_rng(0)) for X in checks.NAMES},
         guard_L={X: 0.0 for X in checks.NAMES}, guard_U={X: 0.0 for X in checks.NAMES})
cls = [ana.classify(S, X)[0] for X in ana.ALTERNATIVES]
check("classification runs (one subject: not evaluable)", all(c == "not evaluable" for c in cls))
print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
