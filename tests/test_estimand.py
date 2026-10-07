#!/usr/bin/env python3
"""Correctness tests of the false-comfort estimator in limiting cases (docs/revision_protocol.md,
Section 7).

No data file is opened. The synthetic cases use the THINGS layout. The count identities are also
checked on every committed result row. It checks the following:
- zero noise in the analyst's and reference blocks: naive labels equal the truth, every
  disagreement is 0, and both extrapolants return 0;
- exact references: the disagreement equals the exact false comfort at every reference size,
  and both extrapolants return it;
- zero acceptance: false comfort is NaN, not 0, and too few subjects give "not evaluable";
- full acceptance: a check that passes everything has FC = 1 - TF;
- ties are broken by object index, deterministically;
- the count identities, within each map and within each subject;
- weighting: proportionally more sub-blocks, or more splits with identical counts, change nothing.
"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

from reach import things, checks, target as tg
from analyses.eeg_confirmatory import run_things_confirmatory as RT
from analyses.mouse_v1 import analyze_allen_confirmatory as AA      # disagreement() and estimates() with RH as an argument
from analyses.checks import analyze_checks as AC

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_090_003)
RH = RT.RH
PRIMARY = [c for c in RT.CELLS if c[0] == "primary"]

# a THINGS-layout signal: 200 images, 80 repetitions in 4 sessions, a 4-dimensional signal in 60 features
n, p = 200, 60
Z = rng.normal(size=(n, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
mu = 2.5 * Z @ rng.normal(size=(4, p)) / np.sqrt(4)
ses = things.sessions()
s, split = 1, 0
b = things.split_blocks(ses, np.random.default_rng(RT.SPLIT_BASE + 1000 * s + split))   # as split_rows draws it


def rows_for(X):
    rows, _ = RT.split_rows(s, split, {"full": X, "early": X[:, :, :40]}, ses,
                            truth={"full": mu, "early": mu[:, :40]}, cells=PRIMARY, R_list=(8,), n_seeds=1)
    return rows


# zero noise in U and H (F keeps its noise, so the whitening can be fitted)
X = mu[:, None, :] + rng.normal(size=(n, 80, p))
X[:, b["U"]] = mu[:, None, :]
X[:, b["H"]] = mu[:, None, :]
rows = rows_for(X)
check("zero noise: naive labels equal the truth at every reference size and against the truth",
      all(r["n10"] == 0 and r["dis"] == 0 for r in rows), f"{len(rows)} rows")
est = AA.estimates(AA.disagreement([r for r in rows if r["Rh"] != "exact"], RH), RH)[(s, "primary", 8)]
check("zero noise: both extrapolants return 0 for FC and q",
      max(abs(v) for v in est["FC"] + est["q"]) < 1e-12, f"FC {est['FC']}, q {est['q']}")

# exact references: U noisy, H noise-free
X = mu[:, None, :] + rng.normal(size=(n, 80, p))
X[:, b["H"]] = mu[:, None, :]
rows = rows_for(X)
ex = [r for r in rows if r["Rh"] == "exact"][0]
same = all(r["n10"] == ex["n10"] * r["nb"] and r["dis"] == ex["dis"] * r["nb"] and r["t1"] == ex["t1"] * r["nb"]
           for r in rows if r["Rh"] != "exact")
check("exact references: every sub-block's counts equal the exact counts", same,
      f"exact FC {ex['n10'] / ex['n1']:.3f} from {ex['n1']} faithful calls")
est = AA.estimates(AA.disagreement([r for r in rows if r["Rh"] != "exact"], RH), RH)[(s, "primary", 8)]
check("exact references: both extrapolants return the exact false comfort",
      max(abs(v - ex["n10"] / ex["n1"]) for v in est["FC"]) < 1e-12, f"estimates {est['FC']}")


# synthetic check rows from label arrays, through the analysis code's own counting
def map_rows(subject, labels, faithful, kp=20):
    """Rows of one map: labels {check: bool[n]}, faithful {Rh: [bool[n] per sub-block]}."""
    out = []
    for Rh, blocks in faithful.items():
        rr = [np.where(T, 1, kp + 1) for T in blocks]             # reach 1 if faithful, k' + 1 if not
        out.append(dict(subject=subject, split=0, R=8, seed=0, Rh=Rh, **checks.check_counts(labels, rr, kp)))
    return out


def random_map(subject, pass_share, tf_by_rh, zero=(), full=()):
    N = rng.random(n) < pass_share
    labels = {X: (np.zeros(n, bool) if X in zero else np.ones(n, bool) if X in full
                  else (N if X == "N" else rng.random(n) < pass_share)) for X in checks.NAMES}
    faithful = {Rh: [rng.random(n) < tf for _ in range(40 // Rh)] for Rh, tf in tf_by_rh.items()}
    return map_rows(subject, labels, faithful)


tf_curve = {8: 0.06, 20: 0.09, 40: 0.11}
syn = [r for subj in range(1, 11) for _ in range(3)
       for r in random_map(subj, 0.2, tf_curve, zero=("BOOT",) + (("SH",) if subj <= 6 else ()), full=("TW",))]
est = AC.estimates(syn, RH)

# zero acceptance
fc_boot = [est[(subj, 8)]["BOOT"]["FC"] for subj in range(1, 11)]
check("zero acceptance: false comfort is NaN, not 0", all(np.isnan(v).all() for v in fc_boot))
rng_b = np.random.default_rng(0)
S = dict(FC={X: AC.summarize(est, 8, X, "FC", rng_b) for X in checks.NAMES},
         guard_L={X: 0.0 for X in checks.NAMES}, guard_U={X: 0.0 for X in checks.NAMES})
check("zero acceptance: 4 of 10 subjects with a pass give 'not evaluable'",
      AC.classify(S, "SH")[0] == "not evaluable" and S["FC"]["SH"]["n_defined"] == 4,
      f"SH defined in {S['FC']['SH']['n_defined']} subjects")
plug = AA.estimates(AA.disagreement([dict(r, cell="primary", n1=0, n10=0)
                                     for r in syn if r["subject"] == 1], RH), RH)[(1, "primary", 8)]
check("zero acceptance of the plug-in check: its false comfort is NaN", np.isnan(plug["FC"]).all())

# full acceptance
tf = AA.estimates(AA.disagreement([dict(r, cell="primary") for r in syn], RH), RH)
worst = max(abs(est[(subj, 8)]["TW"]["FC"][j] - (1 - tf[(subj, "primary", 8)]["TF"][j]))
            for subj in range(1, 11) for j in range(2))
check("full acceptance: FC = 1 - TF for a check that passes everything", worst < 1e-12, f"max |diff| {worst:.1e}")

# ties
D = np.array([[0, 0], [0, 0], [1, 0], [3, 0], [3, 0]], float)      # objects 0 and 1 coincide, as do 3 and 4
rk = tg.rank_matrix(tg.pairwise(D))
check("ties: equal distances are ranked by object index",
      rk[2, 0] < rk[2, 1] and rk[2, 3] < rk[2, 4] and rk[0, 1] == 1 and rk[0, 0] == len(D),
      f"row 2 ranks {rk[2].tolist()}")
E1, _ = tg.displayed(D, 2)
E2, _ = tg.displayed(D, 2)
check("ties: displayed neighbors are deterministic", (E1 == E2).all() and E1[2].tolist() == [0, 1])


# count identities
def identities(r):
    ok = r["t1"] == r["n1"] * r["nb"] - 2 * r["n10"] + r["dis"] and 0 <= r["n10"] <= r["n1"] * r["nb"]
    if "N_pass" in r:
        ok &= r["N_pass"] == r["n1"] and r["N_nx"] == r["n1"] and r["N_a"] == r["n10"] and r["N_b"] == 0
        for X in checks.NAMES:
            ok &= (r[f"{X}_a"] <= r[f"{X}_pass"] * r["nb"] and r[f"{X}_b"] <= r["n10"]
                   and r[f"{X}_nx"] <= min(r[f"{X}_pass"], r["n1"]))
    return ok


check("count identities hold in every synthetic map", all(identities(r) for r in syn), f"{len(syn)} rows")
files = sorted(glob.glob("results/*_confirmatory/*.jsonl") + glob.glob("results/checks/*/*.jsonl"))
files = [f for f in files if os.path.basename(f) in ("rows.jsonl", "semi_synthetic.jsonl")]
bad, total = 0, 0
for f in files:
    for line in open(f):
        total += 1
        bad += not identities(json.loads(line))
check("count identities hold in every committed row", bad == 0 and total > 0,
      f"{total} rows in {len(files)} files, {bad} violations")
worst = 0.0
for key, e in est.items():
    for X in checks.NAMES:
        for j in range(2):
            if np.isfinite(e[X]["FC"][j]):
                worst = max(worst, abs(e[X]["PASS"] - (e[X]["YIELD"][j] + e[X]["PASS"] * e[X]["FC"][j])))
check("within each subject, acceptance = faithful yield + unfaithful yield", worst < 1e-12,
      f"max |diff| {worst:.1e}")

# weighting
doubled = [dict(r, nb=2 * r["nb"], n10=2 * r["n10"], dis=2 * r["dis"], t1=2 * r["t1"], cell="primary") for r in syn]
a = AA.disagreement([dict(r, cell="primary") for r in syn], RH)
d2 = AA.disagreement(doubled, RH)
worst = max(abs(a[k][h][i] - d2[k][h][i]) for k in a for h in RH for i in range(3))
check("weighting: proportionally more sub-blocks leave every disagreement unchanged", worst < 1e-12,
      f"max |diff| {worst:.1e}")
nb_ok = True
for f in files:
    nbs = {}
    for line in open(f):
        r = json.loads(line)
        nbs.setdefault(r["Rh"], set()).add(r["nb"])
    nb_ok &= all(len(v) == 1 for v in nbs.values())
check("weighting: every map has the same number of sub-blocks at each reference size, in every committed file",
      nb_ok)
more = syn + [dict(r, split=1) for r in syn]                # a second split with identical counts
e1, e2 = AC.estimates(syn, RH), AC.estimates(more, RH)
same_est = all(np.allclose(e1[k][X]["FC"], e2[k][X]["FC"], equal_nan=True) for k in e1 for X in checks.NAMES)
L1 = AC.summarize(e1, 8, "N", "FC", np.random.default_rng(5))["L"]
L2 = AC.summarize(e2, 8, "N", "FC", np.random.default_rng(5))["L"]
check("weighting: more splits with identical counts change no estimate and no bound", same_est and L1 == L2,
      f"L {L1:.4f} vs {L2:.4f}")

print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
