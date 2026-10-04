#!/usr/bin/env python3
"""Checks of the simulation-study code (docs/revision_protocol.md, Section 3) on synthetic data.

No data file is opened. Both layouts are checked, each with a synthetic base: a low-dimensional
signal in 60 features, unit noise. It checks the following:
- the truth families: T1 is block H's mean, T2 has rank 5 about its mean, T3 runs between them;
- the residual pool and the drawn means have the variances they should;
- mechanism D shares a component of correlation rho within each image and stratum;
- the zero-noise cell gives naive labels equal to the truth against every reference;
- the stress cell's analyst is more often right than the matched cell's;
- deeper references disagree less, and the counts keep their identities;
- the extrapolated estimates land near the exact value in an easy regime;
- development and validation seeds differ, and a validation run refuses to start before the
  addendum is committed.
Two splits and two regenerations per repetition keep the run short.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

from reach import sim, things, target as tg
import run_validation as RV
import analyze_validation as AV

fails = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        fails.append(name)


rng = np.random.default_rng(3_090_004)


def synthetic(n, reps, p=60, scale=2.5):
    Z = rng.normal(size=(n, 4)) * np.array([1.0, 0.8, 0.6, 0.4])
    mu = scale * Z @ rng.normal(size=(4, p)) / 2.0
    return mu[:, None, :] + rng.normal(size=(n, reps, p))


def cell(name, **kw):
    c = dict(name=name, family="T1", map="tsne", R_list=(8,), lam_F=1.0, lam_U=1.0, lam_H=1.0, rho=0.0)
    c.update(kw)
    return c


def rows_of(base, c, mu, seed, n_splits=2, n_regen=2):
    out = {}
    for r in sim.run_repetition(base, c, mu, seed, n_splits=n_splits, n_regen=n_regen):
        out.setdefault((r["layer"], r["regen"], r["R"]), {})[(r["ref"], r["Rh"])] = np.array(
            [r[f] for f in sim.FIELDS], float)
    return out


for ds, (n, reps) in (("EEG", (200, 80)), ("mouse", (118, 50))):
    feat = synthetic(n, reps)
    base = sim.Base(ds, 1, feat=feat)
    L = base.L
    b0 = L["split"](base.strata, np.random.default_rng(L["split_base"] + 1000))
    check(f"{ds}: T1 is block H's mean in split 0", np.allclose(sim.truth(base, "T1"), feat[:, b0["H"]].mean(axis=1)))
    mu2 = sim.truth(base, "T2")
    check(f"{ds}: T2 has rank 5 about its mean", np.linalg.matrix_rank(mu2 - mu2.mean(axis=0), tol=1e-8) == 5)
    check(f"{ds}: T3 is T2 at alpha 0 and T1 at alpha 1",
          np.allclose(sim.truth(base, "T3", alpha=0.0), mu2) and np.allclose(sim.truth(base, "T3", alpha=1.0), base.mu1))
    U = feat[:, b0["U"]]
    nU = U.shape[1]
    check(f"{ds}: the pool is centered within image and rescaled by sqrt(nU/(nU-1))",
          np.abs(base.pool.mean(axis=1)).max() < 1e-12
          and np.isclose(base.pool.var(axis=1).mean(), U.var(axis=1, ddof=1).mean()))
    draws = np.array([sim.mean_noise(base, rng, 8)[0] for _ in range(4000)])
    ratio = draws.var(axis=0).mean() / (base.pool[0].var(axis=0).mean() / 8)
    check(f"{ds}: a mean of 8 drawn residuals has the pool's variance over 8", abs(ratio - 1) < 0.08, f"ratio {ratio:.3f}")
    g = [sim.relative_gap(base, sim.truth(base, "T3", alpha=a)) for a in (0.0, 0.5, 1.0)]
    srt = np.sort(tg.pairwise((base.mu1 - base.m00) @ base.W0), axis=1)
    check(f"{ds}: the relative gap is computed as in explore_*_snr.py",
          np.isclose(g[2], np.median((srt[:, 19] - srt[:, 7]) / srt[:, 7])), f"gaps at alpha 0, 0.5, 1: {np.round(g, 3)}")

    # mechanism D: a shared component of correlation rho within image and stratum
    rep = sim.Repetition(base, cell("D", rho=0.3), base.mu1, np.random.default_rng(1))
    st = base.strata
    a_, b_ = np.flatnonzero(st == 0)[:2]
    c_ = np.flatnonzero(st == 1)[0]
    E = rep.E
    corr = lambda i, j: np.mean([np.corrcoef(E[k, i], E[k, j])[0, 1] for k in range(base.n)])
    check(f"{ds}: mechanism D correlates trials of a stratum by rho, not across strata",
          abs(corr(a_, b_) - 0.3) < 0.05 and abs(corr(a_, c_)) < 0.05,
          f"within {corr(a_, b_):.3f}, across {corr(a_, c_):.3f}")

    # the zero-noise cell
    z = rows_of(base, cell("zero", lam_U=0.0, lam_H=0.0), base.mu1, (1, 0))
    allz = all(v[1] == 0 and v[2] == 0 for blk in z.values() for v in blk.values())
    check(f"{ds}: zero noise in U and H gives naive labels equal to the truth against every reference", allz,
          f"{sum(len(b) for b in z.values())} pooled counts")

    # counts: identities, maps pooled, deeper references disagree less
    r1 = rows_of(base, cell("l1"), base.mu1, (2, 0))
    B = r1[("B", -1, 8)]
    ident = all(v[3] == v[0] - 2 * v[1] + v[2] for v in B.values())
    check(f"{ds}: pooled counts keep t1 = n1*nb - 2*n10 + dis, over 2 splits x 3 seeds",
          ident and B[("exact", 0)][4] == 6 * n, f"{len(B)} references")
    d = lambda ref, h: B[(ref, h)][1] / B[(ref, h)][0]
    ex = d("exact", 0)
    deep = L["deep"][-1]
    check(f"{ds}: the deepest reference is nearer the exact false comfort than the shallowest",
          abs(d("deep", deep) - ex) < abs(d("plug", L["RH"][0]) - ex),
          f"exact {ex:.3f}, R_h={L['RH'][0]} {d('plug', L['RH'][0]):.3f}, R_h={deep} {d('deep', deep):.3f}")
    rec = AV.record(B, L)
    check(f"{ds}: the estimates land near the exact false comfort in an easy regime",
          all(abs(rec[e] - ex) < 0.1 for e in ("sqrt", "inv", "b", "c")),
          ", ".join(f"{e} {rec[e]:.3f}" for e in ("sqrt", "inv", "b", "c")) + f"; exact {ex:.3f}")
    A = [k for k in r1 if k[0] == "A"]
    check(f"{ds}: Layer A has its fixed maps' exact counts and each regeneration's references",
          sorted(k[1] for k in A) == [-1, 0, 1] and ("deep", deep) in r1[("A", 0, 8)])

    # the stress cell: the analyst's block nearly noise-free
    rs = rows_of(base, cell("S", lam_U=0.1), base.mu1, (3, 0), n_regen=0)
    exS = rs[("B", -1, 8)][("exact", 0)]
    check(f"{ds}: the stress cell's analyst is more often right than the matched cell's",
          exS[1] / exS[0] < ex or exS[1] == 0, f"exact FC {exS[1] / max(exS[0], 1):.3f} vs {ex:.3f}")

    # cv-MDS maps run, with their own naive labels
    rm = rows_of(base, cell("cvmds", map="mds_cv"), base.mu1, (4, 0), n_regen=0)
    check(f"{ds}: cv-MDS cells run", ("exact", 0) in rm[("B", -1, 8)])

# seeds and the validation guard
c0 = RV.cells("EEG")[0]
check("development and validation seeds differ", RV.cell_seed("dev", "EEG", c0, 1, 0)[0] != RV.cell_seed("validation", "EEG", c0, 1, 0)[0])
check("cell seeds do not depend on the cells' order",
      RV.cell_seed("dev", "EEG", c0, 1, 0) == RV.cell_seed("dev", "EEG", dict(c0), 1, 0))
try:
    RV.require_committed(["docs/revision_addendum.md"])
    guarded = os.path.exists("docs/revision_addendum.md")
except SystemExit:
    guarded = True
check("a validation run refuses to start before the addendum is committed", guarded)
A, B_ = sim.halves_of(8, 4)
check("the halves of an EEG reference sub-block of 8 are balanced over sessions",
      sorted(A % 4) == [0, 1, 2, 3] and sorted(B_ % 4) == [0, 1, 2, 3])
check("there are 20 cells per dataset", len(RV.cells("EEG")) == 20 and len(RV.cells("mouse")) == 20)

print(f"\n{len(fails)} failures")
sys.exit(1 if fails else 0)
