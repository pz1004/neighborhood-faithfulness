#!/usr/bin/env python3
"""Scratch pilot, not part of the repository: do comparison-level confidence statements
beat the ball ceiling, and what does a rank-tolerant target do?

Ball side: exactly the repository's balls, interval bound and witness ceiling, on the same
panels and seeds as run_radius_sweep.py (synthetic) and run_real_data_pilot.py (digits), so
the ball numbers reproduce the committed ones.

Comparison side: for query i, displayed member j and any other object l, the bias-corrected
statistic  Delta_ijl = |c_i-c_l|^2 - |c_i-c_j|^2 - (D/R)(s_l^2 - s_j^2)  is unbiased for
|mu_i-mu_l|^2 - |mu_i-mu_j|^2 and, to first order, Gaussian with variance
  4/R [ s_i^2 |mu_j-mu_l|^2 + s_l^2 |mu_i-mu_l|^2 + s_j^2 |mu_i-mu_j|^2 ]  + quadratic terms.
Conservative variant (primary): the variance uses simultaneous upper bounds on sigma^2 and
ball upper bounds on the distances, delta split three ways, and z is Bonferroni over every
unordered comparison n*C(n-1,2), so a data-dependent display costs nothing extra.
Optimistic variant: plain plug-ins and Bonferroni over the comparisons actually used.
Both rest on a Gaussian approximation in the far tail; this is a feasibility probe only.

Target: |E_i  cap  Top_k'(mu_i)|, with k' = k (the protocol) or k' = 2k (rank-tolerant).
Member j is certified IN if at most k'-1 others are not certified farther than j, and
certified OUT if at least k' others are certified closer. Lower bound = #IN, upper
bound = k - #OUT. Three-way map at tolerance eps: reliable (lower >= thr), broken
(upper < thr), unresolved otherwise.
"""
import sys, json, time
REPO = "/home/user/Workspaces/setwise-certificates-neighborhood-fidelity"
sys.path.insert(0, REPO)
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import norm, chi2
from sklearn.datasets import load_digits
from sklearn.manifold import TSNE

from setwise import uncertainty as unc, generators as gen
from setwise.panel import pca2
from setwise.certificates import interval_bounds
from setwise.witness import witness_ceiling, overlap_adverse

N, K, D, DELTA = 1000, 20, 20, 0.05
EPS = (0.10, 0.20, 0.30, 0.50)
THR = {e: int(np.ceil((1 - e) * K)) for e in EPS}
OUT = sys.argv[1] if len(sys.argv) > 1 else "pilot_comparison.json"


def comparison_bounds(X, E_idx, kp, conservative=True):
    n, R, Dm = X.shape
    c = X.mean(1)
    s2 = ((X - c[:, None, :]) ** 2).sum((1, 2)) / (Dm * (R - 1))
    df = Dm * (R - 1)
    dc = cdist(c, c)
    d2 = dc ** 2 - (Dm / R) * (s2[:, None] + s2[None, :])
    k = E_idx.shape[1]
    if conservative:
        s2v = s2 * df / chi2.ppf(DELTA / 3 / n, df)
        _, rho, _ = unc.balls(X, DELTA / 3, n)
        dv2 = (dc + rho[:, None] + rho[None, :]) ** 2
        z = norm.isf(DELTA / 3 / (2 * n * (n - 1) * (n - 2) / 2))
    else:
        s2v = s2
        dv2 = np.maximum(d2, 0.0)
        z = norm.isf(DELTA / (2 * n * k * (n - 2)))
    lo, hi = np.zeros(n, int), np.zeros(n, int)
    ar = np.arange(k)
    for i in range(n):
        J = E_idx[i]
        Dl = d2[i][None, :] - d2[i, J][:, None]                  # >0: l farther than j
        var = 4.0 / R * (s2v[i] * dv2[J, :] + s2v[None, :] * dv2[i][None, :]
                         + s2v[J][:, None] * dv2[i, J][:, None])
        var += (2 * Dm / R ** 2) * (1 + 1 / (R - 1)) * (s2v[None, :] ** 2 + s2v[J][:, None] ** 2)
        var += 4 * Dm * s2v[i] * (s2v[None, :] + s2v[J][:, None]) / R ** 2
        se = np.sqrt(var)
        ok = np.ones((k, n), bool); ok[:, i] = False; ok[ar, J] = False
        not_farther = ((Dl - z * se) <= 0) & ok
        closer = ((Dl + z * se) < 0) & ok
        lo[i] = int((not_farther.sum(1) <= kp - 1).sum())
        hi[i] = k - int((closer.sum(1) >= kp).sum())
    return lo, hi, float(z)


def ball_bounds(c, rho, dc, E_idx, E_mask, kp, rng):
    """Interval lower/upper bounds and the witness ceiling, for Top_k'."""
    L, U = interval_bounds(dc, rho)
    score = np.where(E_mask, U, L); np.fill_diagonal(score, np.inf)
    lo = overlap_adverse(score, E_mask, kp)
    hi = np.empty(len(c), int)
    for i in range(len(c)):
        J = E_idx[i]
        Ui = U[i].copy(); Ui[i] = np.inf
        closer = (Ui[None, :] < L[i, J][:, None]); closer[np.arange(len(J)), J] = False
        hi[i] = len(J) - int((closer.sum(1) >= kp).sum())
    v, _ = witness_ceiling(c, rho, dc, E_idx, E_mask, kp, rng)
    return lo, hi, v


def displayed(emb):
    de = cdist(emb, emb); np.fill_diagonal(de, np.inf)
    E_idx = np.argsort(de, axis=1, kind="stable")[:, :K]
    E_mask = np.zeros((len(emb), len(emb)), bool)
    np.put_along_axis(E_mask, E_idx, True, axis=1); np.fill_diagonal(E_mask, False)
    return E_idx, E_mask


def truth(d_mu, E_idx, kp):
    top = np.argsort(d_mu, axis=1, kind="stable")[:, :kp]
    return np.array([len(set(E_idx[i]) & set(top[i])) for i in range(len(E_idx))])


def summarize(tag, o, sub, **arms):
    rec = dict(tag)
    for e in EPS:
        t = THR[e]
        rec[f"truth@{e}"] = float((o[sub] >= t).mean())
        for name, (lo, hi) in arms.items():
            lo_, hi_ = lo[sub], hi[sub] if hi is not None else None
            rec[f"{name}_reliable@{e}"] = float((lo_ >= t).mean())
            rec[f"{name}_false_reliable@{e}"] = int(((lo_ >= t) & (o[sub] < t)).sum())
            if hi_ is not None:
                rec[f"{name}_broken@{e}"] = float((hi_ < t).mean())
                rec[f"{name}_false_broken@{e}"] = int(((hi_ < t) & (o[sub] >= t)).sum())
    for name, (lo, hi) in arms.items():
        rec[f"{name}_bound_violations"] = int((lo[sub] > o[sub]).sum()
                                             + (0 if hi is None else (hi[sub] < o[sub]).sum()))
    return rec


rows, t0 = [], time.time()

# ---------------------------------------------------------------- synthetic, as run_radius_sweep.py
R = 8
sub = np.arange(0, N, N // 100)[:100]
for fi, fam in enumerate(list(gen.ALL)):
    if fam == "gaussian_control" or fam not in ("folded_curve", "curved_surface",
                                                "thin_bridge", "imbalanced_mixture"):
        continue
    dev = gen.draw_means(fam, N, np.random.default_rng(900_000 + fi))
    d0 = unc.median_kth_neighbor_distance(dev, K)
    for gi, gamma in enumerate((0.15, 0.30)):
        seed = 100_000 + 1000 * fi + 100 * gi
        rng = np.random.default_rng(seed)
        mu = gen.draw_means(fam, N, rng)
        X = unc.acquire(mu, gamma, d0, R, rng)
        c, rho, _ = unc.balls(X, DELTA, N)
        d_mu = cdist(mu, mu); np.fill_diagonal(d_mu, np.inf)
        E_idx, E_mask = displayed(pca2(c))
        dc = cdist(c, c); np.fill_diagonal(dc, np.inf)
        for kp in (K, 2 * K):
            o = truth(d_mu, E_idx, kp)
            blo, bhi, v = ball_bounds(c, rho, dc, E_idx, E_mask, kp, np.random.default_rng(seed + 1))
            clo, chi_, z = comparison_bounds(X, E_idx, kp, conservative=True)
            olo, ohi, zo = comparison_bounds(X, E_idx, kp, conservative=False)
            rec = summarize(dict(part="synthetic", family=fam, gamma=gamma, R=R, kprime=kp,
                                 display="pca2", rho_over_d0=float(rho.mean() / d0),
                                 z_cons=z, z_opt=zo), o, sub,
                            ball_box=(blo, bhi), ball_ceiling=(v, None),
                            cmp_cons=(clo, chi_), cmp_opt=(olo, ohi))
            rows.append(rec)
            print(f"syn {fam:18s} g={gamma} k'={kp}  truth@.2={rec['truth@0.2']:.2f} "
                  f"ballceil@.2={rec['ball_ceiling_reliable@0.2']:.2f} "
                  f"cmp_cons@.2={rec['cmp_cons_reliable@0.2']:.2f} "
                  f"cmp_opt@.2={rec['cmp_opt_reliable@0.2']:.2f}  "
                  f"viol={rec['cmp_cons_bound_violations']}/{rec['cmp_opt_bound_violations']}",
                  flush=True)

# ---------------------------------------------------------------- digits, as run_real_data_pilot.py
X0 = load_digits().data
dev_idx, ev_idx = np.split(np.random.default_rng(20260923).permutation(len(X0)), [len(X0) - N])
keep = X0[dev_idx].std(0) > 0
Xd, Xe = X0[dev_idx][:, keep], X0[ev_idx][:, keep]
loc, scale = Xd.mean(0), Xd.std(0)
Zd = (Xd - loc) / scale
Vt = np.linalg.svd(Zd - Zd.mean(0), full_matrices=False)[2][:D]
mu = ((Xe - loc) / scale - Zd.mean(0)) @ Vt.T
d_mu = cdist(mu, mu); np.fill_diagonal(d_mu, np.inf)
d0 = unc.median_kth_neighbor_distance(mu, K)
allpts = np.arange(N)
TSNE_KW = dict(n_components=2, perplexity=50, init="pca", n_jobs=8)

for s in (42, 123, 456):                        # the display's own ceiling, rank-tolerant
    E_idx, _ = displayed(TSNE(random_state=s, **TSNE_KW).fit_transform(mu))
    for kp in (20, 30, 40):
        o = truth(d_mu, E_idx, kp)
        rec = dict(part="digits_exact", display="tsne2", seed=s, kprime=kp,
                   **{f"truth@{e}": float((o >= THR[e]).mean()) for e in EPS})
        rows.append(rec)
        print(f"digits exact tsne s={s} k'={kp}  " +
              "  ".join(f"truth@{e}={rec[f'truth@{e}']:.3f}" for e in EPS), flush=True)

gamma, seed = 0.15, 710_000
X_all = unc.acquire(mu, gamma, d0, 64, np.random.default_rng(seed))
for R in (8, 16, 32, 64):
    X = X_all[:, :R, :]
    c, rho, _ = unc.balls(X, DELTA, N)
    dc = cdist(c, c); np.fill_diagonal(dc, np.inf)
    for s in (42, 123, 456):
        E_idx, E_mask = displayed(TSNE(random_state=s, **TSNE_KW).fit_transform(c))
        for kp in (K, 2 * K):
            o = truth(d_mu, E_idx, kp)
            blo, bhi, v = ball_bounds(c, rho, dc, E_idx, E_mask, kp, np.random.default_rng(seed + 1))
            clo, chi_, z = comparison_bounds(X, E_idx, kp, conservative=True)
            olo, ohi, zo = comparison_bounds(X, E_idx, kp, conservative=False)
            rec = summarize(dict(part="digits", display="tsne2", seed=s, gamma=gamma, R=R,
                                 kprime=kp, rho_over_d0=float(rho.mean() / d0),
                                 z_cons=z, z_opt=zo), o, allpts,
                            ball_box=(blo, bhi), ball_ceiling=(v, None),
                            cmp_cons=(clo, chi_), cmp_opt=(olo, ohi))
            rows.append(rec)
            print(f"digits R={R:2d} s={s} k'={kp}  truth@.2={rec['truth@0.2']:.3f} "
                  f"ballceil@.2={rec['ball_ceiling_reliable@0.2']:.3f} "
                  f"cmp_cons@.2={rec['cmp_cons_reliable@0.2']:.3f} "
                  f"broken_cons@.2={rec['cmp_cons_broken@0.2']:.3f}  "
                  f"truth@.3={rec['truth@0.3']:.3f} cmp_cons@.3={rec['cmp_cons_reliable@0.3']:.3f}  "
                  f"viol={rec['cmp_cons_bound_violations']}/{rec['cmp_opt_bound_violations']}",
                  flush=True)

json.dump(rows, open(OUT, "w"), indent=1)
print(f"\n{len(rows)} rows in {time.time() - t0:.0f}s -> {OUT}")
