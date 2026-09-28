#!/usr/bin/env python3
"""Scratch: how much rank tolerance buys, on the digits panel of pilot_comparison.py.
Reuses that file's definitions; t-SNE seed 42 only; gamma = 0.15."""
src = open("pilot_comparison.py").read().split("rows, t0 = [], time.time()")[0]
exec(src)
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
X_all = unc.acquire(mu, 0.15, d0, 64, np.random.default_rng(710_000))
TSNE_KW = dict(n_components=2, perplexity=50, init="pca", n_jobs=8)
t = 0.2; thr = THR[t]
print(f"digits, t-SNE seed 42, gamma=0.15, eps={t}: share of points")
print(" R   k'   truth  ballceil  cmp_faithful  cmp_distorted  unresolved")
for R in (8, 32):
    X = X_all[:, :R, :]
    c, rho, _ = unc.balls(X, DELTA, N)
    dc = cdist(c, c); np.fill_diagonal(dc, np.inf)
    E_idx, E_mask = displayed(TSNE(random_state=42, **TSNE_KW).fit_transform(c))
    for kp in (20, 25, 30, 40, 60, 80):
        o = truth(d_mu, E_idx, kp)
        v, _ = witness_ceiling(c, rho, dc, E_idx, E_mask, kp, np.random.default_rng(710_001))
        lo, hi, z = comparison_bounds(X, E_idx, kp, conservative=True)
        f, b = (lo >= thr).mean(), (hi < thr).mean()
        bad = int(((lo >= thr) & (o < thr)).sum() + ((hi < thr) & (o >= thr)).sum())
        print(f"{R:2d}  {kp:3d}   {(o >= thr).mean():.3f}   {(v >= thr).mean():.3f}     "
              f"{f:.3f}         {b:.3f}         {1 - f - b:.3f}   wrong={bad}", flush=True)
