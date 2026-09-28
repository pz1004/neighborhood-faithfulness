#!/usr/bin/env python3
"""Structural tests of plan Section 5.7. Run before any Gate 2 driver; print one line per check.

    python3 tests/test_reach.py            # all checks, including the pilot reproduction
    python3 tests/test_reach.py --fast     # skip the reproduction (no t-SNE)
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np

from reach import target as tg, engines as eng, noise, data
from reach.predecessor.panel import pca2

FAILS = []


def check(name, cond, detail=""):
    print(f"[{'ok' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)
    if not cond:
        FAILS.append(name)


def explicit_reach(d, E_idx, m):
    """Reach by definition: sort others by (distance, id), take the m-th displayed position."""
    n = d.shape[0]
    out = np.empty(n, int)
    for i in range(n):
        others = [l for l in range(n) if l != i]
        order = sorted(others, key=lambda l: (d[i, l], l))
        pos = {l: p + 1 for p, l in enumerate(order)}
        out[i] = sorted(pos[j] for j in E_idx[i])[m - 1]
    return out


def test_definition(rng):
    n, k = 70, 12
    mu = rng.integers(0, 4, size=(n, 3)).astype(float)          # many exact ties
    d = tg.pairwise(mu)
    E_idx = np.array([rng.choice([l for l in range(n) if l != i], k, replace=False)
                      for i in range(n)])
    rank = tg.rank_matrix(d)
    for eps in (0.1, 0.2, 0.5):
        m = tg.m_of(eps, k)
        r = tg.reach_from_ranks(rank, E_idx, m)
        check(f"definition: reach equals explicit sort, eps={eps}",
              (r == explicit_reach(d, E_idx, m)).all())
        top = np.argsort(d, axis=1, kind="stable")
        ok = True
        for kp in (m, k, 2 * k, 3 * k):
            ov = np.array([len(set(E_idx[i]) & set(top[i, :kp])) for i in range(n)])
            ok &= ((r <= kp) == (ov >= m)).all()
        check(f"definition: reach <= k' iff overlap with Top_k' >= m, eps={eps}", ok)


def test_box(rng):
    n, k = 120, 10
    P = rng.normal(size=(n, 4))
    d = tg.pairwise(P)
    E_idx, E_mask = tg.displayed(P[:, :2] + 0.3 * rng.normal(size=(n, 2)), k)
    rank = tg.rank_matrix(d)
    for m in (6, 8, 10):
        truth = tg.reach_from_ranks(rank, E_idx, m)
        lo, hi = tg.box_reach(d, d, E_idx, E_mask, m)
        check(f"box: degenerate box is exact, m={m}", (lo == truth).all() and (hi == truth).all())
        a = rng.uniform(0, 0.5, size=d.shape)
        L, U = np.maximum(0, d - a), d + rng.uniform(0, 0.5, size=d.shape)
        lo, hi = tg.box_reach(L, U, E_idx, E_mask, m)
        check(f"box: contains truth, m={m}", ((lo <= truth) & (truth <= hi)).all())
        inside = True
        for _ in range(20):
            with np.errstate(invalid="ignore"):                  # inf - inf on the diagonal
                dd = L + rng.uniform(size=d.shape) * (U - L)
            np.fill_diagonal(dd, np.inf)
            r = tg.reach_from_ranks(tg.rank_matrix(dd), E_idx, m)
            inside &= ((lo <= r) & (r <= hi)).all()
        check(f"box: contains the reach of 20 random matrices in the box, m={m}", inside)
        adv = np.where(E_mask, U, L); np.fill_diagonal(adv, np.inf)
        fav = np.where(E_mask, L, U); np.fill_diagonal(fav, np.inf)
        check(f"box: members-at-U / outsiders-at-L attains r_hi, m={m}",
              (tg.reach_from_ranks(tg.rank_matrix(adv), E_idx, m) == hi).all())
        check(f"box: members-at-L / outsiders-at-U attains r_lo, m={m}",
              (tg.reach_from_ranks(tg.rank_matrix(fav), E_idx, m) == lo).all())


def test_zero_noise(rng):
    mu, d0 = data.geometry("surface")
    mu = mu[:300]
    X = np.repeat(mu[:, None, :], 8, axis=1)                    # identical replicates
    k = 20
    E_idx, E_mask = tg.displayed(pca2(mu), k)
    rank = tg.rank_matrix(tg.pairwise(mu))
    m = tg.m_of(0.2, k)
    truth = tg.reach_from_ranks(rank, E_idx, m)
    for name, out in (("ball", eng.ball(X)), ("pair_t", eng.pair(X)),
                      ("pair_pool", eng.pair(X, pooled=True))):
        lo, hi = tg.box_reach(out["L"], out["U"], E_idx, E_mask, m)
        check(f"zero noise: {name} reach bounds exact", (lo == truth).all() and (hi == truth).all())
    rl, rh, _ = eng.cmp_member_ranks(X, E_idx)
    lo, hi = tg.member_reach(rl, rh, m)
    check("zero noise: cmp reach bounds exact", (lo == truth).all() and (hi == truth).all())


def test_monotone(rng):
    mu, d0 = data.geometry("surface")
    mu = mu[:300]
    X = noise.acquire(mu, 8, 0.15 * d0 / np.sqrt(mu.shape[1]), "gauss", rng)
    a = eng.pair(X)
    b = eng.pair(X, dZ=0.004, dW=0.0005, dS=0.0005)
    check("monotone: pair_t intervals widen as delta falls",
          (b["L"] <= a["L"] + 1e-12).all() and (b["U"] >= a["U"] - 1e-12).all())
    c1, c2 = eng.ball(X), eng.ball(X, delta=0.005)
    check("monotone: ball intervals widen as delta falls",
          (c2["L"] <= c1["L"] + 1e-12).all() and (c2["U"] >= c1["U"] - 1e-12).all())


def test_reproduction():
    """Pilot of 2026-09-28 (pilot/2026-09-28/pilot_comparison.json), plan Section 5.7."""
    from sklearn.manifold import TSNE
    from reach.predecessor import uncertainty as unc, generators as gen
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rows = json.load(open(os.path.join(root, "pilot/2026-09-28/pilot_comparison.json")))
    K, N = 20, 1000
    m = tg.m_of(0.2, K)

    # synthetic: PCA-2 display, deterministic; the pilot scores a fixed 100-point subsample
    sub = np.arange(0, N, N // 100)[:100]
    for fi, fam in enumerate(gen.ALL):
        if fam == "gaussian_control":
            continue
        dev = gen.draw_means(fam, N, np.random.default_rng(900_000 + fi))
        d0 = unc.median_kth_neighbor_distance(dev, K)
        for gi, gamma in enumerate((0.15, 0.30)):
            seed = 100_000 + 1000 * fi + 100 * gi
            rng = np.random.default_rng(seed)
            mu = gen.draw_means(fam, N, rng)
            X = unc.acquire(mu, gamma, d0, 8, rng)
            E_idx, E_mask = tg.displayed(pca2(X.mean(1)), K)
            b = eng.ball(X)
            blo, bhi = tg.box_reach(b["L"], b["U"], E_idx, E_mask, m)
            rl, rh, _ = eng.cmp_member_ranks(X, E_idx)
            clo, chi_ = tg.member_reach(rl, rh, m)
            for kp in (20, 40):
                rec = [r for r in rows if r.get("part") == "synthetic" and r["family"] == fam
                       and r["gamma"] == gamma and r["kprime"] == kp][0]
                got = dict(ball_box_reliable=(bhi[sub] <= kp).mean(),
                           cmp_cons_reliable=(chi_[sub] <= kp).mean(),
                           cmp_cons_broken=(clo[sub] > kp).mean())
                same = all(abs(got[f] - rec[f + "@0.2"]) < 1e-12 for f in got)
                more = (blo[sub] > kp).mean() >= rec["ball_box_broken@0.2"] - 1e-12
                check(f"reproduction: synthetic {fam} gamma={gamma} k'={kp}", same and more,
                      " ".join(f"{f}={got[f]:.3f}/{rec[f + '@0.2']:.3f}" for f in got))

    # digits: t-SNE display at seed 42 (the pilot used n_jobs=8)
    mu, d0 = data.digits_panel()
    X_all = unc.acquire(mu, 0.15, d0, 64, np.random.default_rng(710_000))
    X = X_all[:, :8, :]
    c = X.mean(1)
    emb = TSNE(n_components=2, perplexity=50, init="pca", n_jobs=8, random_state=42).fit_transform(c)
    E_idx, E_mask = tg.displayed(emb, K)
    b = eng.ball(X)
    blo, bhi = tg.box_reach(b["L"], b["U"], E_idx, E_mask, m)
    rl, rh, _ = eng.cmp_member_ranks(X, E_idx)
    clo, chi_ = tg.member_reach(rl, rh, m)
    p = eng.pair(X)
    plo, phi = tg.box_reach(p["L"], p["U"], E_idx, E_mask, m)
    for kp in (20, 40):
        rec = [r for r in rows if r.get("part") == "digits" and r["R"] == 8 and r["seed"] == 42
               and r["kprime"] == kp][0]
        got = dict(ball_box_reliable=(bhi <= kp).mean(), cmp_cons_reliable=(chi_ <= kp).mean(),
                   cmp_cons_broken=(clo > kp).mean())
        same = all(abs(got[f] - rec[f + "@0.2"]) < 1e-12 for f in got)
        more = (blo > kp).mean() >= rec["ball_box_broken@0.2"] - 1e-12
        check(f"reproduction: digits t-SNE seed 42 R=8 k'={kp}", same and more,
              " ".join(f"{f}={got[f]:.3f}/{rec[f + '@0.2']:.3f}" for f in got))
        print(f"   seen before Gate 2 (plan 5.7): pair_t at this configuration, k'={kp}: "
              f"faithful {(phi <= kp).mean():.3f}, distorted {(plo > kp).mean():.3f}", flush=True)


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(3_090_000)
    test_definition(rng)
    test_box(rng)
    test_zero_noise(rng)
    test_monotone(rng)
    if "--fast" not in sys.argv:
        test_reproduction()
    print(f"\n{len(FAILS)} failures in {time.time() - t0:.0f}s" + (": " + ", ".join(FAILS) if FAILS else ""))
    sys.exit(1 if FAILS else 0)
