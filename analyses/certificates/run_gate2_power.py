#!/usr/bin/env python3
"""Gate 2 power runs (plan Sections 5.5, 5.6 and 5.8).

Digits: gauss and gauss_het noise, gamma in {0.15, 0.30}, R in {8, 16, 32}, 5 draws each,
displays PCA-2 and t-SNE (perplexity 50, PCA init, seeds 17/29/31) fitted to the means.
Synthetic: the predecessor's four primary families, gauss, gamma in {0.15, 0.30}, R = 8,
PCA-2, 5 draws. For every engine and display: faithful / distorted shares and false labels
at eps x k', interval widths, naive-map mislabel rates, wall time, the ball witness ceiling
at k' = 40 (t-SNE displays, gamma = 0.15), and the planning-rule constants c_i.

    python3 run_gate2_power.py [--workers 32]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json, time, argparse, functools, platform
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn
from sklearn.manifold import TSNE

from reach import data, noise, target as tg, engines as eng
from reach.predecessor.panel import pca2
from reach.predecessor.intrinsic import levina_bickel
from reach.predecessor.witness import witness_ceiling
from reach.predecessor import generators as gen

K = 20
EPS = (0.1, 0.2, 0.3)
KPRIME = (20, 30, 40, 60, 80)
TSNE_SEEDS = (17, 29, 31)
TSNE_KW = dict(n_components=2, perplexity=50, init="pca", n_jobs=1)
SEED_BASE = 3_200_000
OUT = "results/gate2/power.jsonl"


def configs():
    C = []
    for model in ("gauss", "gauss_het"):
        for gamma in (0.15, 0.30):
            for R in (8, 16, 32):
                C.append(dict(geometry="digits", model=model, gamma=gamma, R=R, draws=5,
                              displays=["pca2"] + [f"tsne{s}" for s in TSNE_SEEDS]))
    for fam in gen.PRIMARY:
        for gamma in (0.15, 0.30):
            C.append(dict(geometry=fam, model="gauss", gamma=gamma, R=8, draws=5,
                          displays=["pca2"]))
    for i, c in enumerate(C):
        c["index"] = i
    return C


@functools.lru_cache(maxsize=None)
def panel(name):
    mu, d0 = data.geometry(name)
    d_mu = tg.pairwise(mu)
    idim = float(np.nanmedian(levina_bickel(d_mu, K)))
    d0_i = np.sort(d_mu, axis=1)[:, K - 1]
    return mu, d0, d_mu, tg.rank_matrix(d_mu), idim, d0_i


def embed(name, c):
    if name == "pca2":
        return pca2(c)
    return TSNE(random_state=int(name[4:]), **TSNE_KW).fit_transform(c)


def labels(lo, hi, truth, kp):
    f, d, t = hi <= kp, lo > kp, truth <= kp
    return dict(faithful=float(f.mean()), distorted=float(d.mean()),
                false_faithful=int((f & ~t).sum()), false_distorted=int((d & t).sum()))


def task(cfg, draw):
    mu, d0, d_mu, rank, idim, d0_i = panel(cfg["geometry"])
    n, D = mu.shape
    R = cfg["R"]
    rng = np.random.default_rng(SEED_BASE + 10_000 * cfg["index"] + draw)
    X = noise.acquire(mu, R, cfg["gamma"] * d0 / np.sqrt(D), cfg["model"], rng)
    c, s2, _ = eng.moments(X)
    boxes, times = {}, {}
    for name, fn in (("ball", lambda: eng.ball(X)), ("pair_t", lambda: eng.pair(X)),
                     ("pair_pool", lambda: eng.pair(X, pooled=True))):
        t0 = time.time()
        boxes[name] = fn()
        times[name] = time.time() - t0
    half = dict(ball=float(2 * np.median(boxes["ball"]["rho"])),
                pair_t=boxes["pair_t"]["half_width"], pair_pool=boxes["pair_pool"]["half_width"])
    naive = {"plugin": tg.naive_ranks(c, s2, R, D, corrected=False),
             "corrected": tg.naive_ranks(c, s2, R, D, corrected=True)}
    recs = []
    for dname in cfg["displays"]:
        t0 = time.time()
        E_idx, E_mask = tg.displayed(embed(dname, c), K)
        rec = dict(cfg=cfg["index"], draw=draw, display=dname, idim=idim,
                   tsne_seconds=time.time() - t0, engines={}, naive={}, truth={})
        bounds = {}
        for name, out in boxes.items():
            t0 = time.time()
            bounds[name] = {e: tg.box_reach(out["L"], out["U"], E_idx, E_mask, tg.m_of(e, K))
                            for e in EPS}
            times[f"{name}_reach"] = time.time() - t0
        t0 = time.time()
        rl, rh, _ = eng.cmp_member_ranks(X, E_idx)
        bounds["cmp"] = {e: tg.member_reach(rl, rh, tg.m_of(e, K)) for e in EPS}
        times["cmp"] = time.time() - t0
        truths = {e: tg.reach_from_ranks(rank, E_idx, tg.m_of(e, K)) for e in EPS}
        for e in EPS:
            rec["truth"][str(e)] = {str(kp): float((truths[e] <= kp).mean()) for kp in KPRIME}
            for nm, rk in naive.items():
                nr = tg.reach_from_ranks(rk, E_idx, tg.m_of(e, K))
                rec["naive"].setdefault(nm, {})[str(e)] = {
                    str(kp): dict(mislabel=float(((nr <= kp) != (truths[e] <= kp)).mean()),
                                  false_faithful=int(((nr <= kp) & (truths[e] > kp)).sum()),
                                  false_distorted=int(((nr > kp) & (truths[e] <= kp)).sum()))
                    for kp in KPRIME}
        for name, bd in bounds.items():
            er = dict(labels={}, width={})
            for e in EPS:
                lo, hi = bd[e]
                er["labels"][str(e)] = {str(kp): labels(lo, hi, truths[e], kp) for kp in KPRIME}
                er["width"][str(e)] = float(np.median(np.log2(hi / lo)))
            # planning rule (plan 5.8): distance-ratio inflation per unit half-width
            if name in half:
                lo, hi = bd[0.2]
                tr = truths[0.2]
                sel = (tr <= 4 * K) & (hi < n - 1)
                infl = (hi[sel] / K) ** (1 / idim) - (tr[sel] / K) ** (1 / idim)
                cc = infl * d0_i[sel] / half[name]
                er["planning_c_quantiles"] = [float(q) for q in np.percentile(cc, [10, 25, 50, 75, 90])]
                er["half_width_over_d0"] = half[name] / d0
            rec["engines"][name] = er
        if dname.startswith("tsne") and cfg["gamma"] == 0.15:
            b = boxes["ball"]
            dc = tg.pairwise(b["c"])
            v, _ = witness_ceiling(b["c"], b["rho"], dc, E_idx, E_mask, 2 * K,
                                   np.random.default_rng(SEED_BASE + 999))
            rec["ball_ceiling_k40_eps0.2"] = float((v >= tg.m_of(0.2, K)).mean())
        rec["seconds"] = dict(times)
        recs.append(rec)
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=32)
    a = ap.parse_args()
    C = configs()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(configs=C, seed_base=SEED_BASE, k=K, eps=EPS, kprime=KPRIME,
                   tsne=dict(seeds=TSNE_SEEDS, **TSNE_KW), delta=eng.DELTA, split=eng.SPLIT,
                   python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                   sklearn=sklearn.__version__, started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open("results/gate2/manifest_power.json", "w"), indent=1)
    todo = [(c, d) for c in C for d in range(c["draws"])]
    print(f"{len(todo)} tasks, {a.workers} workers", flush=True)
    t0, k = time.time(), 0
    with ProcessPoolExecutor(a.workers) as ex, open(OUT, "w") as f:
        futs = [ex.submit(task, c, d) for c, d in todo]
        for fu in as_completed(futs):
            for rec in fu.result():
                f.write(json.dumps(rec) + "\n")
            f.flush()
            k += 1
            if k % 20 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)} in {time.time() - t0:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
