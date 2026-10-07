#!/usr/bin/env python3
"""Gate 2 validity runs (plan Sections 5.3 and 5.4): family-wise failure of every engine.

Each task is one noise draw of one configuration. For every engine it counts the false
statements the engine makes: every pairwise interval of the box engines (and every ball of
E-ball), every point's reach interval at eps in {0.1, 0.2, 0.3} under two displays (PCA-2 of
the means, and the noisy-kNN display that follows the noise exactly), and E-cmp's elementary
comparisons for the displayed members (all triples on digits200). A draw is a family-wise
failure for an engine if any count is nonzero. Results stream to a JSONL file, one line per
task, so an interrupted run resumes where it stopped.

    python3 run_gate2_validity.py [--workers 32] [--only CFG_INDEX,...]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time, argparse, functools, platform
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn

from reach import data, noise, target as tg, engines as eng
from reach.predecessor.panel import pca2

K, EPS = 20, (0.1, 0.2, 0.3)
SEED_BASE = 3_100_000
OUT = "results/gate2/validity.jsonl"


def configs():
    """Plan Section 5.4, in a fixed order: the index sets the seeds."""
    C = [dict(geometry="digits", model="gauss", R=8, gamma=0.15, draws=1000, primary=True),
         dict(geometry="digits", model="gauss", R=16, gamma=0.15, draws=300)]
    for model in ("gauss_het", "t5", "t3", "rep_t5", "aniso4", "aniso16"):
        for R in (8, 16):
            C.append(dict(geometry="digits", model=model, R=R, gamma=0.15, draws=300))
    for model in ("gauss", "gauss_het", "t5", "aniso4"):
        C.append(dict(geometry="surface", model=model, R=8, gamma=0.15, draws=300))
    for model in ("gauss", "t5"):
        C.append(dict(geometry="digits", model=model, R=8, gamma=0.5, draws=300))
    for model in ("gauss", "t5", "aniso4", "rep_t5"):
        C.append(dict(geometry="digits200", model=model, R=8, gamma=0.15, draws=300,
                      all_triples=True))
    for i, c in enumerate(C):
        c["index"] = i
    return C


@functools.lru_cache(maxsize=None)
def panel(name):
    mu, d0 = data.geometry(name)
    d_mu = tg.pairwise(mu)
    return mu, d0, d_mu, tg.rank_matrix(d_mu)


def reach_violations(lo_hi, rank, E_idx):
    """Points whose true reach falls outside its interval, summed over eps."""
    bad = 0
    for eps, (lo, hi) in lo_hi.items():
        truth = tg.reach_from_ranks(rank, E_idx, tg.m_of(eps, K))
        bad += int(((truth < lo) | (truth > hi)).sum())
    return bad


def task(cfg, draw):
    mu, d0, d_mu, rank = panel(cfg["geometry"])
    n, D = mu.shape
    rng = np.random.default_rng(SEED_BASE + 10_000 * cfg["index"] + draw)
    X = noise.acquire(mu, cfg["R"], cfg["gamma"] * d0 / np.sqrt(D), cfg["model"], rng)
    c = X.mean(1)
    disp = {"pca2": tg.displayed(pca2(c), K), "noisy_knn": tg.displayed(c, K)}
    iu = np.triu_indices(n, 1)
    dt = d_mu[iu]
    rec = dict(cfg=cfg["index"], draw=draw)

    for name, fn in (("ball", lambda: eng.ball(X)), ("pair_t", lambda: eng.pair(X)),
                     ("pair_pool", lambda: eng.pair(X, pooled=True))):
        t0 = time.time()
        out = fn()
        L, U = out["L"], out["U"]
        r = dict(pair_viol=int(((dt < L[iu]) | (dt > U[iu])).sum()), reach_viol=0)
        if name == "ball":
            r["ball_viol"] = int((np.linalg.norm(c - mu, axis=1) > out["rho"]).sum())
        for dn, (E_idx, E_mask) in disp.items():
            lh = {e: tg.box_reach(L, U, E_idx, E_mask, tg.m_of(e, K)) for e in EPS}
            r["reach_viol"] += reach_violations(lh, rank, E_idx)
            lo, hi = lh[0.2]
            r[f"resolved_{dn}"] = float(((hi <= 2 * K) | (lo > 2 * K)).mean())
        r["seconds"] = time.time() - t0
        rec[name] = r

    t0 = time.time()
    r = dict(elem_viol=0, reach_viol=0)
    for dn, (E_idx, E_mask) in disp.items():
        rl, rh, info = eng.cmp_member_ranks(X, E_idx, mu=mu)
        r["elem_viol"] += info["false_statements"]
        lh = {e: tg.member_reach(rl, rh, tg.m_of(e, K)) for e in EPS}
        r["reach_viol"] += reach_violations(lh, rank, E_idx)
        lo, hi = lh[0.2]
        r[f"resolved_{dn}"] = float(((hi <= 2 * K) | (lo > 2 * K)).mean())
    if cfg.get("all_triples"):
        J = np.array([[l for l in range(n) if l != i] for i in range(n)])
        _, _, info = eng.cmp_member_ranks(X, J, mu=mu)
        r["elem_viol_all_triples"] = info["false_statements"]
    r["seconds"] = time.time() - t0
    rec["cmp"] = r
    return rec


def done_tasks(path):
    seen = set()
    if os.path.exists(path):
        for line in open(path):
            if line.strip():
                o = json.loads(line)
                seen.add((o["cfg"], o["draw"]))
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    C = configs()
    only = {int(x) for x in a.only.split(",") if x}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(configs=C, seed_base=SEED_BASE, k=K, eps=EPS, delta=eng.DELTA, split=eng.SPLIT,
                   python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                   sklearn=sklearn.__version__, started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open("results/gate2/manifest_validity.json", "w"), indent=1)
    seen = done_tasks(OUT)
    todo = [(c, d) for c in C if not only or c["index"] in only
            for d in range(c["draws"]) if (c["index"], d) not in seen]
    print(f"{len(todo)} tasks to run ({len(seen)} already done), {a.workers} workers", flush=True)
    t0, n_done = time.time(), 0
    with ProcessPoolExecutor(a.workers) as ex, open(OUT, "a") as f:
        futs = [ex.submit(task, c, d) for c, d in todo]
        for fu in as_completed(futs):
            f.write(json.dumps(fu.result()) + "\n")
            f.flush()
            n_done += 1
            if n_done % 250 == 0 or n_done == len(todo):
                print(f"  {n_done}/{len(todo)} in {time.time() - t0:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
