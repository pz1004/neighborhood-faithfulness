#!/usr/bin/env python3
"""Same-panel comparison for docs/gate2_report.md, Section 5: on the synthetic power panels
(identical seeds to run_gate2_power.py), the predecessor's ball witness ceiling at the exact
target k' = k = 20, eps = 0.2, against E-pair-t's certified faithful share. Prints the table
saved as results/gate2/synthetic_ceiling.txt."""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from analyses.certificates import run_gate2_power as P
from reach import noise, target as tg, engines as eng
from reach.predecessor.panel import pca2
from reach.predecessor.witness import witness_ceiling

K = 20


def one(args):
    cfg, draw = args
    mu, d0, d_mu, rank, idim, d0_i = P.panel(cfg["geometry"])
    rng = np.random.default_rng(P.SEED_BASE + 10_000 * cfg["index"] + draw)
    X = noise.acquire(mu, cfg["R"], cfg["gamma"] * d0 / np.sqrt(mu.shape[1]), cfg["model"], rng)
    c = X.mean(1)
    E_idx, E_mask = tg.displayed(pca2(c), K)
    m = tg.m_of(0.2, K)
    b, p = eng.ball(X), eng.pair(X)
    v, _ = witness_ceiling(b["c"], b["rho"], tg.pairwise(b["c"]), E_idx, E_mask, K,
                           np.random.default_rng(P.SEED_BASE + 999))
    _, hi_b = tg.box_reach(b["L"], b["U"], E_idx, E_mask, m)
    _, hi_p = tg.box_reach(p["L"], p["U"], E_idx, E_mask, m)
    return np.array([(v >= m).mean(), (hi_b <= K).mean(), (hi_p <= K).mean()])


if __name__ == "__main__":
    allv = []
    for cfg in [c for c in P.configs() if c["geometry"] != "digits"]:
        with ProcessPoolExecutor(5) as ex:
            r = np.mean(list(ex.map(one, [(cfg, d) for d in range(cfg["draws"])])), axis=0)
        allv.append(r)
        print(f"{cfg['geometry']:18s} gamma={cfg['gamma']:.2f}: ball ceiling {r[0]:.3f}  ball interval {r[1]:.3f}"
              f"  pair_t {r[2]:.3f}", flush=True)
    a = np.mean(allv, axis=0)
    print(f"pooled over the 8 panels: ball ceiling {a[0]:.3f}  ball interval {a[1]:.3f}  pair_t {a[2]:.3f}")
