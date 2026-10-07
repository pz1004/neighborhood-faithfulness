#!/usr/bin/env python3
"""Need cross-tab for docs/gate2_report.md, Section 6: on points the naive map mislabels, does
E-pair-t correct the label or abstain? Same seeds and configurations as run_gate2_power.py
(digits, gauss, t-SNE displays), recomputed per point because the power runs saved aggregates
only. Prints the table saved as results/gate2/need_crosstab.txt."""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from analyses.certificates import run_gate2_power as P
from reach import noise, target as tg, engines as eng

K, KP, EPS = 20, 40, 0.2


def one(args):
    cfg, draw = args
    mu, d0, d_mu, rank, idim, d0_i = P.panel("digits")
    rng = np.random.default_rng(P.SEED_BASE + 10_000 * cfg["index"] + draw)
    X = noise.acquire(mu, cfg["R"], cfg["gamma"] * d0 / np.sqrt(mu.shape[1]), cfg["model"], rng)
    c, s2, _ = eng.moments(X)
    box = eng.pair(X)
    nrank = tg.naive_ranks(c, s2, cfg["R"], mu.shape[1], corrected=False)
    m = tg.m_of(EPS, K)
    out = np.zeros((2, 3), int)       # rows: naive right / wrong; cols: certified right / abstain / wrong
    for dn in cfg["displays"][1:]:
        E_idx, E_mask = tg.displayed(P.embed(dn, c), K)
        truth = tg.reach_from_ranks(rank, E_idx, m) <= KP
        naive = tg.reach_from_ranks(nrank, E_idx, m) <= KP
        lo, hi = tg.box_reach(box["L"], box["U"], E_idx, E_mask, m)
        cert = np.where(hi <= KP, 1, np.where(lo > KP, 0, -1))
        for w, sel in ((0, naive == truth), (1, naive != truth)):
            out[w, 0] += int(((cert == truth.astype(int)) & sel).sum())
            out[w, 1] += int(((cert == -1) & sel).sum())
            out[w, 2] += int(((cert >= 0) & (cert != truth.astype(int)) & sel).sum())
    return out


if __name__ == "__main__":
    for gamma in (0.15, 0.30):
        for R in (8, 16):
            cfg = [c for c in P.configs() if c["geometry"] == "digits" and c["model"] == "gauss"
                   and c["gamma"] == gamma and c["R"] == R][0]
            with ProcessPoolExecutor(5) as ex:
                tab = sum(ex.map(one, [(cfg, d) for d in range(cfg["draws"])]))
            nr, nw = tab[0].sum(), tab[1].sum()
            print(f"gamma={gamma:.2f} R={R:2d} (k'={KP}, eps={EPS}, 15 t-SNE displays): "
                  f"naive right on {nr} points -> pair_t right {tab[0,0]}, abstains {tab[0,1]}, wrong {tab[0,2]}; "
                  f"naive wrong on {nw} -> pair_t right {tab[1,0]}, abstains {tab[1,1]}, wrong {tab[1,2]}", flush=True)
