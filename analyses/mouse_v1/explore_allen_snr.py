#!/usr/bin/env python3
"""Exploratory, post hoc (not part of docs/note_protocol_allen.md's decision): is there
image-specific structure in the Allen recordings, or is the true faithful share ~0 because
no neighborhood structure exists?

Split 0 of each subject, blocks U + H (40 repetitions) in the F-fitted 20-dim space:
- the share of image pairs with a positive cross-validated squared distance (two 20-rep halves);
- the split-half reliability of the distance matrix (cross-validated distances from two
  disjoint pairs of 10-rep quarters, correlated over pairs);
- the relative gap between each image's 8th and 20th true neighbors (the gap a faithful label
  at k=10, eps=0.2, k'=20 depends on), next to the pairwise noise of 16-rep means.
The same numbers are printed for THINGS-EEG2 alongside, from its own 72 U + H repetitions.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach import allen, things, target as tg
from analyses.mouse_v1 import run_allen_confirmatory as RA
from analyses.eeg_confirmatory import run_things_confirmatory as RT


def cv(A, B):
    """Cross-validated squared distances between image means of two independent blocks."""
    cA, cB = A.mean(1), B.mean(1)
    G = cA @ cB.T
    g = np.diag(G)
    return g[:, None] + g[None, :] - G - G.T


def stats(Z, R_user=16):
    n, R, D = Z.shape
    iu = np.triu_indices(n, 1)
    h = R // 2
    d2 = cv(Z[:, 0::2], Z[:, 1::2])                     # interleaved halves
    q = [Z[:, j::4] for j in range(4)]
    rel = np.corrcoef(cv(q[0], q[1])[iu], cv(q[2], q[3])[iu])[0, 1]
    d = np.sqrt(np.maximum(d2, 0.0))
    np.fill_diagonal(d, np.inf)
    srt = np.sort(d, axis=1)
    ok = srt[:, 7] > 0                                  # the gap is undefined where the 8th cv distance is 0
    E = Z - Z.mean(1, keepdims=True)
    s2 = (E ** 2).sum() / (n * D * (R - 1))
    noise = np.sqrt(2 * D * s2 / R_user)                # E|noise of a 16-rep mean pair|
    return dict(pos=float((d2[iu] > 0).mean()), reliability=float(rel),
                gap=float(np.median((srt[ok, 19] - srt[ok, 7]) / srt[ok, 7])), gap_n=int(ok.sum()),
                d8=float(np.median(srt[:, 7])), noise16=float(noise))


def allen_one(entry_s):
    s, entry = entry_s
    feats, _ = RA.verified_load(entry)
    b = allen.split_blocks(allen.strata(), np.random.default_rng(RA.SPLIT_BASE + 1000 * s + 0))
    W, m0 = RA.fit_variant("white", feats["full"][:, b["F"]])
    Y = (feats["full"] - m0) @ W
    return stats(Y[:, np.concatenate([b["U"], b["H"]])])


def things_one(s):
    X, ch, t = RT.verified_load(s)
    feat = things.features_window(X, t, 0.0, 0.79)
    del X
    b = things.split_blocks(things.sessions(), np.random.default_rng(RT.SPLIT_BASE + 1000 * s + 0))
    W, m0, _ = things.fit_representation(feat[:, b["F"]])
    Y = (feat - m0) @ W
    return stats(Y[:, np.concatenate([b["U"], b["H"]])])


if __name__ == "__main__":
    entries = RA.included()
    with ProcessPoolExecutor(10) as ex:
        A = list(ex.map(allen_one, [(i + 1, e) for i, e in enumerate(entries)]))
        T = list(ex.map(things_one, range(1, 11)))
    for name, rows in (("Allen two-photon, mouse V1 (40 reps)", A), ("THINGS-EEG2 (72 reps)", T)):
        med = lambda k: np.median([r[k] for r in rows])
        rng_ = lambda k: (min(r[k] for r in rows), max(r[k] for r in rows))
        print(f"{name}:")
        print(f"  image pairs with positive cross-validated d^2: {med('pos'):.3f} (range {rng_('pos')[0]:.3f}-{rng_('pos')[1]:.3f})")
        print(f"  split-half reliability of the distance matrix: {med('reliability'):.3f} (range {rng_('reliability')[0]:.3f}-{rng_('reliability')[1]:.3f})")
        print(f"  relative gap, 8th to 20th true neighbor: {med('gap'):.3f} (over images whose 8th cv distance is > 0: "
              f"median {med('gap_n'):.0f} of {118 if 'Allen' in name else 200}); 8th-neighbor distance {med('d8'):.3f}; "
              f"noise of a 16-rep mean pair {med('noise16'):.3f}; ratio {med('d8') / med('noise16'):.2f}")
