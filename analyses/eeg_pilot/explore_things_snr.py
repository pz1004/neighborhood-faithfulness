#!/usr/bin/env python3
"""Exploratory, post hoc (not part of Amendment 2's decision): why E-pair-t resolves nothing
on THINGS-EEG2. Split 0 of each subject, U + H (72 repetitions) in the F-fitted space.

Cross-validated squared distances (inner products of two independent 36-repetition halves)
are unbiased for the true squared distances. From them: how far apart the images truly
are relative to the noise, and how different the 8th and 20th true neighbors of a point
are in distance, which is the gap a certificate at k=10, eps=0.2, k'=20 must resolve.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach import things, engines as eng, target as tg
from analyses.eeg_pilot import run_things_need_pilot as P


def one(s):
    X, ch, times = things.load_test(f"{P.DATA}/sub-{s:02d}/preprocessed_eeg_test.npy")
    feat = things.features(X, times)
    del X
    ses = things.sessions()
    b = things.split_blocks(ses, np.random.default_rng(P.SPLIT_BASE + 1000 * s + 0))
    W, m0, _ = things.fit_representation(feat[:, b["F"]])
    Y = (feat - m0) @ W
    Z = Y[:, np.concatenate([b["U"], b["H"]])]                # 72 reps, round-robin by session
    A, B = Z[:, 0::2].mean(1), Z[:, 1::2].mean(1)             # two 36-rep halves, 9 per session
    dA, dB = A[:, None, :] - A[None], B[:, None, :] - B[None]
    d2cv = (dA * dB).sum(-1)                                   # unbiased for |mu_i - mu_j|^2
    n = len(A)
    iu = np.triu_indices(n, 1)
    E = Z - Z.mean(1, keepdims=True)
    s2 = (E ** 2).sum((1, 2)) / (Z.shape[2] * (Z.shape[1] - 1))
    noise_pair_d2 = Z.shape[2] * 2 * s2.mean() / Z.shape[1]   # E|noise of a 72-rep mean pair|^2
    d = np.sqrt(np.maximum(d2cv, 0))
    np.fill_diagonal(d, np.inf)
    srt = np.sort(d, axis=1)
    gap = np.median((srt[:, 19] - srt[:, 7]) / srt[:, 7])      # 20th vs 8th true neighbor
    box = eng.pair(Z)
    hw = np.median((box["U"] - box["L"])[iu] / 2)
    return dict(subject=s, frac_pairs_d2cv_pos=float((d2cv[iu] > 0).mean()),
                median_true_d=float(np.median(d[iu])), noise_pair_d=float(np.sqrt(noise_pair_d2)),
                median_d8=float(np.median(srt[:, 7])), rel_gap_8_to_20=float(gap),
                pair_t_half_width=float(hw))


if __name__ == "__main__":
    with ProcessPoolExecutor(10) as ex:
        rows = list(ex.map(one, range(1, 11)))
    for r in rows:
        print(r)
    med = lambda k: np.median([r[k] for r in rows])
    print(f"\nmedians over subjects: true distance {med('median_true_d'):.3f}, noise of a 72-rep mean"
          f" pair {med('noise_pair_d'):.3f}; 8th-neighbor distance {med('median_d8'):.3f}, relative gap"
          f" to the 20th {med('rel_gap_8_to_20'):.3f}; E-pair-t half-width at R=72 {med('pair_t_half_width'):.3f};"
          f" pairs with positive cross-validated d^2 {med('frac_pairs_d2cv_pos'):.3f}")


def false_comfort(s):
    """Exploratory: in split 0's semi-synthetic setting (truth = H means, resampled U noise),
    the true faithful share of the user's t-SNE display and the naive map's errors by kind."""
    from sklearn.manifold import TSNE
    X, ch, times = things.load_test(f"{P.DATA}/sub-{s:02d}/preprocessed_eeg_test.npy")
    feat = things.features(X, times)
    del X
    b = things.split_blocks(things.sessions(), np.random.default_rng(P.SPLIT_BASE + 1000 * s + 0))
    W, m0, _ = things.fit_representation(feat[:, b["F"]])
    Y = (feat - m0) @ W
    YU, YH = Y[:, b["U"]], Y[:, b["H"]]
    mu = YH.mean(1)
    rank_t = tg.rank_matrix(tg.pairwise(mu))
    pool = (YU - YU.mean(1, keepdims=True)) * np.sqrt(32 / 31)
    m, out = tg.m_of(0.2, 10), {}
    for R in (8, 16, 32):
        acc = np.zeros(4)
        for d in range(10):
            rng = np.random.default_rng(3_800_000 + 100 * s + d)
            idx = rng.integers(0, 32, size=(200, R))
            c = (mu[:, None, :] + np.take_along_axis(pool, idx[:, :, None], axis=1)).mean(1)
            E_idx, _ = tg.displayed(TSNE(random_state=17, **P.TSNE_KW).fit_transform(c), 10)
            t = tg.reach_from_ranks(rank_t, E_idx, m) <= 20
            nv = tg.reach_from_ranks(tg.rank_matrix(tg.pairwise(c)), E_idx, m) <= 20
            acc += [t.mean(), nv.mean(), (nv & ~t).mean(), (~nv & t).mean()]
        out[R] = acc / 10
    return out


if __name__ == "__main__":
    with ProcessPoolExecutor(10) as ex:
        fc = list(ex.map(false_comfort, range(1, 11)))
    print("\nsemi-synthetic, k'=20: truth faithful / naive faithful / naive false-faithful / naive false-distorted (medians)")
    for R in (8, 16, 32):
        v = np.median([f[R] for f in fc], axis=0)
        print(f"  R={R:2d}: {v[0]:.3f} / {v[1]:.3f} / {v[2]:.3f} / {v[3]:.3f}   "
              f"(share of naive 'faithful' calls that are wrong: {v[2] / v[1]:.2f})")
