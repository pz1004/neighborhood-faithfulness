#!/usr/bin/env python3
"""THINGS-EEG2 need pilot (plan Amendment 2). One process per subject.

Per subject and split: session-stratified blocks F (8 reps), U (32) and H (40); the
representation is fitted on F. The user's naive labels come from prefixes of U
(R = 8, 16, 32), and E-pair-t's resolved share is recorded for those and for U + H (R = 72).
For R <= 32 the script also records:
- disagreement with plug-in labels from disjoint sub-blocks of H (R_h = 8, 20, 40), for the
  extrapolated estimator;
- certified bounds from E-pair-t on all of H.
Split 0 also gets the semi-synthetic need and validity runs and the residual diagnostics.
Writes results/things_pilot/records.jsonl, semi_synthetic.jsonl and diagnostics.json;
analyze_things_need_pilot.py applies the decision rule.

    python3 run_things_need_pilot.py [--subjects 1,...,10] [--splits 10]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json, time, argparse, platform
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn
from sklearn.manifold import TSNE

from reach import things, target as tg, engines as eng
from reach.predecessor.panel import pca2

K, EPS = 10, 0.2
KPRIME = (10, 20, 40)
RU, RH = (8, 16, 32), (8, 20, 40)
TSNE_SEEDS = (17, 29, 31)
TSNE_KW = dict(n_components=2, perplexity=30, init="pca", n_jobs=1)
SPLIT_BASE, SIM_BASE = 3_500_000, 3_600_000
N_SIM_NEED, N_SIM_VALID = 30, 300
DATA, OUT = "data/things_eeg2", "results/things_pilot"


def displays(c):
    out = {"pca2": pca2(c)}
    for s in TSNE_SEEDS:
        out[f"tsne{s}"] = TSNE(random_state=s, **TSNE_KW).fit_transform(c)
    return out


def mean_rank(block):
    return tg.rank_matrix(tg.pairwise(block.mean(axis=1)))


def one_split(s, split, feat, ses):
    rng = np.random.default_rng(SPLIT_BASE + 1000 * s + split)
    b = things.split_blocks(ses, rng)
    W, m0, shrink = things.fit_representation(feat[:, b["F"]])
    Y = (feat - m0) @ W
    YU, YH = Y[:, b["U"]], Y[:, b["H"]]
    m = tg.m_of(EPS, K)
    ref = {Rh: [mean_rank(YH[:, j * Rh:(j + 1) * Rh]) for j in range(40 // Rh)] for Rh in RH}
    boxH = eng.pair(YH)
    recs = []
    for R in RU + (72,):
        Xu = YU[:, :R] if R <= 32 else np.concatenate([YU, YH], axis=1)
        c = Xu.mean(axis=1)
        urank = tg.rank_matrix(tg.pairwise(c))
        box = eng.pair(Xu)
        for dn, emb in displays(c).items():
            E_idx, E_mask = tg.displayed(emb, K)
            ureach = tg.reach_from_ranks(urank, E_idx, m)
            lo, hi = tg.box_reach(box["L"], box["U"], E_idx, E_mask, m)
            rec = dict(subject=s, split=split, R=R, display=dn, lw_shrinkage=shrink,
                       naive_faithful={str(kp): float((ureach <= kp).mean()) for kp in KPRIME},
                       use={str(kp): dict(faithful=float((hi <= kp).mean()),
                                          distorted=float((lo > kp).mean())) for kp in KPRIME})
            if R <= 32:
                dis = {str(kp): {} for kp in KPRIME}
                for Rh in RH:
                    rr = [tg.reach_from_ranks(rk, E_idx, m) for rk in ref[Rh]]
                    for kp in KPRIME:
                        dis[str(kp)][str(Rh)] = float(np.mean([((r <= kp) != (ureach <= kp)).mean()
                                                               for r in rr]))
                rec["disagree"] = dis
                hlo, hhi = tg.box_reach(boxH["L"], boxH["U"], E_idx, E_mask, m)
                cb = {}
                for kp in KPRIME:
                    cf, cd = hhi <= kp, hlo > kp
                    res = cf | cd
                    wrong = res & ((ureach <= kp) != cf)      # on resolved points truth = cf
                    cb[str(kp)] = dict(lo=float(wrong.mean()), hi=float(wrong.mean() + (~res).mean()),
                                       resolved=float(res.mean()))
                rec["certified_bounds"] = cb
            recs.append(rec)
    return recs, Y, b


def semi_synthetic(s, Y, b):
    """Truth = split 0's H means; noise = U residuals resampled within image (Amendment 2)."""
    YU, YH = Y[:, b["U"]], Y[:, b["H"]]
    mu = YH.mean(axis=1)
    d_t = tg.pairwise(mu)
    rank_t = tg.rank_matrix(d_t)
    nU = YU.shape[1]
    pool = (YU - YU.mean(axis=1, keepdims=True)) * np.sqrt(nU / (nU - 1))
    n = mu.shape[0]
    m = tg.m_of(EPS, K)

    def draw(rng, R):
        idx = rng.integers(0, nU, size=(n, R))
        return mu[:, None, :] + np.take_along_axis(pool, idx[:, :, None], axis=1)

    need = []
    for d in range(N_SIM_NEED):
        rng = np.random.default_rng(SIM_BASE + 10_000 * s + d)
        for R in RU:
            Xs = draw(rng, R)
            c = Xs.mean(axis=1)
            emb = TSNE(random_state=TSNE_SEEDS[d % 3], **TSNE_KW).fit_transform(c)
            E_idx, E_mask = tg.displayed(emb, K)
            tr = tg.reach_from_ranks(rank_t, E_idx, m)
            nr = tg.reach_from_ranks(tg.rank_matrix(tg.pairwise(c)), E_idx, m)
            box = eng.pair(Xs)
            lo, hi = tg.box_reach(box["L"], box["U"], E_idx, E_mask, m)
            rec = dict(subject=s, draw=d, R=R)
            for kp in KPRIME:
                rec[f"mislabel_{kp}"] = float(((nr <= kp) != (tr <= kp)).mean())
                rec[f"resolved_{kp}"] = float(((hi <= kp) | (lo > kp)).mean())
                rec[f"false_{kp}"] = int((((hi <= kp) & (tr > kp)) | ((lo > kp) & (tr <= kp))).sum())
            need.append(rec)
    valid = []
    iu = np.triu_indices(n, 1)
    for d in range(N_SIM_VALID):
        rng = np.random.default_rng(SIM_BASE + 10_000 * s + 1000 + d)
        Xs = draw(rng, 8)
        c = Xs.mean(axis=1)
        box = eng.pair(Xs)
        pv = int(((d_t[iu] < box["L"][iu]) | (d_t[iu] > box["U"][iu])).sum())
        rv = 0
        for E_idx, E_mask in (tg.displayed(pca2(c), K), tg.displayed(c, K)):
            for e in (0.1, 0.2, 0.3):
                mm = tg.m_of(e, K)
                lo, hi = tg.box_reach(box["L"], box["U"], E_idx, E_mask, mm)
                tr = tg.reach_from_ranks(rank_t, E_idx, mm)
                rv += int(((tr < lo) | (tr > hi)).sum())
        valid.append(dict(subject=s, draw=d, pair_viol=pv, reach_viol=rv))
    return need, valid


def diagnostics(s, Y, b, ses):
    """G3-diag on split 0's U + H residuals in the 20-dimensional space (Amendment 2, item 8)."""
    pos = np.concatenate([b["U"], b["H"]])
    Z, sz = Y[:, pos], ses[pos]
    n, R, D = Z.shape
    E = Z - Z.mean(axis=1, keepdims=True)
    lam = np.linalg.eigvalsh(np.einsum("nrd,nre->de", E, E) / (n * (R - 1)))
    flat = E.reshape(-1, D)
    kurt = (flat ** 4).mean(0) / (flat ** 2).mean(0) ** 2 - 3
    s2 = (E ** 2).sum(axis=(1, 2)) / (D * (R - 1))
    Ms = np.stack([Z[:, sz == k].mean(axis=1) for k in range(4)], axis=1)
    nk = R // 4
    inter = Ms - Ms.mean(1, keepdims=True) - Ms.mean(0, keepdims=True) + Ms.mean((0, 1), keepdims=True)
    within = sum(((Z[:, sz == k] - Ms[:, k:k + 1]) ** 2).sum() for k in range(4)) / (n * 4 * (nk - 1) * D)
    F_int = nk * (inter ** 2).sum() / ((n - 1) * 3 * D) / within
    F_ses = n * nk * ((Ms.mean(0) - Ms.mean((0, 1))) ** 2).sum() / (3 * D) / within
    return dict(subject=s, mean_residual_variance=float(lam.mean()),
                eig_max_over_min=float(lam.max() / lam.min()),
                eig_max_over_mean=float(lam.max() / lam.mean()),
                excess_kurtosis_median=float(np.median(kurt)), excess_kurtosis_max=float(kurt.max()),
                s2_cv_over_homoscedastic=float((s2.std() / s2.mean()) / np.sqrt(2 / (D * (R - 1)))),
                F_session_x_image=float(F_int), F_session=float(F_ses))


def subject(s, n_splits):
    t0 = time.time()
    X, ch, times = things.load_test(os.path.join(DATA, f"sub-{s:02d}", "preprocessed_eeg_test.npy"))
    assert X.shape == (200, 80, 17, 100), X.shape
    feat = things.features(X, times).astype(np.float64)
    del X
    ses = things.sessions()
    recs = []
    for split in range(n_splits):
        r, Y, b = one_split(s, split, feat, ses)
        recs += r
        if split == 0:
            Y0, b0 = Y, b
    need, valid = semi_synthetic(s, Y0, b0)
    diag = diagnostics(s, Y0, b0, ses)
    diag.update(channels=ch, n_times=int(len(times)), t_first=float(times[0]),
                t_last=float(times[-1]), seconds=time.time() - t0)
    return recs, need, valid, diag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(str(s) for s in range(1, 11)))
    ap.add_argument("--splits", type=int, default=10)
    a = ap.parse_args()
    subs = [int(x) for x in a.subjects.split(",")]
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(subjects=subs, splits=a.splits, k=K, eps=EPS, kprime=KPRIME, RU=RU, RH=RH,
                   blocks=things.BLOCKS, tsne=dict(seeds=TSNE_SEEDS, **TSNE_KW),
                   split_base=SPLIT_BASE, sim_base=SIM_BASE, n_sim_need=N_SIM_NEED,
                   n_sim_valid=N_SIM_VALID, split=eng.SPLIT, delta=eng.DELTA,
                   python=platform.python_version(), numpy=np.__version__,
                   scipy=scipy.__version__, sklearn=sklearn.__version__,
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    t0 = time.time()
    with ProcessPoolExecutor(len(subs)) as ex, \
            open(os.path.join(OUT, "records.jsonl"), "w") as fr, \
            open(os.path.join(OUT, "semi_synthetic.jsonl"), "w") as fs:
        futs = {ex.submit(subject, s, a.splits): s for s in subs}
        diags = []
        for fu in as_completed(futs):
            recs, need, valid, diag = fu.result()
            for r in recs:
                fr.write(json.dumps(r) + "\n")
            for r in need:
                fs.write(json.dumps(dict(kind="need", **r)) + "\n")
            for r in valid:
                fs.write(json.dumps(dict(kind="valid", **r)) + "\n")
            diags.append(diag)
            print(f"  sub-{futs[fu]:02d} done in {diag['seconds']:.0f}s", flush=True)
    json.dump(sorted(diags, key=lambda d: d["subject"]), open(os.path.join(OUT, "diagnostics.json"), "w"),
              indent=1)
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
