#!/usr/bin/env python3
"""Confirmatory analysis of docs/note_protocol.md (plan Amendment 3) on THINGS-EEG2.

One process per subject. Every data file is checked against its SHA-256 before it is
loaded. For each of 10 new splits, the 12 cells of protocol Section 3 are computed. Each
row holds, for one display and one reference size R_h, the counts behind the disagreement
statistics of protocol Section 4. The semi-synthetic estimator check of Section 5 adds rows
whose Rh is "exact", counted against a known truth. analyze_things_confirmatory.py applies
the decision rules of Section 6.

    python3 run_things_confirmatory.py [--subjects 1,...,10]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, hashlib, json, platform, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn
from sklearn.manifold import TSNE, MDS

from reach import things, target as tg
from reach.predecessor.panel import pca2

RU, RH = (8, 16, 32), (8, 20, 40)
N_SPLITS, N_SIM = 10, 20
SPLIT_BASE, SIM_BASE = 4_100_000, 4_200_000
TSNE_SEEDS = (4_000_101, 4_000_102, 4_000_103)
MDS_SEEDS = (4_000_201, 4_000_202, 4_000_203)
WINDOWS = {"full": (0.0, 0.79), "early": (0.05, 0.50)}
BASE = dict(pipe="P1", rep="white", win="full", disp="tsne30", k=10, eps=0.2)
CELLS = [(name, {**BASE, **change}) for name, change in (
    ("primary", {}),
    ("rsa_cv_mds", dict(pipe="P2", disp="mds_cv")),
    ("k5", dict(k=5)), ("k20", dict(k=20)),
    ("eps0.1", dict(eps=0.1)), ("eps0.3", dict(eps=0.3)),
    ("tsne_perp15", dict(disp="tsne15")), ("tsne_perp50", dict(disp="tsne50")),
    ("pca2", dict(disp="pca2")), ("mds", dict(disp="mds")),
    ("early_window", dict(win="early")), ("unwhitened", dict(rep="pca")))]
DATA, OUT = "data/things_eeg2", "results/things_confirmatory"


def seeds_for(disp):
    if disp.startswith("tsne"):
        return TSNE_SEEDS
    if disp.startswith("mds"):
        return MDS_SEEDS
    return (None,)


def _mds(Dm, seed):
    Dm = Dm.copy()
    np.fill_diagonal(Dm, 0.0)
    return MDS(n_components=2, dissimilarity="precomputed", n_init=1,
               random_state=seed).fit_transform(Dm)


def embed(disp, c, seed, d2cv=None):
    if disp.startswith("tsne"):
        return TSNE(n_components=2, perplexity=int(disp[4:]), init="pca", random_state=seed,
                    n_jobs=1).fit_transform(c)
    if disp == "pca2":
        return pca2(c)
    if disp == "mds":
        return _mds(tg.pairwise(c), seed)
    if disp == "mds_cv":
        return _mds(np.sqrt(np.maximum(d2cv, 0.0)), seed)
    raise KeyError(disp)


def counts(N, ref_reaches, kp):
    """Counts behind d_q, d_FC and d_TF for one display, summed over reference sub-blocks."""
    n10 = dis = t1 = 0
    for rr in ref_reaches:
        T = rr <= kp
        n10 += int((N & ~T).sum())
        dis += int((N != T).sum())
        t1 += int(T.sum())
    return dict(n=int(len(N)), n1=int(N.sum()), n10=n10, dis=dis, t1=t1, nb=len(ref_reaches))


def fit_variant(rep, featF):
    if rep == "white":
        W, m0, _ = things.fit_representation(featF)
    else:
        W, m0 = things.fit_pca_representation(featF)
    return W, m0


def split_rows(s, split, feats, ses, truth=None, cells=CELLS, R_list=RU, n_seeds=3):
    """All rows of one split. `truth` (synthetic tests only) maps a window to true means."""
    b = things.split_blocks(ses, np.random.default_rng(SPLIT_BASE + 1000 * s + split))
    rows = []
    for rep, win in sorted({(c["rep"], c["win"]) for _, c in cells}):
        feat = feats[win]
        W, m0 = fit_variant(rep, feat[:, b["F"]])
        Y = (feat - m0) @ W
        YU, YH = Y[:, b["U"]], Y[:, b["H"]]
        ref = {Rh: [tg.rank_matrix(tg.pairwise(YH[:, j * Rh:(j + 1) * Rh].mean(axis=1)))
                    for j in range(40 // Rh)] for Rh in RH}
        trank = None if truth is None else tg.rank_matrix(tg.pairwise((truth[win] - m0) @ W))
        vcells = [(nm, c) for nm, c in cells if (c["rep"], c["win"]) == (rep, win)]
        for R in R_list:
            Xu = YU[:, :R]
            cmean = Xu.mean(axis=1)
            nrank = {"P1": tg.rank_matrix(tg.pairwise(cmean))}
            d2cv = None
            if any(c["pipe"] == "P2" for _, c in vcells):
                d2cv = things.cv_sq_dist(Xu, *things.halves(R))
                nrank["P2"] = tg.rank_matrix(d2cv)
            for disp in sorted({c["disp"] for _, c in vcells}):
                for seed in seeds_for(disp)[:n_seeds]:
                    emb = embed(disp, cmean, seed, d2cv)
                    for nm, c in vcells:
                        if c["disp"] != disp:
                            continue
                        k, m, kp = c["k"], tg.m_of(c["eps"], c["k"]), 2 * c["k"]
                        E_idx, _ = tg.displayed(emb, k)
                        N = tg.reach_from_ranks(nrank[c["pipe"]], E_idx, m) <= kp
                        base = dict(subject=s, split=split, cell=nm, R=R, seed=seed)
                        for Rh in RH:
                            rr = [tg.reach_from_ranks(rk, E_idx, m) for rk in ref[Rh]]
                            rows.append(dict(base, Rh=Rh, **counts(N, rr, kp)))
                        if trank is not None:
                            rows.append(dict(base, Rh="exact",
                                             **counts(N, [tg.reach_from_ranks(trank, E_idx, m)], kp)))
    return rows, b


def semi_synthetic(s, feats, b):
    """Protocol Section 5: truth = H means; user and reference blocks resample U residuals."""
    feat = feats["full"]
    W, m0, _ = things.fit_representation(feat[:, b["F"]])
    Y = (feat - m0) @ W
    YU, YH = Y[:, b["U"]], Y[:, b["H"]]
    mu = YH.mean(axis=1)
    trank = tg.rank_matrix(tg.pairwise(mu))
    nU = YU.shape[1]
    pool = (YU - YU.mean(axis=1, keepdims=True)) * np.sqrt(nU / (nU - 1))
    n = mu.shape[0]
    k, m, kp = 10, tg.m_of(0.2, 10), 20
    rows = []
    for d in range(N_SIM):
        rng = np.random.default_rng(SIM_BASE + 10_000 * s + d)

        def draw(R):
            idx = rng.integers(0, nU, size=(n, R))
            return mu[:, None, :] + np.take_along_axis(pool, idx[:, :, None], axis=1)

        for R in RU:
            c = draw(R).mean(axis=1)
            emb = TSNE(n_components=2, perplexity=30, init="pca", random_state=TSNE_SEEDS[d % 3],
                       n_jobs=1).fit_transform(c)
            E_idx, _ = tg.displayed(emb, k)
            N = tg.reach_from_ranks(tg.rank_matrix(tg.pairwise(c)), E_idx, m) <= kp
            base = dict(subject=s, draw=d, R=R)
            rows.append(dict(base, Rh="exact", **counts(N, [tg.reach_from_ranks(trank, E_idx, m)], kp)))
            for Rh in RH:
                rr = [tg.reach_from_ranks(tg.rank_matrix(tg.pairwise(draw(Rh).mean(axis=1))), E_idx, m)
                      for _ in range(40 // Rh)]
                rows.append(dict(base, Rh=Rh, **counts(N, rr, kp)))
    return rows


def verified_load(s):
    man = {m["subject"]: m for m in json.load(open(os.path.join(DATA, "manifest.json")))}
    path = os.path.join(DATA, f"sub-{s:02d}", "preprocessed_eeg_test.npy")
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if h != man[s]["sha256"]:
        raise IOError(f"sub-{s:02d}: SHA-256 mismatch against the manifest")
    return things.load_test(path)


def subject(s):
    t0 = time.time()
    X, ch, times = verified_load(s)
    assert X.shape == (200, 80, 17, 100), X.shape
    feats = {w: things.features_window(X, times, *WINDOWS[w]).astype(np.float64) for w in WINDOWS}
    del X
    ses = things.sessions()
    rows, b0 = [], None
    for split in range(N_SPLITS):
        r, b = split_rows(s, split, feats, ses)
        rows += r
        b0 = b if split == 0 else b0
    sim = semi_synthetic(s, feats, b0)
    return rows, sim, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subjects", default=",".join(str(s) for s in range(1, 11)))
    a = ap.parse_args()
    subs = [int(x) for x in a.subjects.split(",")]
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(protocol="docs/note_protocol.md", subjects=subs, cells=CELLS, RU=RU, RH=RH,
                   n_splits=N_SPLITS, n_sim=N_SIM, split_base=SPLIT_BASE, sim_base=SIM_BASE,
                   tsne_seeds=TSNE_SEEDS, mds_seeds=MDS_SEEDS, windows=WINDOWS,
                   python=platform.python_version(), numpy=np.__version__,
                   scipy=scipy.__version__, sklearn=sklearn.__version__,
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    t0 = time.time()
    with ProcessPoolExecutor(len(subs)) as ex, \
            open(os.path.join(OUT, "rows.jsonl"), "w") as fr, \
            open(os.path.join(OUT, "semi_synthetic.jsonl"), "w") as fs:
        futs = {ex.submit(subject, s): s for s in subs}
        for fu in as_completed(futs):
            rows, sim, sec = fu.result()
            for r in rows:
                fr.write(json.dumps(r) + "\n")
            for r in sim:
                fs.write(json.dumps(r) + "\n")
            print(f"  sub-{futs[fu]:02d} done in {sec:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
