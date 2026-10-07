#!/usr/bin/env python3
"""Confirmatory analysis of docs/note_protocol_allen.md (plan Amendment 4), second dataset.

Subjects are the included experiments of data/allen_2p/manifest.json, in inclusion order.
Each file is checked against its SHA-256 before it is loaded.

What is reused unchanged from run_things_confirmatory.py: the 12 cells, the display seeds,
the display code and the counting.

What changes (protocol Section 4):
- two time strata;
- blocks F/U/H of 10/16/24 repetitions;
- user budgets R in {8, 16};
- reference sub-blocks R_h in {8, 12, 24}.

    python3 analyses/mouse_v1/run_allen_confirmatory.py
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import hashlib, json, platform, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn
from sklearn.manifold import TSNE

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach import allen, target as tg
from analyses.eeg_confirmatory.run_things_confirmatory import CELLS, TSNE_SEEDS, MDS_SEEDS, seeds_for, embed, counts, fit_variant

RU, RH, H_SIZE = (8, 16), (8, 12, 24), 24
N_SPLITS, N_SIM = 10, 20
SPLIT_BASE, SIM_BASE = 4_500_000, 4_600_000
DATA, OUT = "data/allen_2p", "results/allen_confirmatory"


def split_rows(s, split, feats, st, truth=None, cells=CELLS, R_list=RU, n_seeds=3):
    """All rows of one split; `truth` (synthetic tests only) maps a window to true means."""
    b = allen.split_blocks(st, np.random.default_rng(SPLIT_BASE + 1000 * s + split))
    rows = []
    for rep, win in sorted({(c["rep"], c["win"]) for _, c in cells}):
        feat = feats[win]
        W, m0 = fit_variant(rep, feat[:, b["F"]])
        Y = (feat - m0) @ W
        YU, YH = Y[:, b["U"]], Y[:, b["H"]]
        ref = {Rh: [tg.rank_matrix(tg.pairwise(YH[:, j * Rh:(j + 1) * Rh].mean(axis=1)))
                    for j in range(H_SIZE // Rh)] for Rh in RH}
        trank = None if truth is None else tg.rank_matrix(tg.pairwise((truth[win] - m0) @ W))
        vcells = [(nm, c) for nm, c in cells if (c["rep"], c["win"]) == (rep, win)]
        for R in R_list:
            Xu = YU[:, :R]
            cmean = Xu.mean(axis=1)
            nrank = {"P1": tg.rank_matrix(tg.pairwise(cmean))}
            d2cv = None
            if any(c["pipe"] == "P2" for _, c in vcells):
                from reach.things import cv_sq_dist
                d2cv = cv_sq_dist(Xu, *allen.halves(R))
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
    """Truth = H means; user and reference blocks resample U residuals within image."""
    feat = feats["full"]
    W, m0 = fit_variant("white", feat[:, b["F"]])
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
                      for _ in range(H_SIZE // Rh)]
                rows.append(dict(base, Rh=Rh, **counts(N, rr, kp)))
    return rows


def included():
    man = json.load(open(os.path.join(DATA, "manifest.json")))
    return [m for m in man if m["included"]][:10]


def verified_load(entry):
    path = os.path.join(DATA, f"{entry['experiment_id']}.nwb")
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    if h.hexdigest() != entry["sha256"]:
        raise IOError(f"{entry['experiment_id']}: SHA-256 mismatch against the manifest")
    return allen.load_responses(path)


def subject(s, entry):
    t0 = time.time()
    feats, info = verified_load(entry)
    st = allen.strata()
    rows, b0 = [], None
    for split in range(N_SPLITS):
        r, b = split_rows(s, split, feats, st)
        rows += r
        b0 = b if split == 0 else b0
    sim = semi_synthetic(s, feats, b0)
    return rows, sim, dict(subject=s, experiment_id=entry["experiment_id"], **info,
                           seconds=time.time() - t0)


def main():
    subs = included()
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(protocol="docs/note_protocol_allen.md",
                   subjects=[dict(subject=i + 1, experiment_id=e["experiment_id"], donor_id=e["donor_id"])
                             for i, e in enumerate(subs)],
                   cells=CELLS, RU=RU, RH=RH, blocks=allen.BLOCKS, windows=allen.WINDOWS,
                   n_splits=N_SPLITS, n_sim=N_SIM, split_base=SPLIT_BASE, sim_base=SIM_BASE,
                   tsne_seeds=TSNE_SEEDS, mds_seeds=MDS_SEEDS, python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__,
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    t0 = time.time()
    with ProcessPoolExecutor(len(subs)) as ex, \
            open(os.path.join(OUT, "rows.jsonl"), "w") as fr, \
            open(os.path.join(OUT, "semi_synthetic.jsonl"), "w") as fs:
        futs = {ex.submit(subject, i + 1, e): i + 1 for i, e in enumerate(subs)}
        infos = []
        for fu in as_completed(futs):
            rows, sim, info = fu.result()
            for r in rows:
                fr.write(json.dumps(r) + "\n")
            for r in sim:
                fs.write(json.dumps(r) + "\n")
            infos.append(info)
            print(f"  subject {info['subject']} ({info['experiment_id']}, {info['n_cells']} cells) "
                  f"done in {info['seconds']:.0f}s", flush=True)
    json.dump(sorted(infos, key=lambda d: d["subject"]), open(os.path.join(OUT, "subjects.json"), "w"), indent=1)
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
