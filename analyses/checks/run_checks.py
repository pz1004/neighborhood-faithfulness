#!/usr/bin/env python3
"""Common checks against the truth (docs/checks_protocol.md, plan Amendment 5).

For each subject and each of the confirmatory splits, the primary cell's maps are re-created
with the confirmatory seeds and code. Every check of reach/checks.py is computed for every
map, from the analyst's trials only, and counted against the reference sub-blocks of block
H. The semi-synthetic estimator check adds rows whose Rh is "exact", counted against a known
truth. Every data file is checked against its SHA-256 before it is loaded.
analyze_checks.py applies the protocol's Sections 4 to 6.

    python3 analyses/checks/run_checks.py --dataset things|allen [--workers 36]
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, platform, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach import allen, things, checks, target as tg
from analyses.eeg_confirmatory import run_things_confirmatory as RT
from analyses.mouse_v1 import run_allen_confirmatory as RA

N_SPLITS, N_SIM, SIM_R = 10, 10, 8
K, EPS = 10, 0.2
DATASETS = {
    "things": dict(RU=RT.RU, RH=RT.RH, H_SIZE=40, split_base=RT.SPLIT_BASE,
                   boot=5_100_000, perm=5_110_000, sim=5_200_000, sim_boot=5_210_000,
                   sim_perm=5_220_000, out="results/checks/things"),
    "allen": dict(RU=RA.RU, RH=RA.RH, H_SIZE=RA.H_SIZE, split_base=RA.SPLIT_BASE,
                  boot=5_500_000, perm=5_510_000, sim=5_600_000, sim_boot=5_610_000,
                  sim_perm=5_620_000, out="results/checks/allen"),
}


def layout(ds):
    """(strata of the repetitions, block splitter, halves) of a dataset."""
    if ds == "things":
        return things.sessions(), things.split_blocks, things.halves
    return allen.strata(), allen.split_blocks, allen.halves


def load(ds, s):
    """Subject s's responses in the primary window: objects x repetitions x features."""
    if ds == "things":
        X, _, times = RT.verified_load(s)
        return things.features_window(X, times, *RT.WINDOWS["full"]).astype(np.float64)
    feats, _ = RA.verified_load(RA.included()[s - 1])
    return feats["full"]


def blocks_and_space(ds, s, split, feat):
    """Blocks of one split and the primary cell's representation, as in the confirmatory code.

    Returns the blocks, the representation of blocks U and H, and the fixed map `project`
    from measured features to the representation.
    """
    P = DATASETS[ds]
    st, split_blocks, _ = layout(ds)
    b = split_blocks(st, np.random.default_rng(P["split_base"] + 1000 * s + split))
    W, m0 = RT.fit_variant("white", feat[:, b["F"]])
    project = lambda A: (A - m0) @ W
    Y = project(feat)
    return b, Y[:, b["U"]], Y[:, b["H"]], project


def display_rows(lab_list, refs, RH, base, seeds):
    m, kp = tg.m_of(EPS, K), 2 * K
    rows = []
    for (E, lab), seed in zip(lab_list, seeds):
        for Rh in RH:
            rr = [tg.reach_from_ranks(rk, E, m) for rk in refs[Rh]]
            rows.append(dict(base, seed=seed, Rh=Rh, **checks.check_counts(lab, rr, kp)))
    return rows


def split_task(ds, s, split, feat=None, truth=None, R_list=None, seeds=RT.TSNE_SEEDS):
    """All rows of one split. `truth` (synthetic tests only) gives the true means."""
    P = DATASETS[ds]
    feat = load(ds, s) if feat is None else feat
    _, _, halves = layout(ds)
    b, YU, YH, project = blocks_and_space(ds, s, split, feat)
    H = P["H_SIZE"]
    refs = {Rh: [tg.rank_matrix(tg.pairwise(YH[:, j * Rh:(j + 1) * Rh].mean(axis=1)))
                 for j in range(H // Rh)] for Rh in P["RH"]}
    rows = []
    for R in (R_list or P["RU"]):
        Xu = YU[:, :R]
        c_raw = feat[:, b["U"][:R]].mean(axis=1)
        rb = [np.random.default_rng((P["boot"], s, split, R, j)) for j in range(len(seeds))]
        rp = [np.random.default_rng((P["perm"], s, split, R, j)) for j in range(len(seeds))]
        labs = checks.labels_for_budget(Xu, c_raw, project, halves(R), seeds, rb, rp,
                                        things.cv_sq_dist, K, EPS)
        base = dict(subject=s, split=split, R=R)
        rows += display_rows(labs, refs, P["RH"], base, seeds)
        if truth is not None:
            trank = tg.rank_matrix(tg.pairwise(project(truth)))
            rows += display_rows(labs, {"exact": [trank]}, ["exact"], base, seeds)
    return rows


def sim_task(ds, s, feat=None, n_sim=N_SIM, seeds=RT.TSNE_SEEDS):
    """Semi-synthetic check (protocol Section 5): truth = H means; trials resample U residuals.

    Residuals are resampled by index, in the measured features and in the representation
    alike, so the analyst's raw averages (for scDEED's null) match their trials.
    """
    P = DATASETS[ds]
    feat = load(ds, s) if feat is None else feat
    _, _, halves = layout(ds)
    b, YU, YH, project = blocks_and_space(ds, s, 0, feat)
    mu, mu_raw = YH.mean(axis=1), feat[:, b["H"]].mean(axis=1)
    trank = tg.rank_matrix(tg.pairwise(mu))
    nU = YU.shape[1]
    scale = np.sqrt(nU / (nU - 1))
    pool = (YU - YU.mean(axis=1, keepdims=True)) * scale
    pool_raw = (feat[:, b["U"]] - feat[:, b["U"]].mean(axis=1, keepdims=True)) * scale
    n = mu.shape[0]
    rows = []
    for d in range(n_sim):
        rng = np.random.default_rng((P["sim"], s, d))

        def draw(R, raw=False):
            idx = rng.integers(0, nU, size=(n, R))
            Xd = mu[:, None, :] + np.take_along_axis(pool, idx[:, :, None], axis=1)
            if not raw:
                return Xd
            return Xd, mu_raw + np.take_along_axis(pool_raw, idx[:, :, None], axis=1).mean(axis=1)

        Xu, c_raw = draw(SIM_R, raw=True)
        assert np.allclose(project(c_raw), Xu.mean(axis=1))
        rb = [np.random.default_rng((P["sim_boot"], s, d, j)) for j in range(len(seeds))]
        rp = [np.random.default_rng((P["sim_perm"], s, d, j)) for j in range(len(seeds))]
        labs = checks.labels_for_budget(Xu, c_raw, project, halves(SIM_R), seeds, rb, rp,
                                        things.cv_sq_dist, K, EPS)
        refs = {Rh: [tg.rank_matrix(tg.pairwise(draw(Rh).mean(axis=1))) for _ in range(P["H_SIZE"] // Rh)]
                for Rh in P["RH"]}
        base = dict(subject=s, draw=d, R=SIM_R)
        rows += display_rows(labs, {"exact": [trank]}, ["exact"], base, seeds)
        rows += display_rows(labs, refs, P["RH"], base, seeds)
    return rows


def task(args):
    ds, s, what = args
    t0 = time.time()
    rows = sim_task(ds, s) if what == "sim" else split_task(ds, s, what)
    return args, rows, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    ap.add_argument("--workers", type=int, default=36)
    a = ap.parse_args()
    ds, P = a.dataset, DATASETS[a.dataset]
    subs = list(range(1, 11))
    os.makedirs(P["out"], exist_ok=True)
    json.dump(dict(protocol="docs/checks_protocol.md", dataset=ds, subjects=subs,
                   checks=checks.NAMES, B_boot=checks.B_BOOT, boot_support=checks.BOOT_SUPPORT,
                   scdeed_x=checks.SCDEED_X, scdeed_pct=checks.SCDEED_PCT,
                   n_splits=N_SPLITS, n_sim=N_SIM, sim_R=SIM_R, tsne_seeds=RT.TSNE_SEEDS,
                   **{k: v for k, v in P.items() if k != "out"}, python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__,
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(P["out"], "manifest.json"), "w"), indent=1)
    jobs = [(ds, s, "sim") for s in subs] + [(ds, s, split) for s in subs for split in range(N_SPLITS)]
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex, \
            open(os.path.join(P["out"], "rows.jsonl"), "w") as fr, \
            open(os.path.join(P["out"], "semi_synthetic.jsonl"), "w") as fs:
        futs = [ex.submit(task, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            (_, s, what), rows, sec = fu.result()
            out = fs if what == "sim" else fr
            for r in rows:
                out.write(json.dumps(r) + "\n")
            print(f"  [{i}/{len(jobs)}] subject {s} {'semi-synthetic' if what == 'sim' else f'split {what}'}"
                  f" done in {sec:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {P['out']}")


if __name__ == "__main__":
    main()
