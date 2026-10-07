#!/usr/bin/env python3
"""Per-point rerun of the primary and cv-MDS maps and of every check (docs/revision_protocol.md,
Section 6, items 4, 6 and 7).

For each subject, split and budget, the confirmatory maps are re-created with the committed seeds.
Every check is recomputed as reach/checks.py computes it, keeping the score behind each label,
and each published generator is consumed exactly as before. The run asserts that every committed
count is reproduced exactly: the rows of the comparison of checks (every check, every reference
size), and the confirmatory rows of the cv-MDS cell.

Per map it stores:
- the displayed neighbors and the naive reach;
- per point, every check's score: per-point trustworthiness; the scDEED score, with the cut-off of
  the registered null and of the alternative null; the cross-validated reach; the number of other
  seeds in whose maps the neighborhood recurs; the bootstrap support; recurrence in each half;
- the reach of the displayed neighborhood in every reference sub-block, plug-in at the published
  reference sizes and at candidate (b)'s extra size, and cross-validated (candidate (c)), with the
  continuous overlap |E_i and N_k'(i)| in each;
- for cv-MDS maps, the share of negative cross-validated squared distances, the share of zeros in
  the MDS input and the Kruskal stress-1 of the fit.

The alternative scDEED null permutes the representation's 20 coordinates independently across
images, as the null was first written (Supplementary Section S1.4). It uses its own seeds,
(8,200,000, ...).

    python3 run_points.py --dataset things|allen [--workers 35]
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

from reach import checks, sim, things, target as tg
import run_checks as RC
import run_things_confirmatory as RT

ALT_SEED = 8_200_000
K, EPS = RC.K, RC.EPS
M, KP = tg.m_of(EPS, K), 2 * K
LAYOUT = {"things": sim.LAYOUT["EEG"], "allen": sim.LAYOUT["mouse"]}
CONF = {"things": "results/things_confirmatory/rows.jsonl", "allen": "results/allen_confirmatory/rows.jsonl"}
OUT = "results/points"


def reference_index(L):
    """(family, R_h, sub-block) of every reference sub-block, in the order sim.reference_ranks makes them."""
    return [(fam, Rh, j) for Rh in L["RH"] + (L["extra"],) for fam in ("plug", "cv") for j in range(L["H"] // Rh)]


def ref_arrays(refs, L, E_idx):
    """Reach and continuous overlap of the displayed neighborhoods in every reference sub-block."""
    reach, overlap = [], []
    for fam, Rh, j in reference_index(L):
        rk = refs[(fam, Rh)][j]
        reach.append(tg.reach_from_ranks(rk, E_idx, M))
        overlap.append((np.take_along_axis(rk, E_idx, axis=1) <= KP).sum(axis=1))
    return np.array(reach, np.uint8), np.array(overlap, np.uint8)


def stress1(D, emb):
    """Kruskal stress-1 of an MDS fit of the dissimilarities D (zero diagonal)."""
    iu = np.triu_indices(len(D), 1)
    d = tg.pairwise(emb)[iu]
    return float(np.sqrt(((d - D[iu]) ** 2).sum() / (d ** 2).sum()))


def scores_for_budget(Xu, c_raw, project, halves, seeds, rng_boot, rng_perm, rng_alt):
    """reach/checks.labels_for_budget with every score kept; labels follow from them by the same rules."""
    A, B = halves
    c = Xu.mean(axis=1)
    rank = tg.rank_matrix(tg.pairwise(c))
    cv_rank = tg.rank_matrix(things.cv_sq_dist(Xu, A, B))
    cA, cB = Xu[:, A].mean(axis=1), Xu[:, B].mean(axis=1)
    embs = [checks.tsne(c, s) for s in seeds]
    E = [tg.displayed(e, K)[0] for e in embs]
    out = []
    for j, s in enumerate(seeds):
        naive = tg.reach_from_ranks(rank, E[j], M)
        N = naive <= KP
        tw = checks.trustworthiness_points(rank, E[j], K)
        perm = project(np.column_stack([col[rng_perm[j].permutation(len(col))] for col in c_raw.T]))
        null = checks.scdeed_scores(perm, checks.tsne(perm, s))
        scd = checks.scdeed_scores(c, embs[j])
        cut = np.percentile(null, checks.SCDEED_PCT)
        alt = np.column_stack([col[rng_alt[j].permutation(len(col))] for col in c.T])
        cut_alt = np.percentile(checks.scdeed_scores(alt, checks.tsne(alt, s)), checks.SCDEED_PCT)
        cv = tg.reach_from_ranks(cv_rank, E[j], M)
        rec = [checks.recurs(E[j], embs[i], M, KP) for i in range(len(seeds)) if i != j]
        support = checks.bootstrap_support(Xu, E[j], s, M, KP, rng_boot[j])
        shA = checks.recurs(E[j], checks.tsne(cA, s), M, KP)
        shB = checks.recurs(E[j], checks.tsne(cB, s), M, KP)
        lab = dict(N=N, TW=checks.top_matched(tw, int(N.sum())), SCD=scd >= cut, CV=cv <= KP,
                   SEED=np.logical_and.reduce(rec), BOOT=support >= checks.BOOT_SUPPORT - 1e-9, SH=shA & shB)
        out.append((E[j], lab, dict(naive=naive, tw=tw, scd=scd, scd_cut=cut, scd_cut_alt=cut_alt, cv=cv,
                                    seed_count=np.sum(rec, axis=0), boot=support, shA=shA, shB=shB)))
    return out


def task(args):
    ds, s, split, committed, committed_mds = args
    t0 = time.time()
    P, L = RC.DATASETS[ds], LAYOUT[ds]
    feat = RC.load(ds, s)
    _, _, halves = RC.layout(ds)
    b, YU, YH, project = RC.blocks_and_space(ds, s, split, feat)
    refs = sim.reference_ranks(YH, L)
    pub = {Rh: refs[("plug", Rh)] for Rh in L["RH"]}
    maps, mds, bad = [], [], []
    for R in P["RU"]:
        Xu = YU[:, :R]
        c_raw = feat[:, b["U"][:R]].mean(axis=1)
        rb = [np.random.default_rng((P["boot"], s, split, R, j)) for j in range(3)]
        rp = [np.random.default_rng((P["perm"], s, split, R, j)) for j in range(3)]
        ra = [np.random.default_rng((ALT_SEED, int(ds == "allen"), s, split, R, j)) for j in range(3)]
        for (E_idx, lab, sc), seed in zip(scores_for_budget(Xu, c_raw, project, halves(R), RT.TSNE_SEEDS, rb, rp, ra),
                                          RT.TSNE_SEEDS):
            for Rh in L["RH"]:
                row = checks.check_counts(lab, [tg.reach_from_ranks(rk, E_idx, M) for rk in pub[Rh]], KP)
                if committed.get((R, seed, Rh)) != row:
                    bad.append(("checks", s, split, R, seed, Rh))
            reach, overlap = ref_arrays(refs, L, E_idx)
            maps.append(dict(meta=(split, R, seed), E=E_idx.astype(np.uint8), ref_reach=reach, ref_overlap=overlap,
                             **{k: np.asarray(v) for k, v in sc.items()}))
        # the cv-MDS cell's maps, as the confirmatory run made them
        d2cv = things.cv_sq_dist(Xu, *halves(R))
        cv_rank = tg.rank_matrix(d2cv)
        D = np.sqrt(np.maximum(d2cv, 0.0))
        np.fill_diagonal(D, 0.0)
        iu = np.triu_indices(len(D), 1)
        for seed in RT.MDS_SEEDS:
            emb = RT.embed("mds_cv", Xu.mean(axis=1), seed, d2cv)
            E_idx, _ = tg.displayed(emb, K)
            naive = tg.reach_from_ranks(cv_rank, E_idx, M)
            N = naive <= KP
            for Rh in L["RH"]:
                row = RT.counts(N, [tg.reach_from_ranks(rk, E_idx, M) for rk in pub[Rh]], KP)
                if committed_mds.get((R, seed, Rh)) != row:
                    bad.append(("cv-MDS", s, split, R, seed, Rh))
            reach, overlap = ref_arrays(refs, L, E_idx)
            mds.append(dict(meta=(split, R, seed), E=E_idx.astype(np.uint8), naive=naive, ref_reach=reach,
                            ref_overlap=overlap, stress=stress1(D, emb), neg_share=float((d2cv[iu] < 0).mean()),
                            zero_share=float((D[iu] == 0).mean())))
    return (ds, s, split), maps, mds, bad, time.time() - t0


def stack(items):
    keys = items[0].keys()
    return {k: np.array([it[k] for it in items]) for k in keys}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(RC.DATASETS))
    ap.add_argument("--workers", type=int, default=35)
    a = ap.parse_args()
    ds = a.dataset
    keys = ("n", "nb", "n1", "n10", "dis", "t1") + tuple(f"{X}_{q}" for X in checks.NAMES for q in ("pass", "nx", "a", "b"))
    committed = {}
    for line in open(f"{RC.DATASETS[ds]['out']}/rows.jsonl"):
        r = json.loads(line)
        committed.setdefault((r["subject"], r["split"]), {})[(r["R"], r["seed"], r["Rh"])] = {k: r[k] for k in keys}
    committed_mds = {}
    for line in open(CONF[ds]):
        r = json.loads(line)
        if r["cell"] == "rsa_cv_mds":
            committed_mds.setdefault((r["subject"], r["split"]), {})[(r["R"], r["seed"], r["Rh"])] = \
                {k: r[k] for k in ("n", "n1", "n10", "dis", "t1", "nb")}
    out = f"{OUT}/{ds}"
    os.makedirs(out, exist_ok=True)
    jobs = [(ds, s, split, committed[(s, split)], committed_mds[(s, split)]) for s in range(1, 11) for split in range(10)]
    per = {s: ([], []) for s in range(1, 11)}
    bad_all, t0 = [], time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(task, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            (_, s, split), maps, mds, bad, sec = fu.result()
            per[s][0].extend(maps)
            per[s][1].extend(mds)
            bad_all += bad
            print(f"  [{i}/{len(jobs)}] subject {s} split {split}: {len(maps)} maps, {len(mds)} MDS maps,"
                  f" {len(bad)} counts not reproduced, {sec:.0f}s", flush=True)
    idx = np.array([(0 if f == "plug" else 1, Rh, j) for f, Rh, j in reference_index(LAYOUT[ds])], np.int16)
    for s, (maps, mds) in per.items():
        maps.sort(key=lambda d: d["meta"])
        mds.sort(key=lambda d: d["meta"])
        np.savez_compressed(f"{out}/sub-{s:02d}.npz", ref_index=idx,
                            **{f"tsne_{k}": v for k, v in stack(maps).items()},
                            **{f"mds_{k}": v for k, v in stack(mds).items()})
    json.dump(dict(dataset=ds, not_reproduced=bad_all, seconds=time.time() - t0), open(f"{out}/reproduction.json", "w"),
              indent=1)
    print(f"{len(bad_all)} counts not reproduced; done in {time.time() - t0:.0f}s -> {out}")
    if bad_all:
        raise SystemExit("the per-point rerun did not reproduce the committed counts")


if __name__ == "__main__":
    main()
