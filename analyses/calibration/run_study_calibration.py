#!/usr/bin/env python3
"""Layer C: study-level calibration of the published bounds (docs/revision_protocol.md, Section 3.4;
docs/revision_addendum.md, Section 2).

A population of 30 virtual subjects: the 10 base subjects x lam in {0.8, 1.0, 1.25}, under one truth
family (the scenario).

--oracle computes each virtual subject's true protocol-level false comfort theta_v, and its
mislabel rate, from 200 independent datasets, each with one random split, the three display seeds
and exact labels.

Otherwise it runs the virtual studies. Each study:
- draws 10 virtual subjects with replacement;
- gives each one a new dataset;
- runs the published procedure on it (10 splits x 3 display seeds, plug-in references at the real
  grid);
- in the first 200 studies, also runs the published plasmode check (20 draws on split 0), for the
  guard.
The pooled counts are written. analyze_study_calibration.py forms the per-subject estimates, the
medians, the 10,000-resample subject bootstrap, L and U, the order-statistic bounds and the guard.

Seeds: studies (7,000,000, ...), oracle (7,500,000, ...).

    python3 analyses/calibration/run_study_calibration.py --dataset EEG --scenario T1 --map tsne --oracle
    python3 analyses/calibration/run_study_calibration.py --dataset EEG --scenario T1 --map tsne --studies 1000
    python3 analyses/calibration/run_study_calibration.py --dataset mouse --scenario T1 --R 16 --no-guard --studies 1000
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, functools, gzip, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from sklearn.manifold import TSNE

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach import sim, target as tg
from analyses.eeg_confirmatory import run_things_confirmatory as RT
from analyses.validation import run_validation as RV

LAMS = (0.8, 1.0, 1.25)
POPULATION = [(s, lam) for s in range(1, 11) for lam in LAMS]           # virtual subject v = index here
SUBJECTS, N_ORACLE, N_GUARD, N_SIM = 10, 200, 200, 20
STUDY_SEED, ORACLE_SEED = 7_000_000, 7_500_000
SCEN = {"T1": 1, "T2": 2, "T3": 3}
MAPS = {"tsne": 0, "mds_cv": 1}
OUT = "results/revision/layer_c"
R = 8


@functools.lru_cache(maxsize=1)
def base_of(ds, s):
    return sim.Base(ds, s)


def setup(ds, scenario, s, lam, mp, R=R):
    dec = json.load(open(RV.DECISIONS))
    base = base_of(ds, s)
    mu = sim.truth(base, scenario, r=dec["r"], alpha=dec["alpha"][ds])
    cell = dict(name="layer_c", family=scenario, map=mp, R_list=(R,), lam_F=lam, lam_U=lam, lam_H=lam, rho=0.0)
    return base, mu, cell


def study_unit(ds, scenario, mp, study, j, v, guard, R=R):
    """One virtual subject of one study: the published procedure's pooled counts, and the guard's plasmode."""
    s, lam = POPULATION[v]
    base, mu, cell = setup(ds, scenario, s, lam, mp, R)
    seed = (STUDY_SEED, RV.DS_INDEX[ds], SCEN[scenario], MAPS[mp], R, study, j)
    rep = sim.Repetition(base, cell, mu, np.random.default_rng(seed + (0,)))
    L, acc, first = base.L, {}, None
    for k in range(10):
        b, X, W, m0 = rep.split(np.random.default_rng(seed + (1, k)))
        Y = (X - m0) @ W
        YU, YH = Y[:, b["U"]], Y[:, b["H"]]
        refs = {("plug", Rh): [tg.rank_matrix(tg.pairwise(YH[:, i * Rh:(i + 1) * Rh].mean(axis=1)))
                               for i in range(L["H"] // Rh)] for Rh in L["RH"]}
        refs[("exact", 0)] = [tg.rank_matrix(tg.pairwise((mu - m0) @ W))]
        for E_idx, N in rep.maps(YU, R):
            for key, ranks in refs.items():
                sim.add_counts(acc, (R,) + key, N, ranks, E_idx)
        if k == 0:
            first = (YU, YH)
    counts = {f"{key[1]}|{key[2]}": a.tolist() for key, a in acc.items()}
    g = plasmode(base, first, np.random.default_rng(seed + (2,))) if guard else None
    return dict(study=study, slot=j, v=v, counts=counts, plasmode=g)


def plasmode(base, first, rng):
    """The published plasmode check on split 0 (run_things_confirmatory.semi_synthetic, at R = 8): truth =
    block H's means in the representation, trials resampled from block U's residuals; its pooled counts."""
    YU, YH = first
    L = base.L
    mu = YH.mean(axis=1)
    trank = tg.rank_matrix(tg.pairwise(mu))
    nU = YU.shape[1]
    pool = (YU - YU.mean(axis=1, keepdims=True)) * np.sqrt(nU / (nU - 1))
    n = mu.shape[0]
    acc = {}

    def draw(r):
        idx = rng.integers(0, nU, size=(n, r))
        return mu[:, None, :] + np.take_along_axis(pool, idx[:, :, None], axis=1)

    for d in range(N_SIM):
        c = draw(R).mean(axis=1)
        emb = TSNE(n_components=2, perplexity=30, init="pca", random_state=RT.TSNE_SEEDS[d % 3],
                   n_jobs=1).fit_transform(c)
        E_idx, _ = tg.displayed(emb, sim.K)
        N = tg.reach_from_ranks(tg.rank_matrix(tg.pairwise(c)), E_idx, sim.M) <= sim.KP
        sim.add_counts(acc, ("exact", 0), N, [trank], E_idx)
        for Rh in L["RH"]:
            ranks = [tg.rank_matrix(tg.pairwise(draw(Rh).mean(axis=1))) for _ in range(L["H"] // Rh)]
            sim.add_counts(acc, ("plug", Rh), N, ranks, E_idx)
    return {f"{k[0]}|{k[1]}": a.tolist() for k, a in acc.items()}


def oracle_unit(ds, scenario, mp, v, d, R=R):
    """One independent dataset of virtual subject v: one random split, three maps, exact counts only."""
    s, lam = POPULATION[v]
    base, mu, cell = setup(ds, scenario, s, lam, mp, R)
    seed = (ORACLE_SEED, RV.DS_INDEX[ds], SCEN[scenario], MAPS[mp], R, v, d)
    rep = sim.Repetition(base, cell, mu, np.random.default_rng(seed + (0,)))
    b, X, W, m0 = rep.split(np.random.default_rng(seed + (1,)))
    Y = (X - m0) @ W
    trank = tg.rank_matrix(tg.pairwise((mu - m0) @ W))
    n10 = n1 = dis = nn = 0
    for E_idx, N in rep.maps(Y[:, b["U"]], R):
        T = tg.reach_from_ranks(trank, E_idx, sim.M) <= sim.KP
        n10 += int((N & ~T).sum())
        n1 += int(N.sum())
        dis += int((N != T).sum())
        nn += len(N)
    return dict(v=v, draw=d, n10=n10, n1=n1, dis=dis, n=nn)


def task(args):
    kind, ds, scenario, mp, units, R_, guard = args
    t0 = time.time()
    if kind == "oracle":
        rows = [oracle_unit(ds, scenario, mp, v, d, R_) for v, d in units]
    else:
        rows = [study_unit(ds, scenario, mp, st, j, v, guard and st < N_GUARD, R_) for st, j, v in units]
    return rows, time.time() - t0, len(units)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="EEG", choices=sorted(RV.DS_INDEX))
    ap.add_argument("--scenario", required=True, choices=sorted(SCEN))
    ap.add_argument("--map", default="tsne", choices=sorted(MAPS))
    ap.add_argument("--oracle", action="store_true")
    ap.add_argument("--studies", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=35)
    ap.add_argument("--R", type=int, default=R, help="the analyst's budget")
    ap.add_argument("--no-guard", action="store_true", help="skip the plasmode: the guard applies to false comfort only")
    a = ap.parse_args()
    RV.require_committed([RV.ADDENDUM, RV.DECISIONS])
    tag = f"{a.dataset}_{a.scenario}_{a.map}_R{a.R}"
    os.makedirs(OUT, exist_ok=True)
    if a.oracle:
        units = [(v, d) for v in range(len(POPULATION)) for d in range(N_ORACLE)]
        chunks = [units[i:i + 20] for i in range(0, len(units), 20)]       # one virtual subject per chunk
        jobs = [("oracle", a.dataset, a.scenario, a.map, c, a.R, False) for c in chunks]
        path = f"{OUT}/oracle_{tag}.jsonl.gz"
    else:
        units = []
        for st in range(a.studies):
            pick = np.random.default_rng((STUDY_SEED, RV.DS_INDEX[a.dataset], SCEN[a.scenario], MAPS[a.map], st))
            units += [(st, j, int(v)) for j, v in enumerate(pick.integers(0, len(POPULATION), size=SUBJECTS))]
        units.sort(key=lambda u: (POPULATION[u[2]][0], u[0] >= N_GUARD, u[0], u[1]))   # by base subject: one load per run
        chunks = [units[i:i + 10] for i in range(0, len(units), 10)]
        jobs = [("study", a.dataset, a.scenario, a.map, c, a.R, not a.no_guard) for c in chunks]
        path = f"{OUT}/studies_{tag}.jsonl.gz"
    json.dump(dict(protocol="docs/revision_protocol.md", addendum="docs/revision_addendum.md", dataset=a.dataset,
                   scenario=a.scenario, map=a.map, R=a.R, guard=not a.no_guard, population=POPULATION,
                   oracle=a.oracle, studies=a.studies,
                   n_oracle=N_ORACLE, n_guard=N_GUARD, seeds=dict(study=STUDY_SEED, oracle=ORACLE_SEED),
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(path.replace(".jsonl.gz", "_manifest.json"), "w"), indent=1)
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(a.workers) as ex, gzip.open(path, "wt") as f:
        futs = [ex.submit(task, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            rows, sec, n = fu.result()
            for r in rows:
                f.write(json.dumps(r) + "\n")
            done += n
            if i % 25 == 0 or i == len(jobs):
                print(f"  [{i}/{len(jobs)}] {done} units, last chunk {sec:.0f}s, elapsed {time.time() - t0:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {path}")


if __name__ == "__main__":
    main()
