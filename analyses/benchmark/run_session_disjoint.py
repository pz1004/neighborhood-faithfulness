#!/usr/bin/env python3
"""EEG blocks from separate sessions (docs/revision_protocol.md, Section 6, item 2): a sensitivity.

For each subject and each of the 12 ordered pairs of sessions (a, b):
- F holds 8 random repetitions of session a;
- U holds 16 of session b, in random order; the analyst uses its first R, R in {8, 16};
- H holds the 40 repetitions of the other two sessions, alternating between them, so every
  reference sub-block is balanced over the two.
The control is 12 random splits with the same block sizes, balanced over all four sessions
(2, 4 and 10 repetitions of each). Everything else is the primary cell: the representation
fitted on F, t-SNE with the three display seeds, plug-in references at R_h in {8, 20, 40}. The
references of the revision's candidate estimators are counted too: plug-in at R_h = 10, and
cross-validated (reach/sim.py).

The rows match the confirmatory rows, with the design ('disjoint' or 'control') and its index.
analyze_session_disjoint.py compares them.

    python3 analyses/benchmark/run_session_disjoint.py > results/session_disjoint/run.log
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import itertools, json, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach import sim, things, target as tg
from analyses.eeg_confirmatory import run_things_confirmatory as RT

SEED = 8_100_000
RU, RH, H = (8, 16), RT.RH, 40
PAIRS = list(itertools.permutations(range(4), 2))           # 12 ordered (F session, U session)
OUT = "results/session_disjoint"


def disjoint_blocks(ses, a, b, rng):
    """F from session a, U from session b, H alternating over the other two sessions."""
    F = rng.permutation(np.flatnonzero(ses == a))[:8]
    U = rng.permutation(np.flatnonzero(ses == b))[:16]
    c, d = [x for x in range(4) if x not in (a, b)]
    Hc, Hd = rng.permutation(np.flatnonzero(ses == c)), rng.permutation(np.flatnonzero(ses == d))
    return dict(F=F, U=U, H=np.stack([Hc, Hd], axis=1).reshape(-1))


def control_blocks(ses, rng):
    """The published splitter's layout with 2, 4 and 10 repetitions per session."""
    per = {name: [] for name in ("F", "U", "H")}
    for s_ in range(4):
        pos = rng.permutation(np.flatnonzero(ses == s_))
        o = 0
        for name, m in (("F", 2), ("U", 4), ("H", 10)):
            per[name].append(pos[o:o + m])
            o += m
    return {name: np.stack(per[name], axis=1).reshape(-1) for name in per}


def rows_for(feat, b, base):
    W, m0 = RT.fit_variant("white", feat[:, b["F"]])
    Y = (feat - m0) @ W
    YU, YH = Y[:, b["U"]], Y[:, b["H"]]
    ref = sim.reference_ranks(YH, sim.LAYOUT["EEG"])
    m, kp = tg.m_of(0.2, 10), 20
    rows = []
    for R in RU:
        c = YU[:, :R].mean(axis=1)
        nrank = tg.rank_matrix(tg.pairwise(c))
        for seed in RT.TSNE_SEEDS:
            E_idx, _ = tg.displayed(RT.embed("tsne30", c, seed), 10)
            N = tg.reach_from_ranks(nrank, E_idx, m) <= kp
            for (fam, Rh), ranks in sorted(ref.items()):
                rr = [tg.reach_from_ranks(rk, E_idx, m) for rk in ranks]
                rows.append(dict(base, R=R, seed=seed, ref=fam, Rh=Rh, **RT.counts(N, rr, kp)))
    return rows


def subject(s):
    t0 = time.time()
    X, _, times = RT.verified_load(s)
    feat = things.features_window(X, times, *RT.WINDOWS["full"]).astype(np.float64)
    del X
    ses = things.sessions()
    rows = []
    for i, (a, b) in enumerate(PAIRS):
        blk = disjoint_blocks(ses, a, b, np.random.default_rng((SEED, s, 0, i)))
        rows += rows_for(feat, blk, dict(subject=s, design="disjoint", index=i, f_session=a + 1, u_session=b + 1))
    for i in range(len(PAIRS)):
        blk = control_blocks(ses, np.random.default_rng((SEED, s, 1, i)))
        rows += rows_for(feat, blk, dict(subject=s, design="control", index=i, f_session=0, u_session=0))
    return s, rows, time.time() - t0


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    with ProcessPoolExecutor(10) as ex, open(f"{OUT}/rows.jsonl", "w") as f:
        for fu in as_completed([ex.submit(subject, s) for s in range(1, 11)]):
            s, rows, sec = fu.result()
            for r in rows:
                f.write(json.dumps(r) + "\n")
            print(f"  sub-{s:02d}: {len(rows)} rows in {sec:.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s -> {OUT}")


if __name__ == "__main__":
    main()
