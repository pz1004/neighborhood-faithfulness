#!/usr/bin/env python3
"""The simulation study of docs/revision_protocol.md (plan Amendment 6): Layers A and B.

    python3 run_validation.py --phase dev --alpha        # T3's alpha from the bases; no simulation
    python3 run_validation.py --phase dev                # the development pilot, 2 repetitions per base
    python3 run_validation.py --phase validation         # 10 repetitions per base, after the addendum

One task is one dataset, cell and base subject. A task loads its subject's data, checked against
the recorded SHA-256, builds the base (reach/sim.py) and runs its repetitions. The rows go to
results/revision/<phase>/<dataset>.jsonl.gz and the times to time.jsonl.

Development and validation seeds start at 6,000,000 and 6,500,000 (protocol Section 9). A
validation run refuses to start unless the addendum and the decisions file are committed
(protocol Section 8).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, gzip, json, platform, re, subprocess, time, zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import scipy, sklearn

from reach import sim

SEED_BASE = {"dev": 6_000_000, "validation": 6_500_000}
REPS = {"dev": 2, "validation": 10}               # per base subject: 20 and 100 per cell
DS_INDEX = {"EEG": 0, "mouse": 1}
OUT = "results/revision"
DECISIONS = f"{OUT}/dev/decisions.json"
ADDENDUM = "docs/revision_addendum.md"
ALPHA_GRID = np.round(np.linspace(0.0, 1.0, 21), 2)
LAMS = (0.1, 0.5, 1.0, 2.0)


def cells(ds, lams=LAMS):
    """The cells of protocol Section 3.4, in the order they are run."""
    RU = sim.LAYOUT[ds]["RU"]
    out = []
    for fam in ("T1", "T2", "T3"):
        for lam in lams:
            full = lam == 1.0 and fam in ("T1", "T3")            # the regime of the real data, at every budget
            out.append(dict(name=f"{fam}_l{lam:g}", family=fam, map="tsne", R_list=RU if full else (8,),
                            lam_F=lam, lam_U=lam, lam_H=lam, rho=0.0))
    out.append(dict(name="T1_zero", family="T1", map="tsne", R_list=(8,), lam_F=1.0, lam_U=0.0, lam_H=0.0, rho=0.0))
    for rho in (0.1, 0.3):
        out.append(dict(name=f"T1_D{rho:g}", family="T1", map="tsne", R_list=(8,),
                        lam_F=1.0, lam_U=1.0, lam_H=1.0, rho=rho))
    for fam in ("T1", "T2"):
        out.append(dict(name=f"{fam}_S", family=fam, map="tsne", R_list=(8,), lam_F=1.0, lam_U=0.1, lam_H=1.0, rho=0.0))
    for fam in ("T1", "T2", "T3"):
        out.append(dict(name=f"{fam}_cvmds", family=fam, map="mds_cv", R_list=(8,),
                        lam_F=1.0, lam_U=1.0, lam_H=1.0, rho=0.0))
    return out


def cell_seed(phase, ds, cell, s, rep):
    """The tuple seed of one repetition; the cell enters by a checksum of its name, not its position."""
    return (SEED_BASE[phase], DS_INDEX[ds], zlib.crc32(cell["name"].encode()), s, rep)


def target_gaps():
    """The cross-validated relative gaps the T3 truth is matched to (results/allen_confirmatory/explore_snr.txt)."""
    txt = open("results/allen_confirmatory/explore_snr.txt").read()
    eeg, mouse = txt.split("THINGS-EEG2")[1], txt.split("THINGS-EEG2")[0]
    pick = lambda t: float(re.search(r"relative gap, 8th to 20th true neighbor: ([0-9.]+)", t).group(1))
    return {"EEG": pick(eeg), "mouse": pick(mouse)}


def require_committed(paths):
    """Protocol Section 8: no validation seed before the addendum and the decisions are committed."""
    for f in paths:
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", f], capture_output=True).returncode == 0
        clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", f]).returncode == 0
        if not (tracked and clean):
            raise SystemExit(f"the validation needs {f} committed first (docs/revision_protocol.md, Section 8)")


def gap_task(args):
    ds, s, r = args
    base = sim.Base(ds, s)
    return ds, s, {f"{a:.2f}": sim.relative_gap(base, sim.truth(base, "T3", r=r, alpha=a)) for a in ALPHA_GRID}


def run_alpha(workers, r=5):
    """Choose T3's alpha per dataset: the grid value whose median gap over subjects is nearest the target."""
    targets = target_gaps()
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(gap_task, [(ds, s, r) for ds in DS_INDEX for s in range(1, 11)]))
    out = {}
    for ds in DS_INDEX:
        med = {a: float(np.median([g[a] for d, s, g in res if d == ds])) for a in res[0][2]}
        best = min(med, key=lambda a: (abs(med[a] - targets[ds]), float(a)))
        out[ds] = dict(target=targets[ds], alpha=float(best), median_gap=med,
                       per_subject={s: g for d, s, g in res if d == ds})
        print(f"{ds}: target gap {targets[ds]:.3f}; T2 (alpha 0) {med['0.00']:.3f}, T1 (alpha 1) {med['1.00']:.3f};"
              f" alpha {best} gives {med[best]:.3f}")
    os.makedirs(f"{OUT}/dev", exist_ok=True)
    json.dump(dict(r=r, grid=ALPHA_GRID.tolist(), **out), open(f"{OUT}/dev/alpha.json", "w"), indent=1)


def task(args):
    phase, ds, cell, s, reps, r, alpha = args
    t0 = time.time()
    base = sim.Base(ds, s)
    mu = sim.truth(base, cell["family"], r=r, alpha=alpha)
    rows, times = [], []
    for rep in reps:
        t1 = time.time()
        for row in sim.run_repetition(base, cell, mu, cell_seed(phase, ds, cell, s, rep)):
            rows.append(dict(ds=ds, cell=cell["name"], base=s, rep=rep, **row))
        times.append(dict(ds=ds, cell=cell["name"], base=s, rep=rep, seconds=round(time.time() - t1, 2)))
    return args, rows, times, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=sorted(SEED_BASE), required=True)
    ap.add_argument("--alpha", action="store_true", help="only choose T3's alpha (development)")
    ap.add_argument("--datasets", default="EEG,mouse")
    ap.add_argument("--cells", default="", help="comma-separated subset of cell names")
    ap.add_argument("--workers", type=int, default=35)
    ap.add_argument("--lams", default="", help="development only: noise scales of the screening cells")
    ap.add_argument("--append", action="store_true", help="development only: add rows to an earlier run")
    a = ap.parse_args()
    assert a.phase == "dev" or not (a.lams or a.append), "the validation's grid comes from the decisions file"
    if a.alpha:
        assert a.phase == "dev", "alpha is chosen in the development pilot"
        return run_alpha(a.workers)
    if a.phase == "validation":
        require_committed([ADDENDUM, DECISIONS])
        dec = json.load(open(DECISIONS))
        r, alpha, lams = dec["r"], dec["alpha"], tuple(dec["lams"])
    else:
        al = json.load(open(f"{OUT}/dev/alpha.json"))
        r, alpha = al["r"], {ds: al[ds]["alpha"] for ds in DS_INDEX}
        lams = tuple(float(x) for x in a.lams.split(",")) if a.lams else LAMS
    dss = a.datasets.split(",")
    keep = set(filter(None, a.cells.split(",")))
    jobs = []
    for ds in dss:
        for cell in cells(ds, lams):
            if keep and cell["name"] not in keep:
                continue
            for s in range(1, 11):
                jobs.append((a.phase, ds, cell, s, list(range(REPS[a.phase])), r, alpha[ds]))
    # the costliest first: several budgets, then MDS maps, then EEG before mouse
    jobs.sort(key=lambda j: (-len(j[2]["R_list"]), j[2]["map"] != "mds_cv", j[1] != "EEG"))
    out = f"{OUT}/{a.phase}"
    os.makedirs(out, exist_ok=True)
    mode = "a" if a.append else "w"
    json.dump(dict(protocol="docs/revision_protocol.md", phase=a.phase, seed_base=SEED_BASE[a.phase],
                   reps_per_base=REPS[a.phase], r=r, alpha=alpha, lams=lams, datasets=dss,
                   cells={ds: cells(ds, lams) for ds in dss}, python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__,
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(f"{out}/manifest{'_' + a.lams.replace(',', '_') if a.append else ''}.json", "w"), indent=1, default=list)
    t0 = time.time()
    files = {ds: gzip.open(f"{out}/{ds}.jsonl.gz", mode + "t") for ds in dss}
    with ProcessPoolExecutor(a.workers) as ex, open(f"{out}/time.jsonl", mode) as ft:
        futs = [ex.submit(task, j) for j in jobs]
        for i, fu in enumerate(as_completed(futs), 1):
            (phase, ds, cell, s, reps, _, _), rows, times, sec = fu.result()
            for row in rows:
                files[ds].write(json.dumps(row) + "\n")
            for t in times:
                ft.write(json.dumps(t) + "\n")
            print(f"  [{i}/{len(jobs)}] {ds} {cell['name']} base {s}: {len(reps)} repetitions in {sec:.0f}s",
                  flush=True)
    for f in files.values():
        f.close()
    print(f"done in {time.time() - t0:.0f}s -> {out}")


if __name__ == "__main__":
    main()
