#!/usr/bin/env python3
"""Gate A: each real result's status and reporting class (docs/revision_protocol.md, Section 5, as
applied by docs/revision_addendum.md, Section 2).

Every real cell gets the observable status of Section 4, under the frozen thresholds. Real cells are
dataset x confirmatory cell x budget, for false comfort, the mislabel rate (claim C3) and each common
check. Its matching validation cells are then those of the same dataset and map type, with the same
status, whose median d(R_h,max) and fall lie within 0.05 of the real cell's. The dependent-noise and
stress cells never match. The worst matching cell decides the accuracy. Two classifications follow:
- strict: no candidate estimator qualified (Section 3.5), so every latent estimate is model-dependent
  at best;
- local: Section 5 cell by cell; a threshold claim further needs Layer C.
For every real cell, the finite-reference disagreement at the deepest reference is reported too,
with a 95% subject-bootstrap interval, since it is reportable whatever the class.

Reads the committed confirmatory and checks rows, results/revision/dev/decisions.json and
results/revision/validation/summary.json. Writes results/gate_a/classification.txt and .json.
"""
import collections, json, os
import numpy as np

import analyze_validation as AV
from reach import checks, sim

OUT = "results/gate_a"
REAL = {"EEG": dict(rows="results/things_confirmatory/rows.jsonl", checks="results/checks/things/rows.jsonl",
                    L=sim.LAYOUT["EEG"]),
        "mouse": dict(rows="results/allen_confirmatory/rows.jsonl", checks="results/checks/allen/rows.jsonl",
                      L=sim.LAYOUT["mouse"])}
MAP_OF = {"rsa_cv_mds": "mds_cv", "pca2": "pca2", "mds": "mds"}          # every other cell is a t-SNE map
CLAIMS = {("EEG", "primary", 8, "FC"): "C1", ("EEG", "rsa_cv_mds", 8, "FC"): "C1-RSA",
          ("EEG", "primary", 32, "FC"): "C2", ("mouse", "primary", 8, "FC"): "C1",
          ("mouse", "rsa_cv_mds", 8, "FC"): "C1-RSA"}
CLAIMS.update({(ds, "primary", R, "q"): "C3" for ds, Rs in (("EEG", (8, 16, 32)), ("mouse", (8, 16))) for R in Rs})
B, SEED = 10_000, 9_100_000


def pooled_cells(path):
    """{(cell, R): {subject: counts}} from the confirmatory rows, pooled over splits and seeds."""
    acc = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(
        lambda: np.zeros(len(sim.FIELDS)))))
    for line in open(path):
        r = json.loads(line)
        acc[(r["cell"], r["R"])][r["subject"]][("plug", r["Rh"])] += [
            r["n1"] * r["nb"], r["n10"], r["dis"], r["t1"], r["n"] * r["nb"], r["n1"], r["n"]]
    return acc


def pooled_checks(path):
    """{(check, R): {subject: counts}}: each check's joint counts in the confirmatory counts' layout."""
    acc = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(
        lambda: np.zeros(len(sim.FIELDS)))))
    for line in open(path):
        r = json.loads(line)
        for X in checks.NAMES[1:]:
            acc[(X, r["R"])][r["subject"]][("plug", r["Rh"])] += [
                r[f"{X}_pass"] * r["nb"], r[f"{X}_a"], 0, r["t1"], r["n"] * r["nb"], r[f"{X}_pass"], r["n"]]
    return acc


def deepest(recs, L, rng):
    """Median finite-reference disagreement at the deepest reference, with a 95% subject-bootstrap interval."""
    v = np.array([r["curve"][f"plug|{L['RH'][-1]}"] for r in recs], float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return dict(median=float(np.nanmedian(v)) if len(v) else np.nan, ci=[np.nan, np.nan], n=len(v))
    m = np.median(v[rng.integers(0, len(v), size=(B, len(v)))], axis=1)
    return dict(median=float(np.median(v)), ci=[float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))], n=len(v))


def classify(status, matches, quantity):
    """The strict and local classes of Section 5 from the status and the matching cells' accuracy."""
    pkey = "perf_q" if quantity == "q" else "perf"
    acc = [(m[pkey]["smaller"].get("median_bias", np.nan), m[pkey]["smaller"].get("rmse", np.nan)) for m in matches]
    tight = bool(acc) and all(abs(b) <= AV.BIAS_OK and r <= AV.RMSE_OK for b, r in acc)
    loose = bool(acc) and all(abs(b) <= 0.05 and r <= 0.10 for b, r in acc)
    if status != "estimable":
        return "finite-reference only", "finite-reference only", acc
    if not matches:
        return "model-dependent (no matching validation cell)", "model-dependent (no matching validation cell)", acc
    strict = "model-dependent" if loose else "finite-reference only"
    local = "estimable: threshold claims need Layer C" if tight else ("model-dependent" if loose else "finite-reference only")
    return strict, local, acc


def main():
    dec = json.load(open("results/revision/dev/decisions.json"))
    val = json.load(open("results/revision/validation/summary.json"))
    th = dec["thresholds"]
    rng = np.random.default_rng(SEED)
    os.makedirs(OUT, exist_ok=True)
    lines, out_json = [], {}

    def out(s=""):
        lines.append(s)
        print(s)

    out("GATE A (analyze_gate_a.py): status, matching validation cells and class of every real result")
    out(f"thresholds {th}; accuracy of the smaller extrapolant; worst matching cell decides")
    for ds, P in REAL.items():
        L = P["L"]
        vcells = {k: v for k, v in val[ds]["cells"].items() if not v["diagnostic"]}
        units = [(cell, R, "FC", subj) for (cell, R), subj in pooled_cells(P["rows"]).items()]
        units += [(cell, R, "q", subj) for cell, R, q_, subj in units if cell == "primary"]
        units += [(f"check:{X}", R, "FC", subj) for (X, R), subj in pooled_checks(P["checks"]).items()]
        out(f"\n=== {ds} ===")
        out(f"  {'result':22s} {'R':>3s} {'q':>2s}  {'range (smaller-larger)':>22s}  {'d_max':>6s} {'fall':>6s} {'fd':>7s} {'div':>6s}"
            f"  {'status':25s} {'matches':>7s}  {'worst bias/RMSE':>16s}  {'d at deepest [95%]':>22s}  strict / local")
        for cell, R, qn, subj in sorted(units, key=lambda u: (u[0].startswith("check"), u[0], u[1], u[2])):
            recs = [AV.record(dict(c), L) for s, c in sorted(subj.items())]
            if qn == "q":
                for r in recs:                                    # judge the claim on the q curve's own estimates
                    r["smaller"], r["larger"] = r["q_smaller"], max(r["q_sqrt"], r["q_inv"])
            st = AV.status(recs, th, "plug")
            obs = {k: float(np.nanmedian([r["obs_plug"][k] for r in recs])) for k in ("dmax", "fall", "fd", "div")}
            mp = "tsne" if cell.startswith("check:") else MAP_OF.get(cell, "tsne")
            matches = [v for k, v in vcells.items()
                       if ("cvmds" in k) == (mp == "mds_cv") and mp in ("tsne", "mds_cv") and v["status"] == st
                       and abs(v["obs"]["dmax"] - obs["dmax"]) <= 0.05 and abs(v["obs"]["fall"] - obs["fall"]) <= 0.05]
            names = [k for k, v in vcells.items() if any(v is m for m in matches)]
            strict, local, acc = classify(st, matches, qn)
            worst = max(acc, key=lambda a: abs(a[0])) if acc else (np.nan, np.nan)
            lo, hi = (float(np.nanmedian([r["smaller"] for r in recs])), float(np.nanmedian([r["larger"] for r in recs])))
            fin = deepest(recs, L, rng)
            claim = CLAIMS.get((ds, cell, R, qn), "")
            key = f"{cell}|{R}|{qn}"
            out_json.setdefault(ds, {})[key] = dict(claim=claim, map=mp, status=st, obs=obs, range=[lo, hi],
                                                    matches=names, worst=list(worst), strict=strict, local=local,
                                                    finite_reference=fin,
                                                    note="neighborhood definition differs from the validation's"
                                                    if cell in ("k5", "k20", "eps0.1", "eps0.3") else "")
            out(f"  {(claim + ' ' if claim else '') + cell:22s} {R:3d} {qn:>2s}  {lo:10.3f}-{hi:<10.3f}  {obs['dmax']:6.3f}"
                f" {obs['fall']:6.3f} {obs['fd']:+7.3f} {obs['div']:6.3f}  {st:25s} {len(matches):7d}"
                f"  {worst[0]:+.3f}/{worst[1]:.3f}    {fin['median']:.3f} [{fin['ci'][0]:.3f}, {fin['ci'][1]:.3f}]"
                f"  {strict} / {local}")
    open(f"{OUT}/classification.txt", "w").write("\n".join(lines) + "\n")
    json.dump(out_json, open(f"{OUT}/classification.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
