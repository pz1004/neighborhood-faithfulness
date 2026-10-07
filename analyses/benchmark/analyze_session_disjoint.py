#!/usr/bin/env python3
"""EEG blocks from separate sessions against the control (docs/revision_protocol.md, Section 6, item 2).

For each subject and design, counts are pooled over the 12 block draws and the three display
seeds, as the analyses pool over splits and seeds. For each budget it reports:
- the disagreement with each reference size, and the finite-reference d_40;
- both extrapolants and the estimator chosen by the development pilot;
- the status (frozen thresholds), and the false comfort of a random selection matched per map;
- the paired difference, separate sessions minus control, per subject.
This is a sensitivity analysis, not a measurement of leakage: across sessions the responses can
also change. Reads results/revision/dev/decisions.json for the estimator and the thresholds.

Writes results/session_disjoint/analysis.txt and summary.json.
"""
import collections, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from analyses.validation import analyze_validation as AV
from reach import sim

OUT = "results/session_disjoint"
L = sim.LAYOUT["EEG"]


def pooled(rows):
    """{(subject, design, R): {('plug', Rh): counts, ('rand', Rh): counts}} in sim.FIELDS order."""
    acc = collections.defaultdict(lambda: collections.defaultdict(lambda: np.zeros(len(sim.FIELDS))))
    for r in rows:
        a = acc[(r["subject"], r["design"], r["R"])]
        a[(r["ref"], r["Rh"])] += [r["n1"] * r["nb"], r["n10"], r["dis"], r["t1"], r["n"] * r["nb"], r["n1"], r["n"]]
        if r["ref"] == "plug":
            # a random selection of as many neighborhoods: expected count with T_ref = 0 among n1 random points
            a[("rand", r["Rh"])] += [r["n1"] * r["nb"], r["n1"] * (r["nb"] - r["t1"] / r["n"]), 0, 0, 0, 0, 0]
    return acc


def main():
    dec = json.load(open("results/revision/dev/decisions.json"))
    est, th = dec["estimator"], dec["thresholds"]
    est = "smaller" if est == "published" else est
    fam = AV.FAMILY_OF.get(est, "plug")
    rows = [json.loads(l) for l in open(f"{OUT}/rows.jsonl")]
    acc = pooled(rows)
    lines, summary = [], {}

    def out(s=""):
        lines.append(s)
        print(s)

    out(f"EEG BLOCKS FROM SEPARATE SESSIONS (analyze_session_disjoint.py): estimator '{est}', thresholds {th}")
    out("medians over 10 subjects; d(R_h) is the reference disagreement; 'random' is a matched random selection")
    for R in (8, 16):
        per = {}
        for design in ("control", "disjoint"):
            recs = []
            for s in range(1, 11):
                c = acc[(s, design, R)]
                rec = AV.record(dict(c), L)
                rand = [AV.fc(c, "rand", h) for h in L["RH"]]
                rec["random"] = AV.intercept(L["RH"], rand, "sqrt")
                rec["d40"] = AV.fc(c, "plug", 40)
                recs.append(rec)
            per[design] = recs
            st = AV.status(recs, th, fam)
            med = lambda k: float(np.nanmedian([r[k] for r in recs]))
            curve = [float(np.nanmedian([r["curve"][f"plug|{h}"] for r in recs])) for h in L["RH"]]
            summary[f"{design}|{R}"] = dict(status=st, d=curve, d40=med("d40"), sqrt=med("sqrt"), inv=med("inv"),
                                            chosen=med(est), random=med("random"))
            out(f"  R={R:2d} {design:9s} d = {' '.join(f'{x:.3f}' for x in curve)}   sqrt {med('sqrt'):.3f}"
                f"  inv {med('inv'):.3f}  chosen ({est}) {med(est):.3f}  random {med('random'):.3f}   status: {st}")
        diff = {k: [d[k] - c[k] for d, c in zip(per["disjoint"], per["control"])] for k in ("d40", est, "random")}
        summary[f"difference|{R}"] = {k: dict(median=float(np.median(v)), min=float(np.min(v)), max=float(np.max(v)))
                                      for k, v in diff.items()}
        out(f"  R={R:2d} separate sessions minus control, per subject: "
            + "; ".join(f"{k} median {np.median(v):+.3f} (range {np.min(v):+.3f} to {np.max(v):+.3f})"
                        for k, v in diff.items()))
    open(f"{OUT}/analysis.txt", "w").write("\n".join(lines) + "\n")
    json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
