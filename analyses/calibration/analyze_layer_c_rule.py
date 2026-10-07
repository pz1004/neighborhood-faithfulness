#!/usr/bin/env python3
"""Why the published bound met or missed its level in Layer C: the bootstrap rule against the estimator (descriptive).

Added after the first two Layer C runs were read; it is not part of the protocol's criterion. For every
Layer C run in results/revision/layer_c/, it repeats the published procedure of
analyze_study_calibration.py with the same bootstrap resamples, and also applies that rule to the
sampled subjects' true theta_v instead of their estimates. That separates the rule's own error rate
from the estimator's contribution. It also reports how often L equals x(3), the 3rd-smallest of the
10 subjects' smaller estimates, which audit item E5 took L to be.
Quantities: false comfort, and the mislabel rate q. The guard is left out (it never triggered).

First, the rule's own error rate as a lower bound for a population median, on three continuous
populations of known median (20,000 studies of 10 each; seed range 9,200,000, unused before).

Writes results/revision/layer_c/rule.txt and rule.json.
"""
import collections, glob, gzip, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from analyses.calibration import analyze_study_calibration as AS
from analyses.validation import analyze_validation as AV
from analyses.calibration import run_study_calibration as SC
from reach import sim


SEED, KNOWN_STUDIES = 9_200_000, 20_000


def rate(x):
    return dict(p=float(np.mean(x)), mcse=float(np.sqrt(np.mean(x) * (1 - np.mean(x)) / len(x))))


def known_populations():
    """Pr(L > median) of the published rule (2.5th percentile of B bootstrap medians of 10) on known populations."""
    out = {}
    for i, (name, draw, med) in enumerate([("normal", lambda g: g.normal(size=10), 0.0),
                                           ("uniform", lambda g: g.uniform(size=10), 0.5),
                                           ("exponential", lambda g: g.exponential(size=10), float(np.log(2)))]):
        g = np.random.default_rng((SEED, i))
        above = []
        for _ in range(KNOWN_STUDIES):
            x = draw(g)
            above.append(np.percentile(np.median(x[g.integers(0, 10, size=(AS.B, 10))], axis=1), 2.5) > med)
        out[name] = rate(np.array(above))
    return out


def run(tag):
    ds, scen, mp, R = tag.split("_")
    R = int(R[1:])
    L = sim.LAYOUT[ds]
    orc = collections.defaultdict(list)
    for line in gzip.open(f"{SC.OUT}/oracle_{tag}.jsonl.gz", "rt"):
        r = json.loads(line)
        orc[r["v"]].append((r["n10"], r["n1"], r["dis"], r["n"]))
    sums = {v: np.array(x, float).sum(axis=0) for v, x in orc.items()}
    theta = {v: s[0] / s[1] for v, s in sums.items()}
    theta_q = {v: s[2] / s[3] for v, s in sums.items()}
    target = dict(FC=float(np.median(list(theta.values()))), q=float(np.median(list(theta_q.values()))))
    st = collections.defaultdict(list)
    for line in gzip.open(f"{SC.OUT}/studies_{tag}.jsonl.gz", "rt"):
        r = json.loads(line)
        st[r["study"]].append(r)
    acc = collections.defaultdict(list)
    for study, rows in sorted(st.items()):
        rows.sort(key=lambda r: r["slot"])
        recs = [AV.record(AS.counts(r["counts"]), L) for r in rows]
        rng = np.random.default_rng((SC.STUDY_SEED, 0 if ds == "EEG" else 1, SC.SCEN[scen], SC.MAPS[mp], R, study, 99))
        idx = rng.integers(0, len(recs), size=(AS.B, len(recs)))
        for qn, (k1, k2), truth in (("FC", ("sqrt", "inv"), theta), ("q", ("q_sqrt", "q_inv"), theta_q)):
            e1, e2 = np.array([x[k1] for x in recs]), np.array([x[k2] for x in recs])
            tv = np.array([truth[r["v"]] for r in rows])
            small = np.minimum(e1, e2)
            Lb = np.percentile(np.minimum(np.nanmedian(e1[idx], axis=1), np.nanmedian(e2[idx], axis=1)), 2.5)
            x3 = np.sort(small[np.isfinite(small)])[2]
            acc[qn].append(dict(L=Lb, L_exact=np.percentile(np.median(tv[idx], axis=1), 2.5), x3=x3,
                                x3_exact=np.sort(tv)[2], bias=float(np.nanmedian(small - tv))))
    res = {}
    for qn, a in acc.items():
        g = {k: np.array([x[k] for x in a]) for k in a[0]}
        T = target[qn]
        res[qn] = dict(target=T, L_above=rate(g["L"] > T), L_exact_above=rate(g["L_exact"] > T),
                       x3_above=rate(g["x3"] > T), x3_exact_above=rate(g["x3_exact"] > T),
                       L_is_x3=float(np.mean(np.isclose(g["L"], g["x3"]))),
                       L_below_x3=float(np.mean(g["L"] < g["x3"] - 1e-12)), median_bias=float(np.median(g["bias"])),
                       studies=len(a))
    return res


def main():
    tags = sorted(os.path.basename(p)[len("calibration_"):-len(".json")]
                  for p in glob.glob(f"{SC.OUT}/calibration_*.json"))
    lines = ["LAYER C, the published bound's rule against the estimator (analyze_layer_c_rule.py; descriptive)",
             "Pr(bound > target) for the published L on the estimates and on the exact theta_v, and for x(3);",
             "with exact values, x(3) of 10 errs with probability 0.0547 for a population median"]
    out = {"known_populations": known_populations()}
    lines.append("\nthe rule on known populations, Pr(L > median): " + ", ".join(
        f"{k} {v['p']:.4f} +- {v['mcse']:.4f}" for k, v in out["known_populations"].items()))
    for tag in tags:
        out[tag] = res = run(tag)
        lines.append(f"\n{tag}")
        for qn, r in res.items():
            f = lambda k: f"{r[k]['p']:.4f} +- {r[k]['mcse']:.4f}"
            lines.append(f"  {qn:2s} target {r['target']:.3f}: L on estimates {f('L_above')}, on exact values"
                         f" {f('L_exact_above')}; x(3) {f('x3_above')}, exact {f('x3_exact_above')};"
                         f" L = x(3) in {r['L_is_x3']:.3f}, L < x(3) in {r['L_below_x3']:.3f};"
                         f" median bias of the smaller estimate {r['median_bias']:+.4f}")
    print("\n".join(lines))
    open(f"{SC.OUT}/rule.txt", "w").write("\n".join(lines) + "\n")
    json.dump(out, open(f"{SC.OUT}/rule.json", "w"), indent=1)


if __name__ == "__main__":
    main()
