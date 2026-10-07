#!/usr/bin/env python3
"""Layer C's measures (docs/revision_protocol.md, Sections 3.4 and 3.6), from run_study_calibration.py.

- The truth: each virtual subject's theta_v is pooled over its oracle draws, with a Monte Carlo
  error. The population target Theta is the median of the 30 theta_v; the sample target is the
  median of a study's 10 drawn theta_v.
- Each study's published procedure:
  - per-subject estimates from both extrapolants, and their medians;
  - 10,000 bootstrap resamples of subjects, giving L (the 2.5th percentile of the smaller median) and
    U (the 97.5th percentile of the larger);
  - the guard, from the plasmode of the first 200 studies: L is lowered by the median bias of the
    smaller estimate when that exceeds +0.05.
- The alternative: the order-statistic bounds, the 2nd-smallest smaller estimate and the 2nd-largest
  larger one.
- Reported with Monte Carlo errors: Pr(L > Theta), Pr(U < Theta), coverage of [L, U], width, the rate
  of declaring L > 0.5, and the rate of declaring it erroneously.

    python3 analyses/calibration/analyze_study_calibration.py --dataset EEG --scenario T1 --map tsne
"""
import argparse, collections, gzip, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from analyses.validation import analyze_validation as AV
from analyses.calibration import run_study_calibration as SC
from reach import sim

B = 10_000


def counts(d):
    return {(k.split("|")[0], int(k.split("|")[1])): np.array(v, float) for k, v in d.items()}


def rate(x):
    x = np.asarray(x, float)
    p = float(x.mean())
    return dict(p=p, mcse=float(np.sqrt(p * (1 - p) / len(x))), n=int(len(x)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="EEG")
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--map", default="tsne")
    ap.add_argument("--R", type=int, default=SC.R)
    a = ap.parse_args()
    tag = f"{a.dataset}_{a.scenario}_{a.map}_R{a.R}"
    L = sim.LAYOUT[a.dataset]
    # the truth
    orc = collections.defaultdict(list)
    for line in gzip.open(f"{SC.OUT}/oracle_{tag}.jsonl.gz", "rt"):
        r = json.loads(line)
        orc[r["v"]].append((r["n10"], r["n1"], r["dis"], r["n"]))
    theta, theta_se, theta_q = {}, {}, {}
    for v, rows in orc.items():
        n10, n1, dis, nn = np.array(rows, float).T
        t = n10.sum() / n1.sum()
        theta[v] = t
        theta_se[v] = float(np.sqrt(((n10 - t * n1) ** 2).sum() / (len(n1) * (len(n1) - 1))) / n1.mean())
        theta_q[v] = dis.sum() / nn.sum()
    Theta = float(np.median(list(theta.values())))
    Theta_q = float(np.median(list(theta_q.values())))
    # the studies
    st = collections.defaultdict(list)
    for line in gzip.open(f"{SC.OUT}/studies_{tag}.jsonl.gz", "rt"):
        r = json.loads(line)
        st[r["study"]].append(r)
    res = []
    for study, rows in sorted(st.items()):
        rows.sort(key=lambda r: r["slot"])
        recs = [AV.record(counts(r["counts"]), L) for r in rows]
        e1 = np.array([x["sqrt"] for x in recs])
        e2 = np.array([x["inv"] for x in recs])
        rng = np.random.default_rng((SC.STUDY_SEED, 0 if a.dataset == "EEG" else 1, SC.SCEN[a.scenario],
                                     SC.MAPS[a.map], a.R, study, 99))
        idx = rng.integers(0, len(recs), size=(B, len(recs)))
        m1, m2 = np.nanmedian(e1[idx], axis=1), np.nanmedian(e2[idx], axis=1)
        Lb, Ub = float(np.percentile(np.minimum(m1, m2), 2.5)), float(np.percentile(np.maximum(m1, m2), 97.5))
        guard = None
        if all(r["plasmode"] for r in rows):
            bias = []
            for r in rows:
                p = AV.record(counts(r["plasmode"]), L)
                bias.append(min(p["sqrt"], p["inv"]) - p["exact"])
            guard = float(np.median(bias))
            if guard > 0.05:
                Lb -= guard
        q1, q2 = np.array([x["q_sqrt"] for x in recs]), np.array([x["q_inv"] for x in recs])
        mq1, mq2 = np.nanmedian(q1[idx], axis=1), np.nanmedian(q2[idx], axis=1)
        Lq = float(np.percentile(np.minimum(mq1, mq2), 2.5))
        small = np.sort(np.minimum(e1, e2))
        large = np.sort(np.maximum(e1, e2))
        sample = float(np.median([theta[r["v"]] for r in rows]))
        res.append(dict(study=study, L=Lb, U=Ub, x2=float(small[1]), x9=float(large[-2]), guard=guard, Lq=Lq,
                        sample_q=float(np.median([theta_q[r["v"]] for r in rows])),
                        med_sqrt=float(np.median(e1)), med_inv=float(np.median(e2)), sample=sample))
    Ls, Us = np.array([r["L"] for r in res]), np.array([r["U"] for r in res])
    x2, x9 = np.array([r["x2"] for r in res]), np.array([r["x9"] for r in res])
    smp = np.array([r["sample"] for r in res])
    guards = [r["guard"] for r in res if r["guard"] is not None]
    out = dict(dataset=a.dataset, scenario=a.scenario, map=a.map, Theta=Theta, theta=theta, theta_se=theta_se,
               studies=len(res),
               published=dict(L_above=rate(Ls > Theta), U_below=rate(Us < Theta), coverage=rate((Ls <= Theta) & (Theta <= Us)),
                              width=float(np.mean(Us - Ls)), declares_above_half=rate(Ls > 0.5),
                              erroneous_above_half=rate((Ls > 0.5) & (Theta <= 0.5)),
                              L_above_sample=rate(Ls > smp), U_below_sample=rate(Us < smp)),
               order_statistics=dict(L_above=rate(x2 > Theta), U_below=rate(x9 < Theta),
                                     coverage=rate((x2 <= Theta) & (Theta <= x9)), width=float(np.mean(x9 - x2)),
                                     declares_above_half=rate(x2 > 0.5), L_above_sample=rate(x2 > smp)),
               mislabel=dict(Theta=Theta_q, L_above=rate(np.array([r["Lq"] for r in res]) > Theta_q),
                             declares_above_tenth=rate(np.array([r["Lq"] for r in res]) > 0.10),
                             L_above_sample=rate(np.array([r["Lq"] for r in res]) > np.array([r["sample_q"] for r in res]))),
               guard=dict(studies=len(guards), triggered=int(sum(g > 0.05 for g in guards)),
                          median_bias=float(np.median(guards)) if guards else None),
               medians=dict(sqrt=float(np.median([r["med_sqrt"] for r in res])),
                            inv=float(np.median([r["med_inv"] for r in res]))))
    lines = [f"LAYER C ({a.dataset}, truth {a.scenario}, {a.map} maps, R = {a.R}): {len(res)} virtual studies of"
             f" {SC.SUBJECTS} subjects from a population of {len(theta)}",
             f"  population target Theta = {Theta:.3f} (theta_v {min(theta.values()):.3f}-{max(theta.values()):.3f},"
             f" Monte Carlo error per theta_v up to {max(theta_se.values()):.3f})",
             f"  study medians: sqrt {out['medians']['sqrt']:.3f}, inv {out['medians']['inv']:.3f}"]
    for name in ("published", "order_statistics"):
        o = out[name]
        lines.append(f"  {name:16s} Pr(L > Theta) {o['L_above']['p']:.4f} +- {o['L_above']['mcse']:.4f}"
                     f"  Pr(U < Theta) {o['U_below']['p']:.4f} +- {o['U_below']['mcse']:.4f}"
                     f"  coverage {o['coverage']['p']:.3f}  width {o['width']:.3f}"
                     f"  declares L > 0.5 in {o['declares_above_half']['p']:.3f}"
                     f"  Pr(L > sample median) {o['L_above_sample']['p']:.4f}")
    m = out["mislabel"]
    lines.append(f"  mislabel rate q: Theta {m['Theta']:.3f}; Pr(L > Theta) {m['L_above']['p']:.4f} +- {m['L_above']['mcse']:.4f};"
                 f" declares L > 0.10 in {m['declares_above_tenth']['p']:.3f}")
    g = out["guard"]
    lines.append(f"  guard (first {g['studies']} studies): triggered in {g['triggered']}; median plasmode bias"
                 f" {g['median_bias'] if g['median_bias'] is None else round(g['median_bias'], 4)}")
    lines.append(f"  criterion (addendum): Pr(L > Theta) <= 0.025 + 2 MCSE = "
                 f"{0.025 + 2 * out['published']['L_above']['mcse']:.4f} -> "
                 + ("met" if out["published"]["L_above"]["p"] <= 0.025 + 2 * out["published"]["L_above"]["mcse"] else "NOT met"))
    print("\n".join(lines))
    open(f"{SC.OUT}/calibration_{tag}.txt", "w").write("\n".join(lines) + "\n")
    json.dump(out, open(f"{SC.OUT}/calibration_{tag}.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
