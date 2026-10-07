#!/usr/bin/env python3
"""Apply docs/checks_protocol.md (plan Amendment 5), Sections 4 to 6, to both datasets.

Reads results/checks/{things,allen}/{rows.jsonl, semi_synthetic.jsonl} and the confirmatory
rows the maps must reproduce; writes results/checks/checks_analysis.txt and
checks_summary.json.
"""
import collections, json, os, sys, warnings
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from reach.checks import NAMES
from analyses.checks.run_checks import DATASETS

OUT = "results/checks"
B, BOOT_SEED = 10_000, {"things": 5_300_000, "allen": 5_700_000}
CONFIRMATORY = {"things": "results/things_confirmatory/rows.jsonl",
                "allen": "results/allen_confirmatory/rows.jsonl"}
LABEL = {"things": "THINGS-EEG2, human EEG", "allen": "Allen two-photon, mouse V1"}
ALTERNATIVES = NAMES[1:]
MIN_SUBJECTS, CLAIM_R, GUARD = 5, 8, 0.05
QUANTITIES = ("FC", "CATCH", "YIELD", "FC_NX")


def intercept(x, y):
    A = np.column_stack([np.ones_like(x), x])
    return float(np.linalg.lstsq(A, y, rcond=None)[0][0])


def estimates(rows, RH):
    """Per (subject, R) and check: PASS and PASS_NX, and both extrapolants' FC, CATCH, YIELD
    and FC_NX, the false comfort of X used on top of N (Section 4)."""
    acc = collections.defaultdict(lambda: collections.defaultdict(float))
    for r in rows:
        if r["Rh"] == "exact":
            continue
        A, h = acc[(r["subject"], r["R"])], r["Rh"]
        if h == RH[0]:                               # every display counted once for PASS
            A["n"] += r["n"]
            for X in NAMES:
                A[f"{X}_pass"] += r[f"{X}_pass"]
                A[f"{X}_nx"] += r[f"{X}_nx"]
        A[f"nn{h}"] += r["n"] * r["nb"]
        A[f"c{h}"] += r["n10"]
        for X in NAMES:
            A[f"{X}_a{h}"] += r[f"{X}_a"]
            A[f"{X}_b{h}"] += r[f"{X}_b"]
    xs = (np.array([h ** -0.5 for h in RH]), np.array([1.0 / h for h in RH]))
    est = {}
    for key, A in acc.items():
        rate = lambda name: np.array([A[f"{name}{h}"] / A[f"nn{h}"] for h in RH])
        c_hat = [intercept(x, rate("c")) for x in xs]
        e = {}
        for X in NAMES:
            p, pnx = A[f"{X}_pass"] / A["n"], A[f"{X}_nx"] / A["n"]
            a_hat = [intercept(x, rate(f"{X}_a")) for x in xs]
            b_hat = [intercept(x, rate(f"{X}_b")) for x in xs]
            e[X] = dict(PASS=p, PASS_NX=pnx,
                        FC=tuple(a / p if p > 0 else np.nan for a in a_hat),
                        CATCH=tuple(b / c if c > 0 else np.nan for b, c in zip(b_hat, c_hat)),
                        YIELD=tuple(p - a for a in a_hat),
                        FC_NX=tuple((c - b) / pnx if pnx > 0 else np.nan for b, c in zip(b_hat, c_hat)))
        est[key] = e
    return est


def exact(rows):
    """Per subject and check: the exact FC from the rows counted against a known truth."""
    acc = collections.defaultdict(lambda: np.zeros(2))
    for r in rows:
        if r["Rh"] == "exact":
            for X in NAMES:
                acc[(r["subject"], X)] += [r[f"{X}_a"], r[f"{X}_pass"]]
    return {key: (a / p if p else np.nan) for key, (a, p) in acc.items()}


def summarize(est, R, X, q, rng):
    """Medians over subjects, per-extrapolant 95% intervals, L and U (Section 4)."""
    subs = sorted(s for (s, r) in est if r == R)
    e1 = np.array([est[(s, R)][X][q][0] for s in subs])
    e2 = np.array([est[(s, R)][X][q][1] for s in subs])
    idx = rng.integers(0, len(subs), size=(B, len(subs)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m1, m2 = np.nanmedian(e1[idx], axis=1), np.nanmedian(e2[idx], axis=1)
        return dict(med_sqrt=float(np.nanmedian(e1)), med_inv=float(np.nanmedian(e2)),
                    ci_sqrt=[float(np.nanpercentile(m1, 2.5)), float(np.nanpercentile(m1, 97.5))],
                    ci_inv=[float(np.nanpercentile(m2, 2.5)), float(np.nanpercentile(m2, 97.5))],
                    L=float(np.nanpercentile(np.minimum(m1, m2), 2.5)),
                    U=float(np.nanpercentile(np.maximum(m1, m2), 97.5)),
                    n_defined=int(np.isfinite(e1).sum()), n_subjects=len(subs))


def pass_summary(est, R, X, rng, q="PASS"):
    subs = sorted(s for (s, r) in est if r == R)
    v = np.array([est[(s, R)][X][q] for s in subs])
    m = np.median(v[rng.integers(0, len(subs), size=(B, len(subs)))], axis=1)
    return dict(med=float(np.median(v)), ci=[float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))],
                n_passing=int((v > 0).sum()))


def reproduction(ds, rows):
    """Plug-in counts of every re-created map against the committed confirmatory rows."""
    ref = {}
    for line in open(CONFIRMATORY[ds]):
        r = json.loads(line)
        if r["cell"] == "primary":
            ref[(r["subject"], r["split"], r["R"], r["seed"], r["Rh"])] = (r["n1"], r["n10"], r["dis"], r["t1"])
    mism = [k for k in ((r["subject"], r["split"], r["R"], r["seed"], r["Rh"]) for r in rows)
            if k not in ref]
    new = {(r["subject"], r["split"], r["R"], r["seed"], r["Rh"]): (r["n1"], r["n10"], r["dis"], r["t1"])
           for r in rows}
    diff = [k for k in new if k in ref and new[k] != ref[k]]
    return dict(compared=len(new) - len(mism), missing=len(mism), differing=len(diff),
                missing_in_new=len(set(ref) - set(new)), first_differences=[list(k) for k in diff[:5]])


def semi_synthetic_bias(sim, RH):
    """Section 5: median over subjects of (smaller - exact) and (larger - exact) FC, per check."""
    est, ex = estimates(sim, RH), exact(sim)
    out = {}
    for X in NAMES:
        lo, hi, rows = [], [], []
        for (s, R), e in sorted(est.items()):
            x = ex.get((s, X), np.nan)
            if np.isfinite(x) and np.all(np.isfinite(e[X]["FC"])):
                lo.append(min(e[X]["FC"]) - x)
                hi.append(max(e[X]["FC"]) - x)
                rows.append(dict(subject=s, exact=x, est=list(e[X]["FC"]), PASS=e[X]["PASS"]))
        out[X] = dict(bias_lo=float(np.median(lo)) if lo else np.nan,
                      bias_hi=float(np.median(hi)) if hi else np.nan,
                      n=len(lo), per_subject=rows)
    return out


def classify(S, X):
    """Section 6: the classification of check X at the claim budget, bias guard applied."""
    f = S["FC"][X]
    if f["n_defined"] < MIN_SUBJECTS:
        return "not evaluable", f["L"], f["U"]
    L, U = f["L"] - S["guard_L"][X], f["U"] + S["guard_U"][X]
    if L > 0.5:
        return "leaves most false comfort", L, U
    if U < 0.5:
        return "removes most false comfort", L, U
    return "undetermined", L, U


def main():
    lines, summary = [], {}

    def out(s=""):
        lines.append(s)
        print(s)

    out("COMMON CHECKS AGAINST THE TRUTH (docs/checks_protocol.md): primary cell; medians over")
    out("subjects; pairs are (linear in R_h^-1/2) / (linear in 1/R_h); L, U = conservative bounds")
    classes = {}
    for ds in ("things", "allen"):
        P = DATASETS[ds]
        rows = [json.loads(l) for l in open(f"{P['out']}/rows.jsonl")]
        sim = [json.loads(l) for l in open(f"{P['out']}/semi_synthetic.jsonl")]
        rng = np.random.default_rng(BOOT_SEED[ds])
        rep = reproduction(ds, rows)
        bias = semi_synthetic_bias(sim, P["RH"])
        est = estimates(rows, P["RH"])
        out(f"\n=== {LABEL[ds]} ===")
        out(f"Reproduction of the confirmatory maps: {rep['compared']} rows compared, "
            f"{rep['differing']} differ, {rep['missing']} missing from the confirmatory rows, "
            f"{rep['missing_in_new']} confirmatory rows not re-created")
        out("Estimator check (Section 5), R=8: median bias of FC, smaller / larger estimate")
        guard_L, guard_U = {}, {}
        for X in NAMES:
            b = bias[X]
            guard_L[X] = b["bias_lo"] if (np.isfinite(b["bias_lo"]) and b["bias_lo"] > GUARD) else 0.0
            guard_U[X] = -b["bias_hi"] if (np.isfinite(b["bias_hi"]) and b["bias_hi"] < -GUARD) else 0.0
            flag = ("  guard: L lowered by %.3f" % guard_L[X] if guard_L[X] else "") + \
                   ("  guard: U raised by %.3f" % guard_U[X] if guard_U[X] else "")
            out(f"  {X:5s} {b['bias_lo']:+.3f} / {b['bias_hi']:+.3f}  ({b['n']} subjects){flag}")
        S = {}
        out("\ncheck  R   pass    FC (sqrt/inv)    L      U      catch (sqrt/inv)  yield (sqrt/inv)"
            "   with N: pass  FC (sqrt/inv)")
        for X in NAMES:
            for R in P["RU"]:
                ps, pnx = pass_summary(est, R, X, rng), pass_summary(est, R, X, rng, "PASS_NX")
                q = {qn: summarize(est, R, X, qn, rng) for qn in QUANTITIES}
                S[(X, R)] = dict(PASS=ps, PASS_NX=pnx, **q)
                f, c, y, fn = q["FC"], q["CATCH"], q["YIELD"], q["FC_NX"]
                cat = "   -   " if X == "N" else f"{c['med_sqrt']:.3f}/{c['med_inv']:.3f}"
                out(f"{X:5s} {R:2d}  {ps['med']:.3f}   {f['med_sqrt']:.3f}/{f['med_inv']:.3f}   "
                    f"{f['L']:.3f}  {f['U']:.3f}   {cat:15s}   {y['med_sqrt']:.3f}/{y['med_inv']:.3f}"
                    f"      {pnx['med']:.3f}  {fn['med_sqrt']:.3f}/{fn['med_inv']:.3f}"
                    + (f"   [FC defined in {f['n_defined']} of {f['n_subjects']} subjects]"
                       if f["n_defined"] < f["n_subjects"] else ""))
        claim = dict(FC={X: S[(X, CLAIM_R)]["FC"] for X in NAMES}, guard_L=guard_L, guard_U=guard_U)
        out(f"\nClassification at R = {CLAIM_R} (Section 6):")
        classes[ds] = {}
        for X in ALTERNATIVES:
            cl, L, U = classify(claim, X)
            classes[ds][X] = cl
            out(f"  {X:5s} L = {L:.3f}, U = {U:.3f} -> {cl}")
        summary[ds] = dict(reproduction=rep, bias={X: {k: v for k, v in bias[X].items()} for X in NAMES},
                           guard_L=guard_L, guard_U=guard_U, classification=classes[ds],
                           summary={f"{X}|{R}|{qn}": v for (X, R), d in S.items() for qn, v in d.items()})
    removes = [(ds, X) for ds in classes for X, c in classes[ds].items() if c == "removes most false comfort"]
    leaves_all = [X for X in ALTERNATIVES if all(classes[ds][X] == "leaves most false comfort" for ds in classes)]
    out("\nWhat the manuscript may say (Section 6):")
    if not removes:
        out("  - none of the common checks removed most false comfort in either dataset")
    else:
        out("  - checks that removed most false comfort, reported with pass share and yield: "
            + ", ".join(f"{X} ({ds})" for ds, X in removes))
    if len(leaves_all) == len(ALTERNATIVES):
        out("  - all six checks left most false comfort in both datasets")
    else:
        out("  - checks that left most false comfort in both datasets: " + (", ".join(leaves_all) or "none"))
    open(f"{OUT}/checks_analysis.txt", "w").write("\n".join(lines) + "\n")
    json.dump(dict(classification=classes, removes=removes, leaves_in_both=leaves_all, datasets=summary),
              open(f"{OUT}/checks_summary.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
