#!/usr/bin/env python3
"""Evidence for the review audit (docs/review_audit.md), recomputed from the committed rows.

Nothing is simulated and no data file is opened. Every number comes from the result rows of the
confirmatory runs, the comparison of checks and the pilot, through the analysis scripts' own
estimators. The items follow the audit:

E1  the plasmode check at every budget and for both extrapolants (only R = 8 was analyzed);
E2  how far the reference disagreement falls between the shallowest and the deepest reference;
E3  the finite-depth predictive check: fit the two shallower reference sizes, predict the deepest;
E4  a random selection matched to each map's number of passes, and the enrichment over it;
E5  the lower bound L against the order statistics of the subjects' estimates;
E6  the pooled ratio against the mean of per-map ratios, and the maps that pass nothing;
E7  which maps the plasmode simulated;
E8  how often the pilot's certificates resolved a neighborhood in simulation.

Writes results/review_audit/audit.txt and audit_summary.json.

    python3 analyze_review_audit.py
"""
import collections, json, math, os
import numpy as np

import analyze_allen_confirmatory as AA      # its disagreement() and estimates() take RH
import analyze_checks as AC
from reach.checks import NAMES
from run_things_confirmatory import CELLS

OUT = "results/review_audit"
DATA = {
    "EEG": dict(rows="results/things_confirmatory/rows.jsonl",
                sim="results/things_confirmatory/semi_synthetic.jsonl",
                summary="results/things_confirmatory/confirmatory_summary.json",
                checks="results/checks/things", ds="things", RU=(8, 16, 32), RH=(8, 20, 40)),
    "mouse": dict(rows="results/allen_confirmatory/rows.jsonl",
                  sim="results/allen_confirmatory/semi_synthetic.jsonl",
                  summary="results/allen_confirmatory/allen_summary.json",
                  checks="results/checks/allen", ds="allen", RU=(8, 16), RH=(8, 12, 24)),
}
EXT = ("sqrt", "inv")                 # the extrapolants, linear in R_h^-1/2 and in 1/R_h
G = {"sqrt": lambda h: h ** -0.5, "inv": lambda h: 1.0 / h}
CELL_ORDER = [n for n, _ in CELLS]

lines, summary = [], {}


def out(s=""):
    lines.append(s)
    print(s)


def load(path):
    return [json.loads(l) for l in open(path)]


def stats(v):
    """Median, mean, RMSE, min and max of a vector of errors (NaNs dropped)."""
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if not len(v):
        return dict(n=0)
    return dict(n=len(v), median=float(np.median(v)), mean=float(v.mean()),
                rmse=float(np.sqrt((v ** 2).mean())), min=float(v.min()), max=float(v.max()))


def fmt(s, keys=("median", "rmse", "min", "max")):
    return "  ".join(f"{k} {s[k]:+.3f}" for k in keys) if s.get("n") else "n/a"


def per_subject_curves(rows, RH, cell=None):
    """{(subject, cell, R): {Rh: (d_q, d_FC, d_TF)}} pooled within subject, as the analyses pool."""
    rows = [r for r in rows if r["Rh"] != "exact" and (cell is None or r.get("cell", "primary") == cell)]
    return AA.disagreement([dict(r, cell=r.get("cell", "primary")) for r in rows], RH)


def two_point_error(d, RH, ext):
    """Fit the two shallower reference sizes linearly in g(R_h); predicted minus observed at the deepest."""
    x = [G[ext](h) for h in RH]
    b = (d[1] - d[0]) / (x[1] - x[0])
    return d[0] + b * (x[2] - x[0]) - d[2]


# ---------------------------------------------------------------------------------------------
def e1_plasmode(name, P):
    """Each subject's exact FC, q and TF against both extrapolants' estimates, at every budget."""
    sim = load(P["sim"])
    res = {}
    out(f"\n{name}: plasmode, per budget (errors are estimate minus exact; over subjects)")
    for R in P["RU"]:
        est = AA.estimates(per_subject_curves([r for r in sim if r["R"] == R], P["RH"]), P["RH"])
        ex = collections.defaultdict(lambda: np.zeros(5))       # n10, n1, sum of q, sum of TF, draws
        for r in sim:
            if r["R"] == R and r["Rh"] == "exact":
                ex[r["subject"]] += [r["n10"], r["n1"], r["dis"] / r["n"], r["t1"] / r["n"], 1]
        exact = {s: dict(FC=a[0] / a[1], q=a[2] / a[4], TF=a[3] / a[4]) for s, a in ex.items()}
        res[R] = {}
        for q in ("FC", "q", "TF"):
            errs = {e: [est[(s, "primary", R)][q][j] - exact[s][q] for s in sorted(exact)]
                    for j, e in enumerate(EXT)}
            errs["smaller"] = [min(est[(s, "primary", R)][q]) - exact[s][q] for s in sorted(exact)]
            res[R][q] = {e: stats(v) for e, v in errs.items()}
            res[R][q]["exact"] = stats([exact[s][q] for s in sorted(exact)])
        f = res[R]["FC"]
        out(f"  R={R:2d} FC  exact {f['exact']['min']:.3f}-{f['exact']['max']:.3f}")
        for e in ("sqrt", "inv", "smaller"):
            out(f"         {e:7s} {fmt(f[e])}")
        out(f"       q   sqrt {fmt(res[R]['q']['sqrt'], ('median', 'rmse'))}   inv {fmt(res[R]['q']['inv'], ('median', 'rmse'))}")
    return res


def e2_e3_curves(name, P):
    """Per cell and budget: the disagreement at each reference size, its fall, and the finite-depth error."""
    rows = load(P["rows"])
    curves = per_subject_curves(rows, P["RH"])
    res = {}
    out(f"\n{name}: reference disagreement d_FC (medians over subjects), its fall from the shallowest to")
    out("  the deepest reference, subjects whose curve is saturated (deepest >= 0.95 and fall <= 0.02),")
    out("  and the finite-depth error (two shallower sizes predict the deepest), per extrapolant")
    cells = sorted({c for (_, c, _) in curves}, key=CELL_ORDER.index)
    for cell in cells:
        for R in P["RU"]:
            subs = sorted(s for (s, c, r) in curves if c == cell and r == R)
            d = np.array([[curves[(s, cell, R)][h][1] for h in P["RH"]] for s in subs])
            ok = np.isfinite(d).all(axis=1)
            d = d[ok]
            fall = d[:, 0] - d[:, -1]
            sat = int(((d[:, -1] >= 0.95) & (fall <= 0.02)).sum())
            fd = {e: [two_point_error(row, P["RH"], e) for row in d] for e in EXT}
            res[f"{cell}|{R}"] = dict(d=[float(x) for x in np.median(d, axis=0)], fall=float(np.median(fall)),
                                      saturated=sat, n=int(ok.sum()),
                                      fd={e: stats(v) for e, v in fd.items()})
            out(f"  {cell:13s} R={R:2d}  d = " + " ".join(f"{x:.3f}" for x in np.median(d, axis=0))
                + f"  fall {np.median(fall):.3f}  saturated {sat}/{ok.sum()}"
                + f"  fd sqrt {np.median(fd['sqrt']):+.3f}  inv {np.median(fd['inv']):+.3f}")
    # the same finite-depth check in the plasmode, against which its bias can be read
    sim = load(P["sim"])
    out(f"  plasmode (primary cell): finite-depth error, medians over subjects")
    for R in P["RU"]:
        c = per_subject_curves([r for r in sim if r["R"] == R], P["RH"])
        d = np.array([[c[k][h][1] for h in P["RH"]] for k in sorted(c)])
        fd = {e: [two_point_error(row, P["RH"], e) for row in d] for e in EXT}
        res[f"plasmode|{R}"] = dict(d=[float(x) for x in np.median(d, axis=0)],
                                    fd={e: stats(v) for e, v in fd.items()})
        out(f"  plasmode      R={R:2d}  d = " + " ".join(f"{x:.3f}" for x in np.median(d, axis=0))
            + f"  fd sqrt {np.median(fd['sqrt']):+.3f}  inv {np.median(fd['inv']):+.3f}")
    return res


def extrapolate(rate, RH):
    """Both extrapolants' intercepts of a statistic measured at the reference sizes RH."""
    return tuple(AA.intercept(np.array([G[e](h) for h in RH]), np.asarray(rate)) for e in EXT)


def e4_random_cells(name, P):
    """The plug-in label of every cell against a random selection of as many neighborhoods per map."""
    rows = [r for r in load(P["rows"]) if r["Rh"] != "exact"]
    est = AA.estimates(AA.disagreement(rows, P["RH"]), P["RH"])
    acc = collections.defaultdict(lambda: collections.defaultdict(float))
    for r in rows:
        A = acc[(r["subject"], r["cell"], r["R"])]
        A[f"rand{r['Rh']}"] += r["n1"] * (r["nb"] - r["t1"] / r["n"])   # expected T_ref = 0 among n1 random points
        A[f"den{r['Rh']}"] += r["n1"] * r["nb"]
    res = {}
    out(f"\n{name}: plug-in label vs a random selection matched per map (medians over subjects;")
    out("  pairs are sqrt / inv extrapolants); enrichment = precision / random precision")
    for (cell, R) in sorted({(c, r) for (_, c, r) in acc}, key=lambda k: (CELL_ORDER.index(k[0]), k[1])):
        subs = sorted(s for (s, c, r) in acc if c == cell and r == R)
        fc, fr, en = [], [], []
        for s in subs:
            A = acc[(s, cell, R)]
            if not all(A[f"den{h}"] for h in P["RH"]):
                continue
            rnd = extrapolate([A[f"rand{h}"] / A[f"den{h}"] for h in P["RH"]], P["RH"])
            f = est[(s, cell, R)]["FC"]
            fc.append(f)
            fr.append(rnd)
            en.append(tuple((1 - a) / (1 - b) if 1 - b > 0.005 else np.nan for a, b in zip(f, rnd)))
        fc, fr, en = (np.array(x, float) for x in (fc, fr, en))
        res[f"{cell}|{R}"] = dict(FC=np.nanmedian(fc, axis=0).tolist(), FC_random=np.nanmedian(fr, axis=0).tolist(),
                                  enrichment=[float(np.nanmedian(en[:, j])) if np.isfinite(en[:, j]).any() else None
                                              for j in range(2)])
        e = res[f"{cell}|{R}"]["enrichment"]
        out(f"  {cell:13s} R={R:2d}  FC {np.nanmedian(fc[:, 0]):.3f}/{np.nanmedian(fc[:, 1]):.3f}"
            f"  random {np.nanmedian(fr[:, 0]):.3f}/{np.nanmedian(fr[:, 1]):.3f}"
            f"  enrichment " + "/".join("n.d." if x is None else f"{x:.1f}" for x in e))
    return res


def e4_random_checks(name, P):
    """Every check against a random selection of as many neighborhoods per map, at every budget."""
    rows = [r for r in load(f"{P['checks']}/rows.jsonl") if r["Rh"] != "exact"]
    est = AC.estimates(rows, P["RH"])
    acc = collections.defaultdict(lambda: collections.defaultdict(float))
    for r in rows:
        A = acc[(r["subject"], r["R"])]
        for X in NAMES:
            A[f"{X}_a{r['Rh']}"] += r[f"{X}_a"]
            A[f"{X}_rand{r['Rh']}"] += r[f"{X}_pass"] * (r["nb"] - r["t1"] / r["n"])
            A[f"{X}_den{r['Rh']}"] += r[f"{X}_pass"] * r["nb"]
    res, worst = {}, 0.0
    out(f"\n{name}: checks vs a random selection matched per map (medians over subjects with passes;")
    out("  pairs are sqrt / inv); enrichment = precision / random precision")
    for R in P["RU"]:
        for X in NAMES:
            fc, fr, en = [], [], []
            for s in sorted(s for (s, r) in acc if r == R):
                A = acc[(s, R)]
                if not all(A[f"{X}_den{h}"] for h in P["RH"]):
                    continue
                f = extrapolate([A[f"{X}_a{h}"] / A[f"{X}_den{h}"] for h in P["RH"]], P["RH"])
                worst = max(worst, max(abs(a - b) for a, b in zip(f, est[(s, R)][X]["FC"])))
                rnd = extrapolate([A[f"{X}_rand{h}"] / A[f"{X}_den{h}"] for h in P["RH"]], P["RH"])
                fc.append(f)
                fr.append(rnd)
                en.append(tuple((1 - a) / (1 - b) if 1 - b > 0.005 else np.nan for a, b in zip(f, rnd)))
            if not fc:
                res[f"{X}|{R}"] = dict(n=0)
                out(f"  {X:5s} R={R:2d}  no subject passes in every map pool")
                continue
            fc, fr, en = (np.array(x, float) for x in (fc, fr, en))
            e = [float(np.nanmedian(en[:, j])) if np.isfinite(en[:, j]).any() else None for j in range(2)]
            res[f"{X}|{R}"] = dict(n=len(fc), FC=np.median(fc, axis=0).tolist(),
                                   FC_random=np.median(fr, axis=0).tolist(), enrichment=e)
            out(f"  {X:5s} R={R:2d}  FC {np.median(fc[:, 0]):.3f}/{np.median(fc[:, 1]):.3f}"
                f"  random {np.median(fr[:, 0]):.3f}/{np.median(fr[:, 1]):.3f}"
                f"  enrichment " + "/".join("n.d." if x is None else f"{x:.1f}" for x in e)
                + f"  ({len(fc)} subjects)")
    assert worst < 1e-9, f"recomputed check FC differs from analyze_checks.py by {worst}"
    out(f"  (each check's FC equals analyze_checks.py's to within {worst:.1e})")
    return res


def tail(k, n=10):
    """P(X_(k) > median) for n i.i.d. continuous draws: the one-sided error of X_(k) as a lower bound."""
    return sum(math.comb(n, j) for j in range(k)) / 2 ** n


def e5_order_statistics(name, P):
    """The published L of each claim against the order statistics of the subjects' smaller estimates."""
    rows = [r for r in load(P["rows"]) if r["Rh"] != "exact"]
    est = AA.estimates(AA.disagreement(rows, P["RH"]), P["RH"])
    S = json.load(open(P["summary"]))["summary"]
    claims = [("C1", "primary", 8, "FC"), ("C1-RSA", "rsa_cv_mds", 8, "FC")]
    if 32 in P["RU"]:
        claims.append(("C2", "primary", 32, "FC"))
    claims += [(f"C3 (R={R})", "primary", R, "q") for R in P["RU"]]
    res = {}
    out(f"\n{name}: the lower bound L and the order statistics x(k) of the subjects' smaller estimates")
    out(f"  (as a lower bound for a population median from 10 subjects, x(2) errs with probability"
        f" {tail(2):.4f}, x(3) with {tail(3):.4f}; the nominal one-sided level is 0.025)")
    for lab, cell, R, q in claims:
        v = np.sort([min(est[(s, cell, R)][q]) for s in sorted(s for (s, c, r) in est if c == cell and r == R)])
        L = S[f"{cell}|{R}|{q}"]["L"]
        k = int(np.argmin(np.abs(v - L))) + 1
        res[lab] = dict(L=L, x1=float(v[0]), x2=float(v[1]), x3=float(v[2]), L_is_order_statistic=k
                        if abs(v[k - 1] - L) < 1e-9 else None)
        out(f"  {lab:12s} L {L:.3f}   x(1) {v[0]:.3f}  x(2) {v[1]:.3f}  x(3) {v[2]:.3f}"
            + (f"   L = x({k})" if abs(v[k - 1] - L) < 1e-9 else ""))
    return res


def e5_checks(name, P):
    """The checks' L and U at R = 8 against the order statistics of the smaller and larger estimates."""
    rows = [r for r in load(f"{P['checks']}/rows.jsonl") if r["Rh"] != "exact"]
    est = AC.estimates(rows, P["RH"])
    S = json.load(open("results/checks/checks_summary.json"))["datasets"][P["ds"]]["summary"]
    res = {}
    out(f"  checks at R = 8: L against x(2), x(3) of the smaller estimates; U against x(8), x(9) of the larger")
    for X in NAMES[1:]:
        lo = np.sort([min(est[(s, 8)][X]["FC"]) for s in sorted(s for (s, r) in est if r == 8)
                      if np.isfinite(est[(s, 8)][X]["FC"]).all()])
        hi = np.sort([max(est[(s, 8)][X]["FC"]) for s in sorted(s for (s, r) in est if r == 8)
                      if np.isfinite(est[(s, 8)][X]["FC"]).all()])
        f = S[f"{X}|8|FC"]
        if len(lo) < 10:
            out(f"  {X:5s} {len(lo)} of 10 subjects with passes: order statistics not compared")
            res[X] = dict(n=len(lo))
            continue
        res[X] = dict(n=len(lo), L=f["L"], U=f["U"], lo_x2=float(lo[1]), lo_x3=float(lo[2]),
                      hi_x8=float(hi[-3]), hi_x9=float(hi[-2]))
        out(f"  {X:5s} L {f['L']:.3f} (x(2) {lo[1]:.3f}, x(3) {lo[2]:.3f})   U {f['U']:.3f}"
            f" (x(8) {hi[-3]:.3f}, x(9) {hi[-2]:.3f})")
    return res


def e6_pooling(name, P):
    """At R = 8 and the deepest reference: pooled ratio vs mean of per-map ratios; maps passing nothing."""
    h = P["RH"][-1]
    rows = [r for r in load(f"{P['checks']}/rows.jsonl") if r["R"] == 8 and r["Rh"] == h]
    res = {}
    out(f"\n{name}: R = 8, R_h = {h}: pooled ratio vs mean of per-map ratios (medians over subjects),"
        " largest difference in a subject, and the share of maps passing nothing")
    by = collections.defaultdict(list)
    for r in rows:
        by[r["subject"]].append(r)
    for X in NAMES:
        pooled, mean, diff, zero = [], [], [], []
        for s, rs in sorted(by.items()):
            p = np.array([r[f"{X}_pass"] for r in rs], float)
            a = np.array([r[f"{X}_a"] for r in rs], float)       # nb = 1 at the deepest reference
            zero.append(float((p == 0).mean()))
            if p.sum():
                pooled.append(a.sum() / p.sum())
                mean.append(float(np.mean(a[p > 0] / p[p > 0])))
                diff.append(abs(pooled[-1] - mean[-1]))
        res[X] = dict(pooled=float(np.median(pooled)) if pooled else None,
                      mean_of_maps=float(np.median(mean)) if mean else None,
                      max_abs_diff=float(max(diff)) if diff else None,
                      zero_pass_share_median=float(np.median(zero)), zero_pass_share_max=float(max(zero)))
        out(f"  {X:5s} " + (f"pooled {np.median(pooled):.3f}  per-map {np.median(mean):.3f}  max diff {max(diff):.3f}"
                            if pooled else "no passes                                       ")
            + f"   maps passing nothing: median {np.median(zero):.2f}, max {max(zero):.2f}")
    return res


def e7_plasmode_scope():
    """Which maps the plasmode simulated: its rows carry no cell, so only the primary t-SNE map."""
    paths = [DATA[k]["sim"] for k in DATA] + [f"{DATA[k]['checks']}/semi_synthetic.jsonl" for k in DATA]
    cells = {p: sorted({r.get("cell", "primary (no cell field)") for r in load(p)}) for p in paths}
    out("\nE7: maps simulated by the plasmode checks")
    for p, c in cells.items():
        out(f"  {p}: {', '.join(c)}")
    return cells


def e8_pilot():
    """How many of the pilot's simulated draws resolved any neighborhood, at each tolerance."""
    rows = load("results/things_pilot/semi_synthetic.jsonl")
    need = [r for r in rows if r["kind"] == "need"]
    valid = [r for r in rows if r["kind"] == "valid"]
    res = dict(need_draws=len(need), validity_draws=len(valid),
               resolved_any={k: int(sum(r[k] > 0 for r in need)) for k in ("resolved_10", "resolved_20", "resolved_40")},
               false_any={k: int(sum(r[k] > 0 for r in need)) for k in ("false_10", "false_20", "false_40")},
               pair_violations=int(sum(r["pair_viol"] for r in valid)),
               reach_violations=int(sum(r["reach_viol"] for r in valid)),
               validity_rows_record=sorted(set(valid[0]) - {"kind", "subject", "draw"}))
    out("\nE8: the pilot's certificates in simulation (results/things_pilot/semi_synthetic.jsonl)")
    out(f"  {len(need)} draws recorded issuance: draws resolving any neighborhood at k' = 10, 20, 40: "
        + ", ".join(str(v) for v in res["resolved_any"].values()))
    out(f"  {len(valid)} validity draws record only {', '.join(res['validity_rows_record'])}:"
        f" {res['pair_violations']} pairwise and {res['reach_violations']} reach violations")
    return res


def main():
    out("REVIEW AUDIT: evidence recomputed from the committed rows (analyze_review_audit.py)")
    for name, P in DATA.items():
        summary[name] = dict(E1=e1_plasmode(name, P), E2_E3=e2_e3_curves(name, P),
                             E4_cells=e4_random_cells(name, P), E4_checks=e4_random_checks(name, P),
                             E5=e5_order_statistics(name, P), E5_checks=e5_checks(name, P),
                             E6=e6_pooling(name, P))
    summary["E7"] = e7_plasmode_scope()
    summary["E8"] = e8_pilot()
    os.makedirs(OUT, exist_ok=True)
    open(f"{OUT}/audit.txt", "w").write("\n".join(lines) + "\n")
    json.dump(summary, open(f"{OUT}/audit_summary.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
