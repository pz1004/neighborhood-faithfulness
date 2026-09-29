#!/usr/bin/env python3
"""Apply docs/note_protocol.md (plan Amendment 3), Sections 4 to 6, to the confirmatory rows.

Reads results/things_confirmatory/{rows.jsonl, semi_synthetic.jsonl}; writes
results/things_confirmatory/confirmatory_analysis.txt and confirmatory_summary.json.
"""
import collections, json
import numpy as np

from run_things_confirmatory import CELLS

OUT = "results/things_confirmatory"
RU, RH = (8, 16, 32), (8, 20, 40)
B, BOOT_SEED = 10_000, 4_300_000
QUANTITIES = ("FC", "q", "TF")


def intercept(x, y):
    A = np.column_stack([np.ones_like(x), x])
    return float(np.linalg.lstsq(A, y, rcond=None)[0][0])


def disagreement(rows):
    """{(subject, cell, R): {Rh: (d_q, d_FC, d_TF), 'NF': nf}} from count rows (Section 4)."""
    acc = collections.defaultdict(lambda: np.zeros(5))
    nf = collections.defaultdict(lambda: np.zeros(2))
    for r in rows:
        acc[(r["subject"], r["cell"], r["R"], r["Rh"])] += [r["dis"], r["n10"], r["t1"],
                                                           r["n"] * r["nb"], r["n1"] * r["nb"]]
        if r["Rh"] == RH[0]:                       # each display counted once for NF
            nf[(r["subject"], r["cell"], r["R"])] += [r["n1"], r["n"]]
    out = collections.defaultdict(dict)
    for (s, cell, R, Rh), (dis, n10, t1, nn, n1n) in acc.items():
        out[(s, cell, R)][Rh] = (dis / nn, n10 / n1n if n1n else float("nan"), t1 / nn)
    for key, (n1, n) in nf.items():
        out[key]["NF"] = n1 / n
    return out


def estimates(dis):
    """Per (subject, cell, R): both extrapolants' estimates of FC, q and TF."""
    x1 = np.array([h ** -0.5 for h in RH])
    x2 = np.array([1.0 / h for h in RH])
    est = {}
    for key, d in dis.items():
        e = {"NF": d["NF"]}
        for qi, qn in enumerate(("q", "FC", "TF")):
            y = np.array([d[h][qi] for h in RH])
            e[qn] = (intercept(x1, y), intercept(x2, y)) if np.isfinite(y).all() else (np.nan, np.nan)
        est[key] = e
    return est


def summarize(est, cell, R, qn, rng):
    """Median estimates, per-extrapolant intervals, and the conservative lower bound L."""
    subs = sorted({s for (s, c, r) in est if c == cell and r == R})
    e1 = np.array([est[(s, cell, R)][qn][0] for s in subs])
    e2 = np.array([est[(s, cell, R)][qn][1] for s in subs])
    idx = rng.integers(0, len(subs), size=(B, len(subs)))
    m1, m2 = np.nanmedian(e1[idx], axis=1), np.nanmedian(e2[idx], axis=1)
    return dict(med_sqrt=float(np.nanmedian(e1)), med_inv=float(np.nanmedian(e2)),
                ci_sqrt=[float(np.percentile(m1, 2.5)), float(np.percentile(m1, 97.5))],
                ci_inv=[float(np.percentile(m2, 2.5)), float(np.percentile(m2, 97.5))],
                L=float(np.percentile(np.minimum(m1, m2), 2.5)), n_subjects=len(subs))


def semi_synthetic_bias(sim):
    """Protocol Section 5: median over subjects of min(estimates) - exact, for FC(8) and q(8)."""
    exact = collections.defaultdict(lambda: np.zeros(3))
    for r in sim:
        if r["R"] == 8 and r["Rh"] == "exact":
            exact[r["subject"]] += [r["n10"], r["n1"], r["dis"] / r["n"]]
    ref = [dict(r, cell="primary") for r in sim if r["R"] == 8 and r["Rh"] != "exact"]
    est = estimates(disagreement(ref))
    out = {}
    for qn in ("FC", "q"):
        diffs = []
        for s, (n10, n1, qsum) in exact.items():
            ndraws = sum(1 for r in sim if r["subject"] == s and r["R"] == 8 and r["Rh"] == "exact")
            ex = n10 / n1 if qn == "FC" else qsum / ndraws
            diffs.append(min(est[(s, "primary", 8)][qn]) - ex)
        out[qn] = dict(median_bias=float(np.median(diffs)), per_subject=[float(x) for x in diffs])
    return out


def main():
    rows = [json.loads(l) for l in open(f"{OUT}/rows.jsonl")]
    sim = [json.loads(l) for l in open(f"{OUT}/semi_synthetic.jsonl")]
    lines = []

    def out(s=""):
        lines.append(s)
        print(s)

    rng = np.random.default_rng(BOOT_SEED)
    est = estimates(disagreement(rows))
    order = [n for n, _ in CELLS]
    cells = sorted({c for (_, c, _) in est}, key=order.index)
    bias = semi_synthetic_bias(sim)
    out("CONFIRMATORY ANALYSIS (docs/note_protocol.md): medians over subjects; extrapolants linear in")
    out("R_h^-1/2 and in 1/R_h; L = conservative lower bound (Section 4)")
    out(f"\nEstimator check (Section 5), R=8, primary cell: median bias of the conservative estimator "
        f"FC {bias['FC']['median_bias']:+.3f}, q {bias['q']['median_bias']:+.3f}")
    guard = max(0.0, bias["FC"]["median_bias"]) if bias["FC"]["median_bias"] > 0.05 else 0.0
    if guard:
        out(f"  bias guard active: L for FC(8) is reduced by {guard:.3f} (Section 6)")

    S = {}
    for cell in cells:
        for R in RU:
            for qn in QUANTITIES:
                S[(cell, R, qn)] = summarize(est, cell, R, qn, rng)
            nf = np.median([est[(s, cell, R)]["NF"] for (s, c, r) in est if c == cell and r == R])
            S[(cell, R, "NF")] = float(nf)
    out("\ncell            R   NF     TF (sqrt/inv)   q (sqrt/inv)    FC (sqrt/inv)   FC L")
    for cell in cells:
        for R in RU:
            f, q, t = S[(cell, R, "FC")], S[(cell, R, "q")], S[(cell, R, "TF")]
            out(f"{cell:14s} {R:2d}  {S[(cell, R, 'NF')]:.3f}  {t['med_sqrt']:.3f}/{t['med_inv']:.3f}"
                f"     {q['med_sqrt']:.3f}/{q['med_inv']:.3f}     {f['med_sqrt']:.3f}/{f['med_inv']:.3f}"
                f"    {f['L'] - (guard if R == 8 else 0.0):.3f}")

    L = lambda cell, R, qn: S[(cell, R, qn)]["L"] - (guard if (qn == "FC" and R == 8) else 0.0)
    c1 = L("primary", 8, "FC") > 0.5
    c1r = L("rsa_cv_mds", 8, "FC") > 0.5
    c2 = L("primary", 32, "FC") > 0.25
    c3 = all(L("primary", R, "q") > 0.10 for R in RU)
    holds = [c for c in cells if L(c, 8, "FC") > 0.5]
    general = len(holds) >= len(cells) - 2
    out(f"\nC1  FC(8) > 0.5, primary:          L = {L('primary', 8, 'FC'):.3f} -> {'SUPPORTED' if c1 else 'not supported'}")
    out(f"C1-RSA FC(8) > 0.5, rsa_cv_mds:    L = {L('rsa_cv_mds', 8, 'FC'):.3f} -> {'SUPPORTED' if c1r else 'not supported'}")
    out(f"C2  FC(32) > 0.25, primary:        L = {L('primary', 32, 'FC'):.3f} -> {'SUPPORTED' if c2 else 'not supported'}")
    out(f"C3  q(R) > 0.10 for R = 8, 16, 32: L = " + ", ".join(f"{L('primary', R, 'q'):.3f}" for R in RU)
        + f" -> {'SUPPORTED' if c3 else 'not supported'}")
    out(f"Generality: C1's criterion holds in {len(holds)} of {len(cells)} cells -> "
        + ("stated as general" if general else "restricted to: " + ", ".join(holds)))
    if c1 and c1r:
        verdict = "headline holds for standard and RSA-style maps"
    elif c1:
        verdict = "headline restricted to plug-in maps; the RSA result is reported as its limit"
    else:
        verdict = "no headline: report the estimates and recommend stopping"
    out(f"\nWhat the note may say (Section 6): {verdict}")
    open(f"{OUT}/confirmatory_analysis.txt", "w").write("\n".join(lines) + "\n")
    json.dump(dict(bias=bias, guard=guard, C1=c1, C1_RSA=c1r, C2=c2, C3=c3, cells_holding=holds,
                   general=general, verdict=verdict,
                   summary={f"{c}|{R}|{qn}": v for (c, R, qn), v in S.items()}),
              open(f"{OUT}/confirmatory_summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
