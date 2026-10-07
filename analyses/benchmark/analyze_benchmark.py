#!/usr/bin/env python3
"""The common checks rebuilt around acceptance and yield (docs/revision_protocol.md, Section 6,
items 5 to 8), from the per-point rerun (run_points.py).

For each check and budget, counts are pooled within each subject over splits and seeds. The
estimator and the frozen thresholds come from results/revision/dev/decisions.json. Every result
below is descriptive.
- **Native operating point:** acceptance; the reference disagreement at each size; both
  extrapolants and the chosen estimator; the status; the finite-reference disagreement at the
  deepest reference; the random comparator matched per map; enrichment; yield.
- **Matched acceptance** at 5, 10, 20, 30 and 50% of points per map, ties included, for the checks
  with a continuous score, and for the naive label of the cv-MDS maps.
- **Continuous overlap** with the deepest reference, for accepted and for all neighborhoods.
- **The alternative scDEED null.**
- **The cv-MDS diagnostics.**

Writes results/benchmark/benchmark.txt and summary.json.
"""
import collections, json, math, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from analyses.validation import analyze_validation as AV
from reach import checks, sim

OUT = "results/benchmark"
LEVELS = (0.05, 0.10, 0.20, 0.30, 0.50)
DATA = {"EEG": ("things", sim.LAYOUT["EEG"], 8 + 40), "mouse": ("allen", sim.LAYOUT["mouse"], 10 + 24)}
KP = 20
# accept first: the smallest (+1) or the largest (-1) score
SCORED = {"N": ("naive", 1), "TW": ("tw", -1), "SCD": ("scd", -1), "CV": ("cv", 1), "BOOT": ("boot", -1)}


def load(ds_dir, s):
    z = np.load(f"results/points/{ds_dir}/sub-{s:02d}.npz")
    return {k: z[k] for k in z.files}


def native(Z, X, prefix="tsne"):
    """The check's own labels per map (maps x points)."""
    g = lambda k: Z[f"{prefix}_{k}"]
    N = g("naive") <= KP
    if X == "N":
        return N
    if X == "TW":
        return np.array([checks.top_matched(tw, int(n.sum())) for tw, n in zip(g("tw"), N)])
    if X == "SCD":
        return g("scd") >= g("scd_cut")[:, None]
    if X == "SCD_alt":
        return g("scd") >= g("scd_cut_alt")[:, None]
    if X == "CV":
        return g("cv") <= KP
    if X == "SEED":
        return g("seed_count") == 2
    if X == "BOOT":
        return g("boot") >= checks.BOOT_SUPPORT - 1e-9
    if X == "SH":
        return g("shA") & g("shB")
    raise KeyError(X)


def matched(score, sign, level):
    """Per map, the share `level` of points with the best score, ties at the cut-off included."""
    n = score.shape[1]
    k = max(1, math.ceil(level * n))
    s = sign * score.astype(float)
    cut = np.sort(s, axis=1)[:, k - 1]
    return s <= cut[:, None]


def pooled(A, Z, sel, prefix, L):
    """Counts of the accepted neighborhoods A (maps x points) against every reference, in sim.FIELDS
    order, with the random comparator matched per map, pooled over the selected maps."""
    idx = Z["ref_index"]
    reach = Z[f"{prefix}_ref_reach"][sel]                 # maps x sub-blocks x points
    A = A[sel]
    c = {}
    for fam, name in ((0, "plug"), (1, "cv")):
        for Rh in L["RH"] + (L["extra"],):
            sub = np.flatnonzero((idx[:, 0] == fam) & (idx[:, 1] == Rh))
            T = reach[:, sub, :] <= KP
            nb, n_maps, n = len(sub), A.shape[0], A.shape[1]
            n1 = A.sum()
            c[(name, Rh)] = np.array([n1 * nb, (A[:, None, :] & ~T).sum(), (A[:, None, :] != T).sum(), T.sum(),
                                      n * nb * n_maps, n1, n * n_maps], float)
            if name == "plug":
                share_unfaithful = 1.0 - T.mean(axis=(1, 2))           # per map, over its sub-blocks
                c[("rand", Rh)] = np.array([n1 * nb, (A.sum(axis=1) * nb * share_unfaithful).sum(), 0, 0, 0, 0, 0])
    return c


def subject_record(A, Z, sel, prefix, L, overlap_ref):
    c = pooled(A, Z, sel, prefix, L)
    rec = AV.record(c, L)
    deep = L["RH"][-1]
    rec["random"] = AV.intercept(L["RH"], [AV.fc(c, "rand", h) for h in L["RH"]], "sqrt")
    rec["d_deep"] = AV.fc(c, "plug", deep)
    rec["rand_deep"] = AV.fc(c, "rand", deep)
    rec["accept"] = c[("plug", deep)][5] / c[("plug", deep)][6]
    ov = Z[f"{prefix}_ref_overlap"][sel][:, overlap_ref, :].astype(float) / 10.0     # maps x sub-blocks x points
    Am = A[sel]
    rec["overlap_accepted"] = float(ov[:, :, :][np.broadcast_to(Am[:, None, :], ov.shape)].mean()) if Am.any() else np.nan
    rec["overlap_all"] = float(ov.mean())
    return rec


def summarize(recs, est, th, fam, L):
    finite = [r for r in recs if np.isfinite(r.get(est, np.nan))]
    med = lambda k, rr=recs: float(np.nanmedian([r[k] for r in rr])) if rr else np.nan
    e = med(est)
    rnd = med("random")
    out = dict(status=AV.status(recs, th, fam) if len(finite) >= len(recs) / 2 else "insufficient acceptance",
               subjects_with_acceptance=int(sum(r["acc"] > 0 for r in recs)), accept=med("accept"),
               d=[float(np.nanmedian([r["curve"][f"plug|{h}"] for r in recs])) for h in L["RH"]],
               d_deep=med("d_deep"), rand_deep=med("rand_deep"), sqrt=med("sqrt"), inv=med("inv"),
               chosen=e, random=rnd, overlap_accepted=med("overlap_accepted"), overlap_all=med("overlap_all"))
    out["enrichment"] = float(np.nanmedian([(1 - r[est]) / (1 - r["random"]) for r in finite
                                            if 1 - r["random"] > 0.005])) if finite else np.nan
    out["enrichment_finite"] = float(np.nanmedian([(1 - r["d_deep"]) / (1 - r["rand_deep"]) for r in recs
                                                   if np.isfinite(r["d_deep"]) and 1 - r["rand_deep"] > 0.005]))
    out["yield"] = float(np.nanmedian([r["accept"] * (1 - r[est]) for r in finite])) if finite else np.nan
    return out


def main():
    dec = json.load(open("results/revision/dev/decisions.json"))
    est, th = dec["estimator"], dec["thresholds"]
    est = "smaller" if est == "published" else est
    fam = AV.FAMILY_OF.get(est, "plug")
    os.makedirs(OUT, exist_ok=True)
    lines, summary = [], dict(estimator=est, thresholds=th)

    def out(s=""):
        lines.append(s)
        print(s)

    out(f"BENCHMARK OF THE COMMON CHECKS (analyze_benchmark.py): estimator '{est}'; medians over subjects")
    for name, (ds_dir, L, fixed) in DATA.items():
        Zs = {s: load(ds_dir, s) for s in range(1, 11)}
        idx = Zs[1]["ref_index"]
        deep = np.flatnonzero((idx[:, 0] == 0) & (idx[:, 1] == L["RH"][-1]))
        summary[name] = {}
        out(f"\n=== {name}: F + H = {fixed} repetitions beside the analyst's R ===")
        out(f"  {'selection':16s} {'R':>3s} {'total':>5s} {'accept':>6s}  {'d at R_h':>20s}  {'chosen':>6s} {'random':>6s}"
            f" {'enrich':>6s} {'yield':>6s}  {'d_deep':>6s} {'enr_fin':>7s}  {'ovl acc/all':>11s}  status")
        for R in L["RU"]:
            sels = [(X, "tsne", X) for X in checks.NAMES + ("SCD_alt",)] + [("cv-MDS N", "mds", "N")]
            for X, sc in SCORED.items():
                sels += [(f"{X} @{int(l * 100)}%", "tsne", (X, l)) for l in LEVELS]
            sels += [(f"cv-MDS N @{int(l * 100)}%", "mds", ("N", l)) for l in LEVELS]
            for label, prefix, what in sels:
                recs = []
                for s, Z in Zs.items():
                    sel = Z[f"{prefix}_meta"][:, 1] == R
                    if isinstance(what, tuple):
                        X, level = what
                        key, sign = SCORED[X]
                        score = Z[f"{prefix}_naive"] if (prefix == "mds" or X == "N") else Z[f"{prefix}_{key}"]
                        A = matched(score, sign, level)
                    else:
                        A = native(Z, what, prefix)
                    recs.append(subject_record(A, Z, sel, prefix, L, deep))
                v = summarize(recs, est, th, fam, L)
                summary[name][f"{label}|{R}"] = v
                out(f"  {label:16s} {R:3d} {fixed + R:5d} {v['accept']:6.3f}  {' '.join(f'{x:.3f}' for x in v['d']):>20s}"
                    f"  {v['chosen']:6.3f} {v['random']:6.3f} {v['enrichment']:6.2f} {v['yield']:6.3f}"
                    f"  {v['d_deep']:6.3f} {v['enrichment_finite']:7.2f}  {v['overlap_accepted']:.2f}/{v['overlap_all']:.2f}"
                    f"  {v['status']}")
        # cv-MDS diagnostics
        diag = {}
        for R in L["RU"]:
            per = [(Z["mds_neg_share"][Z["mds_meta"][:, 1] == R].mean(), Z["mds_zero_share"][Z["mds_meta"][:, 1] == R].mean(),
                    np.median(Z["mds_stress"][Z["mds_meta"][:, 1] == R])) for Z in Zs.values()]
            diag[R] = dict(zip(("negative_share", "zero_share", "stress1"), np.median(np.array(per), axis=0).tolist()))
            out(f"  cv-MDS R={R}: negative cross-validated squared distances {diag[R]['negative_share']:.3f},"
                f" zeros in the MDS input {diag[R]['zero_share']:.3f}, Kruskal stress-1 {diag[R]['stress1']:.3f}"
                " (medians over subjects); naive labels rank the signed d^2, references use plug-in distances")
        summary[name]["cv_mds_diagnostics"] = diag
    open(f"{OUT}/benchmark.txt", "w").write("\n".join(lines) + "\n")
    json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
