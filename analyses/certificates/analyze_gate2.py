#!/usr/bin/env python3
"""Gate 2 analysis: applies the pre-registered criteria of plan Section 6 to the saved runs.

Reads results/gate2/validity.jsonl and results/gate2/power.jsonl; writes
results/gate2/gate2_analysis.txt and results/gate2/gate2_summary.json.
"""
import json, collections, os, sys
import numpy as np
from scipy.stats import beta, binom

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # the repository root
from analyses.certificates import run_gate2_validity as V
from analyses.certificates import run_gate2_power as P

ENGINES = ("ball", "pair_t", "pair_pool", "cmp")
EXACT_UNDER = {"ball": ("gauss", "gauss_het"), "pair_t": ("gauss", "gauss_het"),
               "pair_pool": ("gauss",)}
ROBUST_SET = ("gauss_het", "t5", "rep_t5", "aniso4")
lines = []


def out(s=""):
    lines.append(s)
    print(s)


def cp_upper(x, n, level=0.95):
    return 1.0 if x == n else float(beta.ppf(level, x + 1, n - x))


def failed(e, r):
    if e == "cmp":
        return r["elem_viol"] + r["reach_viol"] + r.get("elem_viol_all_triples", 0) > 0
    return r["pair_viol"] + r["reach_viol"] + r.get("ball_viol", 0) > 0


# ------------------------------------------------------------------ validity
C = V.configs()
rows = collections.defaultdict(list)
for line in open(V.OUT):
    o = json.loads(line)
    rows[o["cfg"]].append(o)

out("GATE 2 VALIDITY: family-wise failure rate per engine (draws with any false statement)")
out("cfg geometry   noise      R  gamma draws | " + " | ".join(f"{e:>16s}" for e in ENGINES))
fwer = {}
for c in C:
    rs = rows[c["index"]]
    n = len(rs)
    cells = []
    for e in ENGINES:
        x = sum(failed(e, r[e]) for r in rs)
        fwer[(c["index"], e)] = (x, n)
        cells.append(f"{x:4d}/{n:<4d} ub {cp_upper(x, n):.3f}" if n else "      -        ")
    out(f"{c['index']:3d} {c['geometry']:10s} {c['model']:9s} {c['R']:2d}  {c['gamma']:.2f} {c['draws']:5d} | "
        + " | ".join(cells))
out("  (x/n = failing draws / draws; ub = one-sided 95% Clopper-Pearson upper bound)")

out("\nReach claims only (points whose true reach left its interval), failing draws:")
for c in C:
    rs = rows[c["index"]]
    out(f"{c['index']:3d} {c['geometry']:10s} {c['model']:9s} R={c['R']:2d} gamma={c['gamma']:.2f}: "
        + "  ".join(f"{e}={sum(r[e]['reach_viol'] > 0 for r in rs)}" for e in ENGINES))

out("\nG2-V (exactness under the engine's own model)")
out("  original clauses: rate <= 0.05 and ub <= 0.075; amended (plan Amendment 1): exact one-sided")
out("  binomial test of 'rate <= 0.05' not rejected at level 0.01, i.e. P(X >= x | p = 0.05) >= 0.01")
g2v_orig = g2v = True
for e, models in EXACT_UNDER.items():
    for c in C:
        if c["model"] in models:
            x, n = fwer[(c["index"], e)]
            orig = n > 0 and x / n <= 0.05 and cp_upper(x, n) <= 0.075
            pval = float(binom.sf(x - 1, n, 0.05)) if x > 0 else 1.0
            amend = n > 0 and pval >= 0.01
            g2v_orig &= orig
            g2v &= amend
            out(f"  {e:10s} cfg {c['index']:2d} {c['geometry']:9s} {c['model']:9s} R={c['R']:2d} g={c['gamma']:.2f}: "
                f"{x:3d}/{n:<4d} rate {x / n:.3f} ub {cp_upper(x, n):.3f} p {pval:.3f}  original {'ok' if orig else 'FAIL'}"
                f"  amended {'ok' if amend else 'FAIL'}")
out(f"  G2-V under the original clauses {'PASSES' if g2v_orig else 'FAILS'}; "
    f"under Amendment 1 {'PASSES' if g2v else 'FAILS'}")

out("\nG2-R (robustness set on digits, R=8, gamma=0.15: rate <= 0.10 on each)")
robust = {}
for e in ENGINES:
    ok, cells = True, []
    for model in ROBUST_SET:
        c = [c for c in C if c["geometry"] == "digits" and c["model"] == model and c["R"] == 8
             and c["gamma"] == 0.15][0]
        x, n = fwer[(c["index"], e)]
        ok &= n > 0 and x / n <= 0.10
        cells.append(f"{model}={x}/{n}")
    robust[e] = ok
    out(f"  {e:10s} {'robust' if ok else 'not robust'}: " + "  ".join(cells))

# ------------------------------------------------------------------ power
if not os.path.exists(P.OUT) or "--validity-only" in sys.argv:
    open("results/gate2/validity_analysis.txt", "w").write("\n".join(lines) + "\n")
    sys.exit(0)
PC = P.configs()
prow = collections.defaultdict(list)
for line in open(P.OUT):
    o = json.loads(line)
    prow[o["cfg"]].append(o)


def pcfg(**kw):
    return [c for c in PC if all(c.get(k) == v for k, v in kw.items())]


def mean_over(recs, f):
    v = [f(r) for r in recs]
    return float(np.mean(v)) if v else float("nan")


def resolved(r, e, eps="0.2", kp="40"):
    L = r["engines"][e]["labels"][eps][kp]
    return L["faithful"] + L["distorted"]


def tsne(recs):
    return [r for r in recs if r["display"].startswith("tsne")]


out("\nGATE 2 POWER: digits, t-SNE displays (3 seeds x 5 draws), eps=0.2")
out("primary metric: resolved share (faithful + distorted) at k'=40")
prim = {}
for model in ("gauss", "gauss_het"):
    for gamma in (0.15, 0.30):
        for R in (8, 16, 32):
            c = pcfg(geometry="digits", model=model, gamma=gamma, R=R)[0]
            rs = tsne(prow[c["index"]])
            vals = {e: mean_over(rs, lambda r: resolved(r, e)) for e in ENGINES}
            prim[(model, gamma, R)] = vals
            out(f"  {model:9s} gamma={gamma:.2f} R={R:2d}: " + "  ".join(f"{e}={v:.3f}" for e, v in vals.items())
                + f"   truth faithful@40={mean_over(rs, lambda r: r['truth']['0.2']['40']):.3f}")

cands = [e for e in ENGINES if robust[e] and e != "cmp"] or []
if cands:
    default = max(cands, key=lambda e: prim[("gauss", 0.15, 16)][e])
    scope = "robust"
else:
    default = max(("ball", "pair_t"), key=lambda e: prim[("gauss", 0.15, 16)][e])
    scope = "Gaussian only (no engine robust): residual diagnostic mandatory before Gate 3"
out(f"\nDefault engine (G2-R): {default}  [{scope}]")
p8, p16 = prim[("gauss", 0.15, 8)][default], prim[("gauss", 0.15, 16)][default]
g2p = p8 >= 0.20 and p16 >= 0.30
out(f"G2-P: {default} primary metric {p8:.3f} at R=8 (floor 0.20), {p16:.3f} at R=16 (floor 0.30)"
    f" -> {'PASSES' if g2p else 'FAILS'}")
out(f"GATE 2 {'PASSES' if (g2v and g2p) else 'FAILS'} (G2-V and G2-P)")

out("\nLabels by engine, digits t-SNE, gauss, gamma=0.15, eps=0.2: faithful / distorted / unresolved")
for R in (8, 16, 32):
    c = pcfg(geometry="digits", model="gauss", gamma=0.15, R=R)[0]
    rs = tsne(prow[c["index"]])
    for kp in ("20", "30", "40", "60", "80"):
        cells = []
        for e in ENGINES:
            f = mean_over(rs, lambda r: r["engines"][e]["labels"]["0.2"][kp]["faithful"])
            d = mean_over(rs, lambda r: r["engines"][e]["labels"]["0.2"][kp]["distorted"])
            cells.append(f"{e}={f:.3f}/{d:.3f}/{1 - f - d:.3f}")
        out(f"  R={R:2d} k'={kp:>2s} truth={mean_over(rs, lambda r: r['truth']['0.2'][kp]):.3f}  " + "  ".join(cells))
    out(f"  R={R:2d} ball witness ceiling (faithful, k'=40) = {mean_over(rs, lambda r: r['ball_ceiling_k40_eps0.2']):.3f}"
        f"; median log2(r_hi/r_lo): " + "  ".join(
            f"{e}={mean_over(rs, lambda r: r['engines'][e]['width']['0.2']):.2f}" for e in ENGINES))

out("\nFalse labels over every power record (must be 0 on the event):")
tot = collections.Counter()
for recs in prow.values():
    for r in recs:
        for e in ENGINES:
            for eps, d in r["engines"][e]["labels"].items():
                for kp, L in d.items():
                    tot[e] += L["false_faithful"] + L["false_distorted"]
out("  " + "  ".join(f"{e}={tot[e]}" for e in ENGINES))

out("\nG2-N (need): naive map mislabel rate at k'=40, eps=0.2, digits t-SNE (plug-in / bias-corrected)")
need = {}
for model in ("gauss", "gauss_het"):
    for gamma in (0.15, 0.30):
        for R in (8, 16, 32):
            c = pcfg(geometry="digits", model=model, gamma=gamma, R=R)[0]
            rs = tsne(prow[c["index"]])
            pl = mean_over(rs, lambda r: r["naive"]["plugin"]["0.2"]["40"]["mislabel"])
            co = mean_over(rs, lambda r: r["naive"]["corrected"]["0.2"]["40"]["mislabel"])
            ff = mean_over(rs, lambda r: r["naive"]["plugin"]["0.2"]["40"]["false_faithful"])
            need[(model, gamma, R)] = pl
            out(f"  {model:9s} gamma={gamma:.2f} R={R:2d}: plug-in {pl:.3f}  corrected {co:.3f}"
                f"  (plug-in false-faithful points per display: {ff:.1f})")

out("\nPCA-2 displays, digits gauss gamma=0.15, eps=0.2, k'=40: truth faithful and engine resolved")
for R in (8, 16, 32):
    c = pcfg(geometry="digits", model="gauss", gamma=0.15, R=R)[0]
    rs = [r for r in prow[c["index"]] if r["display"] == "pca2"]
    out(f"  R={R:2d} truth={mean_over(rs, lambda r: r['truth']['0.2']['40']):.3f}  "
        + "  ".join(f"{e}={mean_over(rs, lambda r: resolved(r, e)):.3f}" for e in ENGINES))

out("\nSynthetic families, PCA-2, gauss, R=8, eps=0.2: faithful share at k'=20 and k'=40")
for gamma in (0.15, 0.30):
    for fam in ("folded_curve", "thin_bridge", "imbalanced_mixture", "curved_surface"):
        c = pcfg(geometry=fam, gamma=gamma)[0]
        rs = prow[c["index"]]
        cells = []
        for kp in ("20", "40"):
            cells.append(f"k'={kp}: truth={mean_over(rs, lambda r: r['truth']['0.2'][kp]):.2f} " + " ".join(
                f"{e}={mean_over(rs, lambda r: r['engines'][e]['labels']['0.2'][kp]['faithful']):.2f}" for e in ENGINES))
        out(f"  gamma={gamma:.2f} {fam:18s} " + "   ".join(cells))

out("\nPlanning rule (plan 5.8): median of per-point c (quartiles), eps=0.2, points with reach <= 4k")
for label, sel in (("digits t-SNE gauss g=0.15", lambda R: tsne(prow[pcfg(geometry='digits', model='gauss', gamma=0.15, R=R)[0]['index']])),
                   ("digits t-SNE gauss g=0.30", lambda R: tsne(prow[pcfg(geometry='digits', model='gauss', gamma=0.30, R=R)[0]['index']]))):
    for R in (8, 16, 32):
        rs = sel(R)
        cells = []
        for e in ("ball", "pair_t", "pair_pool"):
            q = np.median([r["engines"][e]["planning_c_quantiles"] for r in rs], axis=0)
            cells.append(f"{e}={q[2]:.2f} [{q[1]:.2f},{q[3]:.2f}]")
        out(f"  {label} R={R:2d}: " + "  ".join(cells))
for gamma in (0.15, 0.30):
    cells = []
    for fam in ("folded_curve", "thin_bridge", "imbalanced_mixture", "curved_surface"):
        rs = prow[pcfg(geometry=fam, gamma=gamma)[0]["index"]]
        q = np.median([r["engines"]["pair_t"]["planning_c_quantiles"] for r in rs], axis=0)
        cells.append(f"{fam}={q[2]:.2f} [{q[1]:.2f},{q[3]:.2f}] (d={rs[0]['idim']:.1f})")
    out(f"  synthetic PCA-2 pair_t gamma={gamma:.2f}: " + "  ".join(cells))

out("\nWall time per call, digits R=8 (median seconds): " + "  ".join(
    f"{k}={np.median([r['seconds'][k] for r in prow[pcfg(geometry='digits', model='gauss', gamma=0.15, R=8)[0]['index']]]):.3f}"
    for k in ("ball", "pair_t", "pair_pool", "cmp")))

open("results/gate2/gate2_analysis.txt", "w").write("\n".join(lines) + "\n")
json.dump(dict(fwer={f"{k[0]}:{k[1]}": v for k, v in fwer.items()}, robust=robust, default=default,
               scope=scope, g2v=g2v, g2v_original=g2v_orig, g2p=g2p, primary={f"{k[0]}:{k[1]}:{k[2]}": v for k, v in prim.items()},
               need={f"{k[0]}:{k[1]}:{k[2]}": v for k, v in need.items()}),
          open("results/gate2/gate2_summary.json", "w"), indent=1)
