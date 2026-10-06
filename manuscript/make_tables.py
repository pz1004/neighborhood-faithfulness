#!/usr/bin/env python3
"""Tables and in-text numbers for manuscript.tex, supplementary.tex and cover_letter.tex, drawn only
from committed results. Nothing is re-measured here.

The script writes into those files between "% BEGIN make_tables.py <name>" and
"% END make_tables.py <name>" lines, and nowhere else:

  tab_design              Table 1: datasets and design
  tab_validation          Table 2: the estimator in simulation, selected cells of Layer B
  tab_status              Table 3: each pre-specified claim's status and reporting class
  tab_primary             Table 4: the plug-in check, reported by class
  tab_benchmark           Table 5: common checks against the truth and a random selection
  tab_snr                 Table 6: exploratory signal-to-noise
  tab_supp_*              Supplementary Tables, numbered in the order the text cites them (the
                          list is in supplementary.tex). S2, the pre-registration trail, is written
                          in supplementary.tex.
  numbers                 every number the text quotes, as macros, so text and tables cannot drift;
                          each file's block defines the macros that file uses

With --check nothing is written: the script shows the difference and exits with an error if any
block differs from what the committed results give.

Sources:
- the pre-specified analyses: results/things_confirmatory/{rows.jsonl, confirmatory_summary.json};
  results/allen_confirmatory/{rows.jsonl, allen_summary.json, explore_snr.txt,
  subjects_described.json}; results/checks/checks_summary.json;
- the pilot's noise diagnostics: results/things_pilot/diagnostics.json;
- the analyses added after the first results (plan Amendment 6): results/revision/{validation/summary.json,
  dev/alpha.json}; results/revision/layer_c/{calibration_*.json, rule.json, oracle_*.jsonl.gz};
  results/gate_a/classification.json; results/benchmark/summary.json;
  results/session_disjoint/summary.json; results/manifest/mouse_overlap.json;
  results/review_audit/audit_summary.json.
"""
import difflib, gzip, json, re, sys
from pathlib import Path

import numpy as np

if sys.argv[1:] not in ([], ["--check"]):
    sys.exit("usage: make_tables.py [--check]")
CHECK_ONLY = sys.argv[1:] == ["--check"]
ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import analyze_things_confirmatory as AT      # noqa: E402  (the registered estimators)
import analyze_allen_confirmatory as AA       # noqa: E402

TH = json.load(open(ROOT / "results/things_confirmatory/confirmatory_summary.json"))
AL = json.load(open(ROOT / "results/allen_confirmatory/allen_summary.json"))
MICE = json.load(open(ROOT / "results/allen_confirmatory/subjects_described.json"))
CELL_ORDER = ["primary", "rsa_cv_mds", "k5", "k20", "eps0.1", "eps0.3", "tsne_perp15", "tsne_perp50",
              "pca2", "mds", "early_window", "unwhitened"]
CELL_LABEL = {"primary": "primary (t-SNE, perplexity 30)", "rsa_cv_mds": "cv-MDS", "k5": "$k=5$",
              "k20": "$k=20$", "eps0.1": "$\\varepsilon=0.1$", "eps0.3": "$\\varepsilon=0.3$",
              "tsne_perp15": "t-SNE perplexity 15", "tsne_perp50": "t-SNE perplexity 50",
              "pca2": "PCA-2", "mds": "MDS, Euclidean", "early_window": "alternative window",
              "unwhitened": "unwhitened"}
CHECK_LABEL = {"N": "Plug-in", "TW": "Trustworthiness", "SCD": "scDEED", "CV": "Cross-validated",
               "SEED": "Seed recurrence", "BOOT": "Bootstrap recurrence", "SH": "Split-half recurrence"}
CHECKS = ["N", "TW", "SCD", "CV", "SEED", "BOOT", "SH"]


def get(S, cell, R, q):
    return S["summary"][f"{cell}|{R}|{q}"]


def pct(x, d=0):
    return f"{100 * x:.{d}f}"


def span(lo, hi, d=0):
    a, b = pct(lo, d), pct(hi, d)
    return a if a == b else f"{a}--{b}"


def rng(S, cell, R, q, d=0):
    """Range spanned by the two extrapolants' medians, in percent."""
    v = get(S, cell, R, q)
    lo, hi = sorted((v["med_sqrt"], v["med_inv"]))
    return span(lo, hi, d)


def L(S, cell, R, q, d=0):
    return pct(get(S, cell, R, q)["L"], d)


def num(x, d=3, plus=True):
    """A signed number with a typeset minus; a value that rounds to zero is shown unsigned."""
    v = f"{x:+.{d}f}" if plus else f"{x:.{d}f}"
    if float(v) == 0:
        return f"{0:.{d}f}"
    return f"$-${v[1:]}" if v.startswith("-") else v


TABLES = {}


def write(name, text):
    """Keep a table's body; it goes into its block at the end."""
    TABLES[name] = text


# ---------------------------------------------------------------------------------------------
# Table 1: design
# ---------------------------------------------------------------------------------------------
cells = [m["n_cells"] for m in MICE]
n_f = sum(m["sex"] == "F" for m in MICE)
write("tab_design", "\\begin{tabular}{>{\\raggedright\\arraybackslash}p{0.26\\linewidth}>{\\raggedright\\arraybackslash}p{0.32\\linewidth}>{\\raggedright\\arraybackslash}p{0.32\\linewidth}}\n\\toprule\n"
      " & Human EEG (THINGS-EEG2) & Mouse V1 (Allen Brain Observatory) \\\\\n\\midrule\n"
      "Subjects & 10 adults (8 female, 2 male) & "
      f"10 mice ({10 - n_f} male, {n_f} female), one experiment each \\\\\n"
      "Stimuli & 200 object images & 118 natural scenes \\\\\n"
      "Repetitions per image & 80, over 4 sessions & 50, early and late strata of 25 \\\\\n"
      "Features per trial & 17 channels $\\times$ 80 samples (0--0.79\\,s) & "
      f"{min(cells)}--{max(cells)} cells, mean $\\Delta F/F$ over 0.5\\,s \\\\\n"
      "Blocks $F$ / $U$ / $H$ & 8 / 32 / 40 & 10 / 16 / 24 \\\\\n"
      "Analyst's budgets $R$ & 8, 16, 32 & 8, 16 \\\\\n"
      "Repetitions an assessment uses, $F+R+H$ & $48+R$ of 80 & $34+R$ of 50 \\\\\n"
      "Reference sizes $R_h$ (sub-blocks) & 8 (5), 20 (2), 40 (1) & 8 (3), 12 (2), 24 (1) \\\\\n"
      "Random splits per subject & 10 & 10 \\\\\n"
      "Representation, fitted on $F$ & \\multicolumn{2}{l}{PCA to 100, Ledoit--Wolf whitening, PCA to $D=20$} \\\\\n"
      "\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Table 6: exploratory signal-to-noise
# ---------------------------------------------------------------------------------------------
snr = (ROOT / "results/allen_confirmatory/explore_snr.txt").read_text()
blocks = dict(zip(("allen", "things"), snr.split("THINGS-EEG2")))
grab = lambda block, pat: float(re.search(pat, block).group(1))
S2 = {key: dict(pos=grab(b, r"positive cross-validated d\^2: ([\d.]+)"),
                rel=grab(b, r"distance matrix: ([\d.]+)"),
                d8=grab(b, r"8th-neighbor distance ([\d.]+)"),
                noise=grab(b, r"16-rep mean pair ([\d.]+)"),
                ratio=grab(b, r"ratio ([\d.]+)"),
                n_pos8=grab(b, r"median (\d+) of")) for key, b in blocks.items()}
write("tab_snr", "\\begin{tabular}{lrr}\n\\toprule\n & Human EEG & Mouse V1 \\\\\n\\midrule\n"
      f"Image pairs with a positive cross-validated distance & {S2['things']['pos']:.3f} & {S2['allen']['pos']:.3f} \\\\\n"
      f"Split-half reliability of the distance matrix & {S2['things']['rel']:.3f} & {S2['allen']['rel']:.3f} \\\\\n"
      f"Distance to the 8th nearest image, $d_{{(8)}}$ & {S2['things']['d8']:.2f} & {S2['allen']['d8']:.2f} \\\\\n"
      f"Noise between two 16-trial means, $\\sigma_{{16}}$ & {S2['things']['noise']:.2f} & {S2['allen']['noise']:.2f} \\\\\n"
      f"Ratio $d_{{(8)}}/\\sigma_{{16}}$ & {S2['things']['ratio']:.2f} & {S2['allen']['ratio']:.2f} \\\\\n"
      "\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Gate A: statuses and reporting classes (results/gate_a/, results/revision/layer_c/)
# ---------------------------------------------------------------------------------------------
GA = json.load(open(ROOT / "results/gate_a/classification.json"))
LCDIR = ROOT / "results/revision/layer_c"
LC = {p.stem[len("calibration_"):]: json.load(open(p)) for p in sorted(LCDIR.glob("calibration_*.json"))}
RULE = json.load(open(LCDIR / "rule.json"))
STATUS = {"estimable": "estimable", "extrapolation unresolved": "unresolved", "reference saturated": "saturated",
          "insufficient acceptance": "too few passes"}
# Layer C ran for the results whose matching cells were accurate (docs/revision_results.md, Section 4):
# EEG C1 on false comfort, mouse V1 C3 on the mislabel rate at both of its budgets; both truths must pass.
LC_RUNS = {"EEG C1": (["EEG_T1_tsne_R8", "EEG_T3_tsne_R8"], "FC"),
           "mouse C3": (["mouse_T1_tsne_R8", "mouse_T3_tsne_R8", "mouse_T1_tsne_R16", "mouse_T3_tsne_R16"], "q")}
LC_UNIT = {("EEG", "primary|8|FC"): "EEG C1", ("mouse", "primary|8|q"): "mouse C3",
           ("mouse", "primary|16|q"): "mouse C3"}


def lc_rate(tag, q):
    o = LC[tag]
    return o["published"]["L_above"] if q == "FC" else o["mislabel"]["L_above"]


def allowed(r):
    return 0.025 + 2 * r["mcse"]


def met(r):
    return r["p"] <= allowed(r)


CALIBRATED = {u: all(met(lc_rate(t, q)) for t in tags) for u, (tags, q) in LC_RUNS.items()}
assert CALIBRATED == {"EEG C1": True, "mouse C3": False}


def klass(ds, key, reading="local"):
    """The reporting class of a real result (manuscript Section 2.3.4)."""
    g = GA[ds][key][reading]
    if g.startswith("estimable"):
        u = LC_UNIT.get((ds, key))
        return "calibrated" if (u and CALIBRATED[u]) else "model-dependent"
    return {"model-dependent": "model-dependent", "finite-reference only": "finite-reference",
            "model-dependent (no matching validation cell)": "model-dependent, unvalidated"}[g]


def lc_text(ds, key):
    u = LC_UNIT.get((ds, key))
    if not u:
        return "not run"
    tags, q = LC_RUNS[u]
    missed = sorted({t.split("_")[1] for t in tags if not met(lc_rate(t, q))})
    return "met under T1 and T3" if not missed else "missed under " + " and ".join(missed)


def finite(ds, key, d=0):
    """The disagreement with the deepest reference, median over subjects with its 95% interval."""
    f = GA[ds][key]["finite_reference"]
    return f"{pct(f['median'], d)} ({pct(f['ci'][0], d)}--{pct(f['ci'][1], d)})"


def matching(ds, key):
    g = GA[ds][key]
    if g["status"] != "estimable":
        return "--"
    if not g["matches"]:
        return "none"
    return f"{len(g['matches'])} ({num(g['worst'][0])} / {g['worst'][1]:.3f})"


# Table 3: the claims
CLAIM_ROWS = [("Human EEG", "EEG", "C1: false comfort, t-SNE", "primary|8|FC"),
              ("", "EEG", "C1-RSA: false comfort, cv-MDS", "rsa_cv_mds|8|FC"),
              ("", "EEG", "C2: false comfort, t-SNE", "primary|32|FC"),
              ("", "EEG", "C3: mislabel rate", "primary|8|q"),
              ("", "EEG", "", "primary|16|q"),
              ("", "EEG", "", "primary|32|q"),
              ("Mouse V1", "mouse", "C1: false comfort, t-SNE", "primary|8|FC"),
              ("", "mouse", "C1-RSA: false comfort, cv-MDS", "rsa_cv_mds|8|FC"),
              ("", "mouse", "C3: mislabel rate", "primary|8|q"),
              ("", "mouse", "", "primary|16|q")]
P = lambda w: f">{{\\raggedright\\arraybackslash}}p{{{w}}}"
body = []
for name, ds, claim, key in CLAIM_ROWS:
    if name:
        if body:
            body.append("\\addlinespace")
        body.append(f"\\multicolumn{{7}}{{l}}{{\\textit{{{name}}}}} \\\\")
    R = key.split("|")[1]
    body.append(f"\\quad {claim} & {R} & {STATUS[GA[ds][key]['status']]} & {matching(ds, key)} & "
                f"{lc_text(ds, key)} & {klass(ds, key)} & {klass(ds, key, 'strict')} \\\\")
write("tab_status", "\\setlength{\\tabcolsep}{4pt}\n"
      f"\\begin{{tabular}}{{{P('3.8cm')}rl{P('2.2cm')}{P('2.0cm')}{P('2.0cm')}{P('2.0cm')}}}\n\\toprule\n"
      " & & & Matching cells & Study-level & \\multicolumn{2}{c}{Class} \\\\\n"
      "\\cmidrule(lr){6-7}\n"
      "Claim & Trials & Status & (worst bias / RMSE) & calibration & reported & strict reading \\\\\n\\midrule\n"
      + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Table 4: the plug-in check, reported by class
# ---------------------------------------------------------------------------------------------
BM = json.load(open(ROOT / "results/benchmark/summary.json"))
rows = []
for name, ds, S, Rs in (("Human EEG", "EEG", TH, (8, 16, 32)), ("Mouse V1", "mouse", AL, (8, 16))):
    d = 1 if ds == "mouse" else 0
    rows.append(f"\\multicolumn{{9}}{{l}}{{\\textit{{{name}}}}} \\\\")
    for pipe, cell, bkey in (("t-SNE", "primary", "N"), ("cv-MDS", "rsa_cv_mds", "cv-MDS N")):
        for R in Rs:
            key = f"{cell}|{R}|FC"
            k = klass(ds, key)
            est = rng(S, cell, R, "FC") if k != "finite-reference" else "--"
            bound = L(S, cell, R, "FC") if k == "calibrated" else "--"
            rows.append(f"\\quad {pipe if R == Rs[0] else ''} & {R} & {pct(S['summary'][f'{cell}|{R}|NF'], 1)} & "
                        f"{STATUS[GA[ds][key]['status']]} & {k} & {est} & {bound} & {finite(ds, key, d)} & "
                        f"{pct(BM[ds][f'{bkey}|{R}']['rand_deep'], d)} \\\\")
    rows.append("\\addlinespace")
write("tab_primary", "\\setlength{\\tabcolsep}{4pt}\n"
      f"\\begin{{tabular}}{{lrrl{P('2.3cm')}rrll}}\n\\toprule\n"
      " & & Called & & & False & Bound & \\multicolumn{2}{c}{Disagreement, deepest reference} \\\\\n"
      "\\cmidrule(lr){8-9}\n"
      "Map & Trials & faithful & Status & Class & comfort & $L$ & check (95\\% interval) & random \\\\\n"
      "\\midrule\n" + "\n".join(rows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Table 2 and Supplementary Tables: the simulation study (Layers A, B and C)
# ---------------------------------------------------------------------------------------------
VAL = json.load(open(ROOT / "results/revision/validation/summary.json"))


def vcell(ds, key):
    return VAL[ds]["cells"][key]


def perf(c, est, k):
    return c["perf"].get(est, {}).get(k, np.nan)


def cell_label(key):
    """A Layer B cell's label, e.g. 'T1_l0.5' -> 'T1, lambda = 0.5'."""
    truth, rest = key.split("_", 1)
    if rest == "S":
        return f"stress, {truth}"
    if rest == "zero":
        return f"{truth}, no noise"
    if rest == "cvmds":
        return f"{truth}, cv-MDS"
    if rest.startswith("D"):
        return f"shared noise $\\rho={rest[1:]}$, {truth}"
    return f"{truth}, $\\lambda={rest[1:]}$"


def vrow(ds, key):
    c = vcell(ds, key)
    name, R = key.split("|")
    fmt = lambda e: f"{num(perf(c, e, 'median_bias'))} ({perf(c, e, 'rmse'):.3f})"
    return (f"{cell_label(name)} & {R} & {c['exact']:.3f} & {STATUS[c['status']]} & {fmt('smaller')} & "
            f"{fmt('larger')} & {pct(c['bracket'])}")


VAL_ROWS = {"EEG": ["T1_l1|8", "T2_l1|8", "T3_l1|8", "T1_l0.5|8", "T1_l1|16", "T3_l1|16", "T1_l1|32", "T3_l1|32",
                    "T1_l2|8", "T1_cvmds|8", "T1_S|8", "T1_D0.3|8"],
            "mouse": ["T1_l1|8", "T1_l1|16", "T3_l1|16", "T1_l0.25|8", "T2_S|8", "T1_D0.3|8"]}
body = []
for ds, name in (("EEG", "Human EEG"), ("mouse", "Mouse V1")):
    body.append(f"\\multicolumn{{7}}{{l}}{{\\textit{{{name}}}}} \\\\")
    body += [f"\\quad {vrow(ds, k)} \\\\" for k in VAL_ROWS[ds]]
    body.append("\\addlinespace")
write("tab_validation", "\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{lrrlllr}\n\\toprule\n"
      " & & Exact false & & \\multicolumn{2}{c}{Bias (RMSE)} & Range \\\\\n\\cmidrule(lr){5-6}\n"
      "Simulated cell & Trials & comfort & Status & smaller & $1/R_h$ & covers (\\%) \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# Supplementary: every Layer B cell
body = []
for ds, name in (("EEG", "Human EEG"), ("mouse", "Mouse V1")):
    body.append(f"\\multicolumn{{9}}{{l}}{{\\textit{{{name}}}}} \\\\")
    for key, c in VAL[ds]["cells"].items():
        q = c["perf_q"].get("smaller", {})
        body.append(f"\\quad {vrow(ds, key)} & {STATUS[c['status_q']]} & {num(q.get('median_bias', np.nan), 4)} \\\\")
    body.append("\\addlinespace")
write("tab_supp_validation", "\\setlength{\\tabcolsep}{2.5pt}\\renewcommand{\\arraystretch}{0.94}\n\\begin{tabular}{lrrlllrlr}\n\\toprule\n"
      " & & Exact false & & \\multicolumn{2}{c}{Bias (RMSE)} & Range & \\multicolumn{2}{c}{Mislabel rate} \\\\\n"
      "\\cmidrule(lr){5-6}\\cmidrule(lr){8-9}\n"
      "Simulated cell & Trials & comfort & Status & smaller & $1/R_h$ & covers (\\%) & status & bias \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# Supplementary: Layer A
body = []
for ds, name, grid in (("EEG", "Human EEG", (8, 20, 40, 80, 160)), ("mouse", "Mouse V1", (8, 12, 24, 48, 96))):
    body.append(f"\\multicolumn{{9}}{{l}}{{\\textit{{{name}}}: $R_h$ = " + ", ".join(map(str, grid)) + "} \\\\")
    for key, a in VAL[ds]["layer_a"].items():
        body.append(f"\\quad {cell_label(key)} & {a['exact']:.3f} & " + " & ".join(f"{a['curve'][str(h)]:.3f}" for h in grid)
                    + f" & {num(a['grid_sqrt'])} & {num(a['deep_sqrt'])} \\\\")
    body.append("\\addlinespace")
write("tab_supp_layer_a", "\\begin{tabular}{lrrrrrrrr}\n\\toprule\n"
      " & Exact false & \\multicolumn{5}{c}{Disagreement $d_{\\FC}(R_h)$} & \\multicolumn{2}{c}{Error of the $R_h^{-1/2}$ fit} \\\\\n"
      "\\cmidrule(lr){3-7}\\cmidrule(lr){8-9}\n"
      "Simulated cell & comfort & \\multicolumn{5}{c}{by reference size, as listed} & real sizes & all sizes \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")


def oracle_range(tag, q):
    """The range of the virtual subjects' true values (pooled over their oracle draws)."""
    acc = {}
    for line in gzip.open(LCDIR / f"oracle_{tag}.jsonl.gz", "rt"):
        r = json.loads(line)
        acc.setdefault(r["v"], np.zeros(4))
        acc[r["v"]] += [r["n10"], r["n1"], r["dis"], r["n"]]
    v = [a[0] / a[1] if q == "FC" else a[2] / a[3] for a in acc.values()]
    return min(v), max(v)


# Supplementary: Layer C
body = []
for label, (tags, q) in LC_RUNS.items():
    for tag in tags:
        o = LC[tag]
        ds, truth, _, R = tag.split("_")
        r = lc_rate(tag, q)
        lo, hi = oracle_range(tag, q)
        target = o["Theta"] if q == "FC" else o["mislabel"]["Theta"]
        rule = RULE[tag][q]
        declares = o["published"]["declares_above_half"]["p"] if q == "FC" else o["mislabel"]["declares_above_tenth"]["p"]
        ostat = f"{o['order_statistics']['L_above']['p']:.3f}" if q == "FC" else "--"
        body.append(f"{label if tag == tags[0] else ''} & {truth} & {R[1:]} & {target:.3f} & {lo:.2f}--{hi:.2f} & "
                    f"{r['p']:.3f} $\\pm$ {r['mcse']:.3f} & {allowed(r):.3f} & {'met' if met(r) else 'missed'} & "
                    f"{pct(declares, 1)} & {ostat} & {rule['L_exact_above']['p']:.3f} & {num(rule['median_bias'])} \\\\")
    body.append("\\addlinespace")
write("tab_supp_layer_c", "\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{llrrlllrrrrr}\n\\toprule\n"
      " & & & Target & & & & & Declares & & Rule, exact & Median \\\\\n"
      "Claim & Truth & Trials & $\\Theta$ & Range & $\\Pr(L>\\Theta)$ & Allowed & Criterion & (\\%) & $x_{(2)}$ & values & bias \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# Supplementary: every real result's status and class
GA_LABEL = dict(CELL_LABEL, **{f"check:{X}": CHECK_LABEL[X] for X in CHECKS[1:]})
CK = json.load(open(ROOT / "results/checks/checks_summary.json"))
for tab, ds, ckds in (("tab_supp_status_eeg", "EEG", "things"), ("tab_supp_status_mouse", "mouse", "allen")):
    body = []
    order = [c for c in CELL_ORDER] + [f"check:{X}" for X in CHECKS[1:]]
    for unit in order:
        keys = sorted((k for k in GA[ds] if k.split("|")[0] == unit), key=lambda k: (int(k.split("|")[1]), k.split("|")[2]))
        for j, key in enumerate(keys):
            g = GA[ds][key]
            cell, R, q = key.split("|")
            if q == "q":
                lab = "\\quad mislabel rate"
            else:
                lab = GA_LABEL[cell] if j == 0 else ""
            pre = "--"
            if unit.startswith("check:") and R == "8":
                pre = {"leaves most false comfort": "leaves most", "removes most false comfort": "removes most",
                       "undetermined": "undetermined", "not evaluable": "not evaluable"}[CK["classification"]["things" if ds == "EEG" else "allen"][unit[6:]]]
            fr = g["finite_reference"]
            fin = "--" if not np.isfinite(fr["median"]) else f"{pct(fr['median'], 1)} ({pct(fr['ci'][0], 1)}--{pct(fr['ci'][1], 1)})"
            body.append(f"{lab} & {R} & {STATUS[g['status']]} & {matching(ds, key)} & {klass(ds, key)} & {fin} & {pre} \\\\")
        body.append("\\addlinespace")
    write(tab, "\n".join(body[:-1]) + "\n")                  # rows of a longtable (supplementary.tex)

# ---------------------------------------------------------------------------------------------
# Supplementary Tables: per subject, plasmode, every variation
# ---------------------------------------------------------------------------------------------
def per_subject(path, mod):
    """Per-subject disagreement and estimates of the primary cell at 8 trials, with the registered estimator."""
    rr = [r for r in map(json.loads, open(path)) if r["cell"] == "primary" and r["R"] == 8]
    dis = mod.disagreement(rr)
    return dis, mod.estimates(dis)


dis_th, est_th = per_subject(ROOT / "results/things_confirmatory/rows.jsonl", AT)
dis_al, est_al = per_subject(ROOT / "results/allen_confirmatory/rows.jsonl", AA)
d_mouse = [dis_al[(s, "primary", 8)][24][1] for s in range(1, 11)]
body = []
for s in range(1, 11):
    a, b = est_th[(s, "primary", 8)], est_al[(s, "primary", 8)]
    mouse = MICE[s - 1]
    body.append(f"{s} & {pct(a['NF'], 1)} & {pct(a['FC'][0], 1)} & {pct(a['FC'][1], 1)} & {pct(dis_th[(s, 'primary', 8)][40][1], 1)} & "
                f"{pct(b['NF'], 1)} & {pct(d_mouse[s - 1], 1)} & {mouse['sex']} & {mouse['age_days']} & {mouse['n_cells']} \\\\")
write("tab_supp_subjects", "\\begin{tabular}{rrrrrrrrrr}\n\\toprule\n"
      " & \\multicolumn{4}{c}{Human EEG} & \\multicolumn{5}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-5}\\cmidrule(lr){6-10}\n"
      "Subject & NF & FC, $R_h^{-1/2}$ & FC, $1/R_h$ & $d_{40}$ & NF & $d_{24}$ & Sex & Age (days) & Cells \\\\\n\\midrule\n"
      + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

body = []
for s in range(10):
    body.append(f"{s + 1} & {TH['bias']['FC']['per_subject'][s]:+.3f} & {TH['bias']['q']['per_subject'][s]:+.3f} & "
                f"{AL['bias']['FC']['per_subject'][s]:+.3f} & {AL['bias']['q']['per_subject'][s]:+.3f} \\\\")
body.append("\\midrule")
body.append(f"Median & {TH['bias']['FC']['median_bias']:+.3f} & {TH['bias']['q']['median_bias']:+.3f} & "
            f"{AL['bias']['FC']['median_bias']:+.3f} & {AL['bias']['q']['median_bias']:+.3f} \\\\")
write("tab_supp_plasmode", "\\begin{tabular}{rrrrr}\n\\toprule\n"
      " & \\multicolumn{2}{c}{Human EEG} & \\multicolumn{2}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n"
      "Subject & FC & $q$ & FC & $q$ \\\\\n\\midrule\n" + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

for tab, ds, S, Rs, d in (("tab_supp_cells_eeg", "EEG", TH, (8, 16, 32), 1), ("tab_supp_cells_mouse", "mouse", AL, (8, 16), 1)):
    body = []
    for cell in CELL_ORDER:
        for R in Rs:
            body.append(f"{CELL_LABEL[cell] if R == Rs[0] else ''} & {R} & {pct(S['summary'][f'{cell}|{R}|NF'], 1)} & "
                        f"{rng(S, cell, R, 'TF', 1)} & {rng(S, cell, R, 'q', 1)} & {rng(S, cell, R, 'FC', d)} & "
                        f"{L(S, cell, R, 'FC', d)} & {STATUS[GA[ds][f'{cell}|{R}|FC']['status']]} \\\\")
        body.append("\\addlinespace")
    write(tab, "\\begin{tabular}{lrrrrrrl}\n\\toprule\n"
          " & & Called & Truly & Mislabel & False & Bound & \\\\\n"
          "Variation & Trials & faithful & faithful & rate & comfort & $L$ & Status \\\\\n\\midrule\n"
          + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# The common checks: pre-specified outputs (Supplementary Tables) and the benchmark (Table 5)
# ---------------------------------------------------------------------------------------------
DS = {"things": ("Human EEG", (8, 16, 32)), "allen": ("Mouse V1", (8, 16))}


def ck(ds, X, R, q):
    return CK["datasets"][ds]["summary"][f"{X}|{R}|{q}"]


def signed(x, d=0):
    """Percent with a typeset minus sign; a value that rounds to zero is shown unsigned."""
    v = pct(x, d)
    if float(v) == 0:
        return pct(0.0, d)
    return f"$-${v[1:]}" if v.startswith("-") else v


def evaluable(ds, X, R):
    return ck(ds, X, R, "FC")["n_defined"] >= 5


def ck_rng(ds, X, R, q, d=0):
    """Range of the two extrapolants' medians; 'n.e.' where the check is not evaluable."""
    v = ck(ds, X, R, q)
    if q in ("FC", "FC_NX") and not evaluable(ds, X, R):
        return "n.e."
    if not (np.isfinite(v["med_sqrt"]) and np.isfinite(v["med_inv"])):
        return "--"
    lo, hi = sorted((v["med_sqrt"], v["med_inv"]))
    a, b = signed(lo, d), signed(hi, d)
    if a == b:
        return a
    return f"{a} to {b}" if lo < 0 else f"{a}--{b}"


def ck_bound(ds, X, R, which, d=0):
    v = ck(ds, X, R, "FC")
    return signed(v[which], d) if (evaluable(ds, X, R) and np.isfinite(v[which])) else "n.e."


def ck_pass(ds, X, R, q="PASS", d=1):
    return pct(ck(ds, X, R, q)["med"], d)


for tab, ds in (("tab_supp_checks_eeg", "things"), ("tab_supp_checks_mouse", "allen")):
    name, Rs = DS[ds]
    body = []
    for X in CHECKS:
        for R in Rs:
            body.append(f"{X if R == Rs[0] else ''} & {R} & {ck_pass(ds, X, R)} & "
                        f"{ck_rng(ds, X, R, 'FC', 1)} & {ck_bound(ds, X, R, 'L', 1)} & {ck_bound(ds, X, R, 'U', 1)} & "
                        f"{'--' if X == 'N' else ck_rng(ds, X, R, 'CATCH', 1)} & {ck_rng(ds, X, R, 'YIELD', 1)} & "
                        f"{ck_pass(ds, X, R, 'PASS_NX')} & {'--' if X == 'N' else ck_rng(ds, X, R, 'FC_NX', 1)} \\\\")
        body.append("\\addlinespace")
    write(tab, "\\begin{tabular}{lrrrrrrrrr}\n\\toprule\n"
          " & & & False & & & Plug-in false & & \\multicolumn{2}{c}{With plug-in} \\\\\n"
          "\\cmidrule(lr){9-10}\n"
          "Check & Trials & Passes & comfort & $L$ & $U$ & comfort flagged & Yield & Passes & False comfort \\\\\n\\midrule\n"
          + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

body = []
for X in CHECKS:
    cells_ = []
    for ds in DS:
        b = CK["datasets"][ds]["bias"][X]
        cells_.append(f"{b['bias_lo']:+.3f} & {b['bias_hi']:+.3f} & {b['n']}"
                      if b["n"] else "-- & -- & 0")
    body.append(f"{CHECK_LABEL[X]} & " + " & ".join(cells_) + " \\\\")
write("tab_supp_checks_plasmode", "\\begin{tabular}{lrrrrrr}\n\\toprule\n"
      " & \\multicolumn{3}{c}{Human EEG} & \\multicolumn{3}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
      "Check & Smaller & Larger & Subjects & Smaller & Larger & Subjects \\\\\n\\midrule\n"
      + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

# Table 5: the benchmark at 8 trials
BENCH = [("N", "Plug-in"), ("TW", "Trustworthiness"), ("SCD", "scDEED"), ("SCD_alt", "scDEED, coordinate null"),
         ("CV", "Cross-validated"), ("SEED", "Seed recurrence"), ("BOOT", "Bootstrap recurrence"),
         ("SH", "Split-half recurrence"), ("cv-MDS N", "cv-MDS label")]
SCORED = {"N", "TW", "SCD", "CV", "BOOT", "cv-MDS N"}


def bm(ds, X, R=8):
    return BM[ds][f"{X}|{R}"]


body = []
for ds, name in (("EEG", "Human EEG"), ("mouse", "Mouse V1")):
    d = 1 if ds == "mouse" else 0
    body.append(f"\\multicolumn{{9}}{{l}}{{\\textit{{{name}}}}} \\\\")
    for X, lab in BENCH:
        e = bm(ds, X)
        ok = e["status"] == "estimable"
        five = BM[ds][f"{X} @5%|8"] if X in SCORED else None
        five_txt = "--" if five is None else (pct(five["chosen"]) if five["status"] == "estimable" else STATUS[five["status"]])
        has = np.isfinite(e["d_deep"]) and e["accept"] > 0
        latent = f"{pct(e['chosen'])} & {pct(e['random'])} & {e['enrichment']:.1f}" if ok else "-- & -- & --"
        deep = f"{pct(e['d_deep'], d)} & {pct(e['rand_deep'], d)}" if has else "-- & --"
        body.append(f"\\quad {lab} & {pct(e['accept'], 1)} & {STATUS[e['status']]} & {latent} & {five_txt} & {deep} \\\\")
    body.append("\\addlinespace")
write("tab_benchmark", "\\setlength{\\tabcolsep}{4pt}\n"
      f"\\begin{{tabular}}{{lr{P('1.6cm')}rrrrrr}}\n\\toprule\n"
      " & & & Estimated & Random & Enrich- & At 5\\% & \\multicolumn{2}{c}{Disagreement} \\\\\n"
      "\\cmidrule(lr){8-9}\n"
      "Check & Passes & Status & false comfort & selection & ment & acceptance & check & random \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# Supplementary: matched acceptance
LEVELS = (5, 10, 20, 30, 50)
body = []
for ds, name in (("EEG", "Human EEG: estimated false comfort"), ("mouse", "Mouse V1: disagreement with the 24-trial reference")):
    body.append(f"\\multicolumn{{6}}{{l}}{{\\textit{{{name}}}}} \\\\")
    for X, lab in BENCH:
        if X not in SCORED:
            continue
        vals = []
        for lv in LEVELS:
            e = BM[ds][f"{X} @{lv}%|8"]
            if ds == "EEG":
                vals.append(pct(e["chosen"]) if e["status"] == "estimable" else f"({pct(e['d_deep'])})")
            else:
                vals.append(pct(e["d_deep"], 1))
        body.append(f"\\quad {lab} & " + " & ".join(vals) + " \\\\")
    body.append("\\addlinespace")
write("tab_supp_matched", "\\begin{tabular}{lrrrrr}\n\\toprule\n"
      " & \\multicolumn{5}{c}{Share of points accepted per map} \\\\\n\\cmidrule(lr){2-6}\n"
      "Score & 5\\% & 10\\% & 20\\% & 30\\% & 50\\% \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# Supplementary: EEG blocks from separate sessions
SD = json.load(open(ROOT / "results/session_disjoint/summary.json"))
body = []
for R in (8, 16):
    for design, lab in (("control", "balanced (control)"), ("disjoint", "separate sessions")):
        e = SD[f"{design}|{R}"]
        body.append(f"{R if design == 'control' else ''} & {lab} & " + ", ".join(f"{x:.3f}" for x in e["d"])
                    + f" & {e['chosen']:.3f} & {e['random']:.3f} & {STATUS[e['status']]} \\\\")
    dd = SD[f"difference|{R}"]["smaller"]
    body.append(f" & difference per subject & & {num(dd['median'])} ({num(dd['min'])} to {num(dd['max'])}) & & \\\\")
    body.append("\\addlinespace")
write("tab_supp_session", "\\setlength{\\tabcolsep}{5pt}\n\\begin{tabular}{rllrrl}\n\\toprule\n"
      " & & Disagreement at & Smaller & Random & \\\\\n"
      "Trials & Blocks & $R_h$ = 8, 20, 40 & extrapolant & selection & Status \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# Supplementary: cv-MDS diagnostics
cvd = {ds: BM[ds]["cv_mds_diagnostics"] for ds in ("EEG", "mouse")}
write("tab_supp_cvmds", "\\begin{tabular}{lrrrrr}\n\\toprule\n"
      " & \\multicolumn{3}{c}{Human EEG} & \\multicolumn{2}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}\n"
      " & 8 trials & 16 & 32 & 8 trials & 16 \\\\\n\\midrule\n"
      "Negative cross-validated squared distances (\\%) & " + " & ".join(
          pct(cvd[ds][str(R)]["negative_share"], 1) for ds, R in (("EEG", 8), ("EEG", 16), ("EEG", 32), ("mouse", 8), ("mouse", 16)))
      + " \\\\\nKruskal stress-1 of the MDS fit & " + " & ".join(
          f"{cvd[ds][str(R)]['stress1']:.2f}" for ds, R in (("EEG", 8), ("EEG", 16), ("EEG", 32), ("mouse", 8), ("mouse", 16)))
      + " \\\\\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Numbers quoted in the text
# ---------------------------------------------------------------------------------------------
M = {
    "EegFCeight": rng(TH, "primary", 8, "FC"), "EegLeight": L(TH, "primary", 8, "FC"),
    "EegNFeight": pct(TH["summary"]["primary|8|NF"]), "EegTFeight": rng(TH, "primary", 8, "TF"),
    "EegRsaFCeight": rng(TH, "rsa_cv_mds", 8, "FC"), "EegRsaNFeight": pct(TH["summary"]["rsa_cv_mds|8|NF"], 1),
    "EegCells": str(len(TH["cells_holding"])), "MouseCells": str(len(AL["cells_holding"])),
    "EegBias": f"${TH['bias']['FC']['median_bias']:+.3f}$",
    "MouseNFeight": pct(AL["summary"]["primary|8|NF"]),
    "MouseRsaNFeight": pct(AL["summary"]["rsa_cv_mds|8|NF"], 1),
    "MouseBias": f"${AL['bias']['FC']['median_bias']:+.3f}$",
    "EegUnwhitenedNF": pct(TH["summary"]["unwhitened|8|NF"]),
    "MouseUnwhitenedNF": pct(AL["summary"]["unwhitened|8|NF"]),
    "RatioEeg": f"{S2['things']['ratio']:.2f}", "RatioMouse": f"{S2['allen']['ratio']:.2f}",
    "RelEeg": f"{S2['things']['rel']:.2f}", "RelMouse": f"{S2['allen']['rel']:.2f}",
    "PosMouse": pct(S2["allen"]["pos"]), "MouseNonPos": str(118 - int(S2["allen"]["n_pos8"])),
    "MouseCellsMin": str(min(cells)), "MouseCellsMax": str(max(cells)),
    "MouseMale": str(10 - n_f), "MouseFemale": str(n_f),
    "MouseAgeMin": str(min(m["age_days"] for m in MICE)), "MouseAgeMax": str(max(m["age_days"] for m in MICE)),
    "MouseExperimentIDs": ", ".join(str(m["experiment_id"]) for m in MICE[:-1]) + f" and {MICE[-1]['experiment_id']}",
}
# the smaller extrapolant, the calibrated estimate of C1
c1 = get(TH, "primary", 8, "FC")
M["EegFCeightSmaller"] = pct(min(c1["med_sqrt"], c1["med_inv"]))
M["EegFCeightLarger"] = pct(max(c1["med_sqrt"], c1["med_inv"]))
AU = json.load(open(ROOT / "results/review_audit/audit_summary.json"))
M["EegXtwoEight"] = pct(AU["EEG"]["E5"]["C1"]["x2"])
M["EegQeight"] = rng(TH, "primary", 8, "q", 1)
M["EegQsixteen"] = rng(TH, "primary", 16, "q", 1)
M["MouseQeight"] = rng(AL, "primary", 8, "q", 1)
M["MouseQsixteen"] = rng(AL, "primary", 16, "q", 1)
# finite-reference values, with intervals
fr = lambda ds, key, i=None, d=0: pct(GA[ds][key]["finite_reference"]["median"] if i is None
                                      else GA[ds][key]["finite_reference"]["ci"][i], d)
for tag, key in (("Eight", "primary|8|FC"), ("Sixteen", "primary|16|FC"), ("ThirtyTwo", "primary|32|FC")):
    M[f"EegD{tag}"], M[f"EegD{tag}Lo"], M[f"EegD{tag}Hi"] = fr("EEG", key), fr("EEG", key, 0), fr("EEG", key, 1)
M["EegRsaD"] = fr("EEG", "rsa_cv_mds|8|FC")
for tag, key in (("Eight", "primary|8|FC"), ("Sixteen", "primary|16|FC")):
    M[f"MouseD{tag}"], M[f"MouseD{tag}Lo"], M[f"MouseD{tag}Hi"] = (fr("mouse", key, d=1), fr("mouse", key, 0, 1),
                                                                    fr("mouse", key, 1, 1))
for tag, key in (("Eight", "primary|8|FC"), ("Sixteen", "primary|16|FC")):
    p = GA["mouse"][key]["finite_reference"]["pooled"]
    M[f"MouseConfirmed{tag}"] = f"{p['confirmed']:.0f}"
    M[f"MouseCalls{tag}"] = f"{p['calls']:,.0f}".replace(",", "{,}")
    M[f"MouseConfirmedShare{tag}"] = f"{100 * p['confirmed'] / p['calls']:.2f}" if tag == "Eight" else f"{100 * p['confirmed'] / p['calls']:.1f}"
M["MouseDminSubject"] = pct(min(d_mouse), 1)
# the random comparator and enrichment (results/benchmark/)
M["EegRandEight"] = pct(bm("EEG", "N")["random"])
M["EegEnrichEight"] = f"{bm('EEG', 'N')['enrichment']:.1f}"
for tag, R in (("Eight", 8), ("Sixteen", 16), ("ThirtyTwo", 32)):
    M[f"EegRandDeep{tag}"] = pct(bm("EEG", "N", R)["rand_deep"])
M["EegRsaEnrich"] = f"{bm('EEG', 'cv-MDS N')['enrichment']:.1f}"
M["EegRsaRand"] = pct(bm("EEG", "cv-MDS N")["random"])
M["MouseRandDeep"] = pct(bm("mouse", "N")["rand_deep"], 2)
six = [bm("EEG", X) for X in ("TW", "SCD", "CV", "SEED")]
assert all(e["status"] == "estimable" for e in six)
M["BmEnrichLo"], M["BmEnrichHi"] = f"{min(e['enrichment'] for e in six):.1f}", f"{max(e['enrichment'] for e in six):.1f}"
for X in ("TW", "SCD", "SCD_alt", "CV", "SEED", "SH", "BOOT"):
    e = bm("EEG", X)
    t = X.replace("_alt", "Alt")
    M[f"Bm{t}Pass"] = pct(e["accept"], 1)
    M[f"Bm{t}FC"] = pct(e["chosen"])
    M[f"Bm{t}Enrich"] = f"{e['enrichment']:.1f}"
    M[f"Bm{t}D"] = pct(e["d_deep"])
five = [BM["EEG"][f"{X} @5%|8"] for X in SCORED]
five_ok = [e["chosen"] for e in five if e["status"] == "estimable"]
M["BmFiveLo"], M["BmFiveHi"] = pct(min(five_ok)), pct(max(five_ok))
M["BmFiveUnresolved"] = str(sum(e["status"] != "estimable" for e in five))
fifty = [BM["EEG"][f"{X} @50%|8"]["chosen"] for X in SCORED if X != "cv-MDS N"]
M["BmFiftyLo"], M["BmFiftyHi"] = pct(min(fifty)), pct(max(fifty))
M["BmOverlapAll"], M["BmOverlapAcc"] = f"{bm('EEG', 'N')['overlap_all']:.2f}", f"{bm('EEG', 'N')['overlap_accepted']:.2f}"
M["BmOverlapAllMouse"], M["BmOverlapAccMouse"] = f"{bm('mouse', 'N')['overlap_all']:.2f}", f"{bm('mouse', 'N')['overlap_accepted']:.2f}"
M["BmChanceMouse"], M["BmChanceEeg"] = f"{20 / 117:.2f}", f"{20 / 199:.2f}"
M["BmSHDThirtyTwo"], M["BmNDThirtyTwo"] = pct(bm("EEG", "SH", 32)["d_deep"]), pct(bm("EEG", "N", 32)["d_deep"])
M["BmBOOTPassThirtyTwo"], M["BmBOOTDThirtyTwo"] = pct(bm("EEG", "BOOT", 32)["accept"], 1), pct(bm("EEG", "BOOT", 32)["d_deep"])
M["BmSHPassThirtyTwo"] = pct(bm("EEG", "SH", 32)["accept"], 1)
M["ChkMousePassCV"], M["ChkMousePassSCD"] = pct(bm("mouse", "CV")["accept"], 1), pct(bm("mouse", "SCD")["accept"], 1)
# reproduction of the committed counts
for ds, tag in (("things", "Eeg"), ("allen", "Mouse")):
    rep = CK["datasets"][ds]["reproduction"]
    M[f"Chk{tag}RepRows"] = str(rep["compared"])
biases = [CK["datasets"][ds]["bias"][X]["bias_lo"] for ds in DS for X in CHECKS
          if CK["datasets"][ds]["bias"][X]["n"] >= 5]                   # evaluable in the simulations
M["ChkBiasMin"], M["ChkBiasMax"] = f"${min(biases):+.3f}$", f"${max(biases):+.3f}$"
# the cv-MDS diagnostics
M["CvNegEeg"], M["CvNegMouse"] = pct(cvd["EEG"]["8"]["negative_share"], 1), pct(cvd["mouse"]["8"]["negative_share"])
M["CvStressEeg"], M["CvStressMouse"] = f"{cvd['EEG']['8']['stress1']:.2f}", f"{cvd['mouse']['8']['stress1']:.2f}"
# pooled ratio against the mean over maps (audit item E6), at the deepest reference
e6 = [abs(v["pooled"] - v["mean_of_maps"]) for ds in ("EEG", "mouse") for v in AU[ds]["E6"].values()
      if isinstance(v.get("mean_of_maps"), float) and isinstance(v.get("pooled"), float) and v["pooled"] > 0]
M["PoolDiffMax"] = f"{max(e6):.3f}"
# the pilot's certificates in simulation (audit item E8)
M["CertDraws"] = str(AU["E8"]["need_draws"])
assert sum(AU["E8"]["resolved_any"].values()) == 0
# the session-disjoint analysis
M["SdControl"], M["SdDisjoint"] = pct(SD["control|8"]["chosen"]), pct(SD["disjoint|8"]["chosen"])
M["SdDiff"], M["SdDiffSixteen"] = num(SD["difference|8"]["smaller"]["median"], 2), num(SD["difference|16"]["smaller"]["median"], 2)
# the mouse overlap accounting
OV = json.load(open(ROOT / "results/manifest/mouse_overlap.json"))
M["OvPerPair"] = f"{np.median([o['u_h_per_image_pair'] for o in OV]):.2f}"
assert {o["overlap_frames"]["full"] for o in OV} == {7.0} and {o["window_frames"]["full"] for o in OV} == {15}
assert {o["overlap_frames"]["early"] for o in OV} == {0.0}

# the simulation study, Layer B
def cells_where(ds, pred):
    return [(k, c) for k, c in VAL[ds]["cells"].items() if pred(k, c)]


real = [vcell("EEG", f"T{t}_l1|8") for t in (1, 2, 3)]
M["ValBiasLo"] = num(min(perf(c, "smaller", "median_bias") for c in real))
M["ValBiasHi"] = num(max(perf(c, "smaller", "median_bias") for c in real))
M["ValRmseMax"] = f"{max(perf(c, 'smaller', 'rmse') for c in real):.3f}"
M["ValInvLo"] = f"{min(perf(c, 'larger', 'median_bias') for c in real):.2f}"
M["ValInvHi"] = f"{max(perf(c, 'larger', 'median_bias') for c in real):.2f}"
moderate = [c for ds in ("EEG", "mouse") for k, c in cells_where(ds, lambda k, c: not c["diagnostic"] and k.endswith("|8")
                                                                and "cvmds" not in k and 0.2 <= c["exact"] <= 0.7)]
M["ValModLo"] = f"{min(-perf(c, 'smaller', 'median_bias') for c in moderate):.2f}"
M["ValModHi"] = f"{max(-perf(c, 'smaller', 'median_bias') for c in moderate):.2f}"
M["ValModInv"] = f"{max(abs(perf(c, 'larger', 'median_bias')) for c in moderate):.2f}"
sixteen = [vcell("EEG", f"T{t}_l1|16") for t in (1, 3)]
M["ValSixteenInvLo"] = f"{min(perf(c, 'larger', 'median_bias') for c in sixteen):.2f}"
M["ValSixteenInvHi"] = f"{max(perf(c, 'larger', 'median_bias') for c in sixteen):.2f}"
t32 = [vcell("EEG", f"T{t}_l1|32") for t in (1, 3)]
M["ValThirtyTwoSmLo"] = f"{min(perf(c, 'smaller', 'median_bias') for c in t32):.2f}"
M["ValThirtyTwoSmHi"] = f"{max(perf(c, 'smaller', 'median_bias') for c in t32):.2f}"
M["ValThirtyTwoInvLo"] = f"{min(perf(c, 'larger', 'median_bias') for c in t32):.2f}"
M["ValThirtyTwoInvHi"] = f"{max(perf(c, 'larger', 'median_bias') for c in t32):.2f}"
br = [c["bracket"] for k, c in cells_where("EEG", lambda k, c: not c["diagnostic"] and c["status"] != "reference saturated"
                                                         and "zero" not in k)]
M["ValBracketLo"], M["ValBracketHi"] = pct(min(br)), pct(max(br))
m16 = [vcell("mouse", f"T{t}_l1|16") for t in (1, 3)]
M["ValMouseSixteenLo"] = f"{min(perf(c, 'smaller', 'median_bias') for c in m16):.2f}"
M["ValMouseSixteenHi"] = f"{max(perf(c, 'larger', 'median_bias') for c in m16):.2f}"
M["ValMouseBracketLo"], M["ValMouseBracketHi"] = pct(min(c["bracket"] for c in m16)), pct(max(c["bracket"] for c in m16))
bys = lambda ds, st: VAL[ds]["by_status"][st]["smaller"]
M["ValUnresRmseLo"] = f"{min(bys(ds, 'extrapolation unresolved')['rmse'] for ds in ('EEG', 'mouse')):.2f}"
M["ValUnresRmseHi"] = f"{max(bys(ds, 'extrapolation unresolved')['rmse'] for ds in ('EEG', 'mouse')):.2f}"
M["ValEstRmseEeg"] = f"{bys('EEG', 'estimable')['rmse']:.3f}"
M["ValEstBiasEeg"] = num(bys("EEG", "estimable")["median_bias"])
M["ValSatRmse"] = f"{max(bys(ds, 'reference saturated')['rmse'] for ds in ('EEG', 'mouse')):.3f}"
sat = [c["exact"] for ds in ("EEG", "mouse") for k, c in cells_where(ds, lambda k, c: c["status"] == "reference saturated"
                                                                      and not c["diagnostic"])]
M["ValSatExactMin"] = f"{min(sat):.3f}"
stress = [perf(c, "smaller", "median_bias") for ds in ("EEG", "mouse") for k, c in cells_where(ds, lambda k, c: "_S|" in k)]
assert all(c["status"] != "estimable" for ds in ("EEG", "mouse") for k, c in cells_where(ds, lambda k, c: "_S|" in k))
M["ValStressLo"], M["ValStressHi"] = f"{min(stress):.2f}", f"{max(stress):.2f}"
M["ValDepEeg"] = num(perf(vcell("EEG", "T1_D0.3|8"), "smaller", "median_bias"), 2)
M["ValDepMouse"] = num(perf(vcell("mouse", "T1_D0.3|8"), "smaller", "median_bias"), 2)
qreal = [abs(vcell(ds, f"T{t}_l1|8")["perf_q"]["smaller"]["median_bias"]) for ds in ("EEG", "mouse") for t in (1, 2, 3)]
M["ValQBias"] = f"{max(qreal):.3f}"
qmod = [-c["perf_q"]["smaller"]["median_bias"] for c in moderate]
M["ValQMod"] = f"{max(qmod):.2f}"
# Layer A
la = VAL["EEG"]["layer_a"]["T1_l1"]
M["LaCurveFirst"], M["LaCurveLast"] = f"{la['curve']['8']:.3f}", f"{la['curve']['160']:.3f}"
M["LaExact"], M["LaGrid"], M["LaDeep"] = f"{la['exact']:.3f}", f"{abs(la['grid_sqrt']):.3f}", num(la["deep_sqrt"])
WORD = {"T1": "Tone", "T2": "Ttwo", "T3": "Tthree"}                 # macro names take letters only
for ds, tags in (("EEG", ("T1_S",)), ("mouse", ("T1_S", "T2_S"))):
    for t in tags:
        a = VAL[ds]["layer_a"][t]
        nm = f"La{'Eeg' if ds == 'EEG' else 'Mouse'}{WORD[t.replace('_S', '')]}"
        M[f"{nm}Grid"], M[f"{nm}Deep"] = num(a["grid_sqrt"]), num(a["deep_sqrt"])
# Layer C
for t in ("T1", "T3"):
    o = LC[f"EEG_{t}_tsne_R8"]
    r = o["published"]["L_above"]
    w = WORD[t]
    M[f"LcEeg{w}"], M[f"LcEeg{w}Se"], M[f"LcEeg{w}Allowed"] = f"{r['p']:.3f}", f"{r['mcse']:.3f}", f"{allowed(r):.3f}"
    M[f"LcEeg{w}Theta"] = f"{o['Theta']:.3f}"
    M[f"LcEeg{w}Declares"] = pct(o["published"]["declares_above_half"]["p"], 1)
    M[f"LcEeg{w}Xtwo"] = f"{o['order_statistics']['L_above']['p']:.3f}"
    M[f"LcEeg{w}XtwoDeclares"] = pct(o["order_statistics"]["declares_above_half"]["p"], 1)
    M[f"LcEeg{w}Exact"] = pct(RULE[f"EEG_{t}_tsne_R8"]["FC"]["L_exact_above"]["p"], 1)
    M[f"LcEeg{w}Bias"] = num(RULE[f"EEG_{t}_tsne_R8"]["FC"]["median_bias"])
    M[f"LcEeg{w}Guard"] = str(o["guard"]["triggered"])
    for R in (8, 16):
        tag = f"mouse_{t}_tsne_R{R}"
        r = LC[tag]["mislabel"]["L_above"]
        M[f"LcMouse{w}R{'Eight' if R == 8 else 'Sixteen'}"] = f"{r['p']:.3f}"
        M[f"LcMouse{w}R{'Eight' if R == 8 else 'Sixteen'}Allowed"] = f"{allowed(r):.3f}"
known = RULE["known_populations"]
M["LcRuleLo"], M["LcRuleHi"] = pct(min(v["p"] for v in known.values()), 1), pct(max(v["p"] for v in known.values()), 1)
iqr = []
for t in ("T1",):
    for R in (8, 16):
        acc = {}
        for line in gzip.open(LCDIR / f"oracle_mouse_{t}_tsne_R{R}.jsonl.gz", "rt"):
            r = json.loads(line)
            acc.setdefault(r["v"], np.zeros(2))
            acc[r["v"]] += [r["dis"], r["n"]]
        v = [a[0] / a[1] for a in acc.values()]
        iqr.append(np.percentile(v, 75) - np.percentile(v, 25))
M["LcMouseIqr"] = f"{max(iqr):.3f}"
lo3, hi3 = zip(*(oracle_range(f"mouse_T3_tsne_R{R}", "q") for R in (8, 16)))
M["LcMouseTthreeLo"], M["LcMouseTthreeHi"] = pct(min(lo3)), pct(max(hi3))

# Supplementary numbers: the pilot's noise diagnostics, the simulation study's constants, Layer C details
PILOT = json.load(open(ROOT / "results/things_pilot/diagnostics.json"))
spread = [x["s2_cv_over_homoscedastic"] for x in PILOT]
M["PilotSpread"], M["PilotSpreadLo"], M["PilotSpreadHi"] = (f"{np.median(spread):.2f}", f"{min(spread):.2f}",
                                                            f"{max(spread):.2f}")
ALPHA = json.load(open(ROOT / "results/revision/dev/alpha.json"))
M["AlphaEeg"], M["AlphaMouse"] = f"{ALPHA['EEG']['alpha']:.2f}", f"{ALPHA['mouse']['alpha']:.0f}"
M["GapEeg"], M["GapMouse"] = f"{ALPHA['EEG']['target']:.3f}", f"{ALPHA['mouse']['target']:.3f}"
M["LcThetaSeMax"] = f"{max(max(o['theta_se'].values()) for o in LC.values()):.3f}"
isx3 = [RULE[t][q]["L_is_x3"] for t in RULE if t != "known_populations" for q in ("FC", "q")]
M["LcIsXthreeLo"], M["LcIsXthreeHi"] = pct(min(isx3)), pct(max(isx3))
M["OvShare"] = pct(np.median([o["u_h_share"] for o in OV]))
M["CertValidityDraws"] = f"{AU['E8']['validity_draws']:,}".replace(",", "{,}")

# ---------------------------------------------------------------------------------------------
# Write the tables and numbers into their blocks
# ---------------------------------------------------------------------------------------------
BLOCK = re.compile(r"^% BEGIN make_tables\.py (\S+)\n.*?^% END make_tables\.py \1\n", re.S | re.M)
placed, stale = [], []
for doc in ("manuscript.tex", "supplementary.tex", "cover_letter.tex"):
    old = (OUT / doc).read_text()
    outside = BLOCK.sub("", old)

    def fill(m):
        name = m.group(1)
        placed.append(name)
        if name == "numbers":                       # the macros this file uses, in the order above
            text = "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()
                           if re.search(rf"\\{k}(?![A-Za-z])", outside))
        else:
            assert name in TABLES, f"{doc}: no table called {name}"
            text = TABLES[name]
        return f"% BEGIN make_tables.py {name}\n{text}% END make_tables.py {name}\n"
    new = BLOCK.sub(fill, old)
    if new != old:
        stale.append(doc)
        if CHECK_ONLY:
            sys.stdout.writelines(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                                       doc, doc + " from the results"))
        else:
            (OUT / doc).write_text(new)
    used = set(re.findall(r"\\([A-Za-z]+)", outside))
    missing = [k for k in used if re.fullmatch(r"(Eeg|Mouse|Val|La|Lc|Bm|Chk|Sd|Ov|Cv|Cert|Pool|Ratio|Rel|Pos|Pilot|Alpha|Gap)[A-Z][A-Za-z]*", k)
               and k not in M]
    assert not missing, f"{doc} uses undefined numbers: {sorted(missing)}"
assert sorted(n for n in placed if n != "numbers") == sorted(TABLES), \
    f"each table needs exactly one block: {sorted(set(TABLES) ^ set(n for n in placed if n != 'numbers'))}"
if CHECK_ONLY and stale:
    sys.exit("blocks differ from the committed results in " + ", ".join(stale) + "; run make_tables.py")
print(("rewrote blocks in " + ", ".join(stale)) if stale and not CHECK_ONLY else
      "every table and number matches the committed results")
