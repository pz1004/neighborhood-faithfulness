#!/usr/bin/env python3
"""Tables and in-text numbers for manuscript/manuscript.tex and supplementary.tex, drawn only from
committed results. Nothing is re-measured here.

  tables/tab_design.tex   Table 1: datasets and design
  tables/tab_primary.tex  Table 2: called faithful, truly faithful, mislabel rate, false comfort
  tables/tab_checks.tex   Table 3: common checks against the truth at 8 trials
  tables/tab_snr.tex      Table 4: exploratory signal-to-noise
  tables/tab_s*.tex       Supplementary Tables S1-S6
  numbers.tex             every number the text quotes, as macros, so text and tables cannot drift

Sources:
- results/things_confirmatory/{rows.jsonl, confirmatory_summary.json};
- results/allen_confirmatory/{rows.jsonl, allen_summary.json, explore_snr.txt,
  subjects_described.json};
- results/checks/checks_summary.json.
"""
import json, re, sys
from pathlib import Path

import numpy as np

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


def write(name, text):
    (OUT / "tables" / name).write_text(text)


# ---------------------------------------------------------------------------------------------
# Table 1: design
# ---------------------------------------------------------------------------------------------
cells = [m["n_cells"] for m in MICE]
n_f = sum(m["sex"] == "F" for m in MICE)
write("tab_design.tex", "\\begin{tabular}{p{0.26\\linewidth}p{0.32\\linewidth}p{0.32\\linewidth}}\n\\toprule\n"
      " & Human EEG (THINGS-EEG2) & Mouse V1 (Allen Brain Observatory) \\\\\n\\midrule\n"
      "Subjects & 10 adults (8 female, 2 male) & "
      f"10 mice ({10 - n_f} male, {n_f} female), one experiment each \\\\\n"
      "Stimuli & 200 object images & 118 natural scenes \\\\\n"
      "Repetitions per image & 80, over 4 sessions & 50, early and late strata of 25 \\\\\n"
      "Features per trial & 17 channels $\\times$ 80 samples (0--0.79\\,s) & "
      f"{min(cells)}--{max(cells)} cells, mean $\\Delta F/F$ over 0.5\\,s \\\\\n"
      "Blocks $F$ / $U$ / $H$ & 8 / 32 / 40 & 10 / 16 / 24 \\\\\n"
      "Analyst's budgets $R$ & 8, 16, 32 & 8, 16 \\\\\n"
      "Reference sizes $R_h$ (sub-blocks) & 8 (5), 20 (2), 40 (1) & 8 (3), 12 (2), 24 (1) \\\\\n"
      "Random splits per subject & 10 & 10 \\\\\n"
      "Representation, fitted on $F$ & \\multicolumn{2}{l}{PCA to 100, Ledoit--Wolf whitening, PCA to $D=20$} \\\\\n"
      "\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Table 2: primary results
# ---------------------------------------------------------------------------------------------
rows = []
for name, S, Rs in (("Human EEG", TH, (8, 16, 32)), ("Mouse V1", AL, (8, 16))):
    d = 1 if name == "Mouse V1" else 0
    for pipe, cell in (("t-SNE", "primary"), ("cv-MDS", "rsa_cv_mds")):
        for R in Rs:
            rows.append(f"{name} & {pipe} & {R} & {pct(S['summary'][f'{cell}|{R}|NF'], 1)} & "
                        f"{rng(S, cell, R, 'TF', 1)} & {rng(S, cell, R, 'q', 1)} & "
                        f"{rng(S, cell, R, 'FC', d)} & {L(S, cell, R, 'FC', d)} \\\\")
        rows.append("\\addlinespace")
write("tab_primary.tex", "\\begin{tabular}{llrrrrrr}\n\\toprule\n"
      " & & & Called & Truly & Mislabel & False & Bound \\\\\n"
      "Data & Map & Trials & faithful & faithful & rate & comfort & $L$ \\\\\n"
      "\\midrule\n" + "\n".join(rows[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Table 4: exploratory signal-to-noise
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
write("tab_snr.tex", "\\begin{tabular}{lrr}\n\\toprule\n & Human EEG & Mouse V1 \\\\\n\\midrule\n"
      f"Image pairs with a positive cross-validated distance & {S2['things']['pos']:.3f} & {S2['allen']['pos']:.3f} \\\\\n"
      f"Split-half reliability of the distance matrix & {S2['things']['rel']:.3f} & {S2['allen']['rel']:.3f} \\\\\n"
      f"Distance to the 8th nearest image, $d_{{(8)}}$ & {S2['things']['d8']:.2f} & {S2['allen']['d8']:.2f} \\\\\n"
      f"Noise between two 16-trial means, $\\sigma_{{16}}$ & {S2['things']['noise']:.2f} & {S2['allen']['noise']:.2f} \\\\\n"
      f"Ratio $d_{{(8)}}/\\sigma_{{16}}$ & {S2['things']['ratio']:.2f} & {S2['allen']['ratio']:.2f} \\\\\n"
      "\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Supplementary Tables S1, S2: every cell and budget
# ---------------------------------------------------------------------------------------------
for tab, S, Rs, d in (("tab_s1_cells_eeg.tex", TH, (8, 16, 32), 1), ("tab_s2_cells_mouse.tex", AL, (8, 16), 1)):
    body = []
    for cell in CELL_ORDER:
        for R in Rs:
            body.append(f"{CELL_LABEL[cell] if R == Rs[0] else ''} & {R} & {pct(S['summary'][f'{cell}|{R}|NF'], 1)} & "
                        f"{rng(S, cell, R, 'TF', 1)} & {rng(S, cell, R, 'q', 1)} & {rng(S, cell, R, 'FC', d)} & "
                        f"{L(S, cell, R, 'FC', d)} \\\\")
        body.append("\\addlinespace")
    write(tab, "\\begin{tabular}{lrrrrrr}\n\\toprule\n"
          " & & Called & Truly & Mislabel & False & Bound \\\\\n"
          "Variation & Trials & faithful & faithful & rate & comfort & $L$ \\\\\n\\midrule\n"
          + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Supplementary Table S3: per-subject estimates, primary cell, 8 trials
# ---------------------------------------------------------------------------------------------
def per_subject(path, mod):
    """Per-subject estimates of the primary cell at 8 trials, with the registered estimator."""
    rr = [r for r in map(json.loads, open(path)) if r["cell"] == "primary" and r["R"] == 8]
    return mod.estimates(mod.disagreement(rr))


est_th = per_subject(ROOT / "results/things_confirmatory/rows.jsonl", AT)
est_al = per_subject(ROOT / "results/allen_confirmatory/rows.jsonl", AA)
body = []
for s in range(1, 11):
    a, b = est_th[(s, "primary", 8)], est_al[(s, "primary", 8)]
    fmt = lambda e: f"{pct(e['NF'], 1)} & {pct(e['FC'][0], 1)} & {pct(e['FC'][1], 1)}"
    body.append(f"{s} & {fmt(a)} & {fmt(b)} & {MICE[s - 1]['n_cells']} \\\\")
write("tab_s3_subjects.tex", "\\begin{tabular}{rrrrrrrr}\n\\toprule\n"
      " & \\multicolumn{3}{c}{Human EEG} & \\multicolumn{4}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\n"
      "Subject & NF & FC, $R_h^{-1/2}$ & FC, $1/R_h$ & NF & FC, $R_h^{-1/2}$ & FC, $1/R_h$ & Cells \\\\\n\\midrule\n"
      + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Supplementary Table S4: plasmode bias per subject, primary cell, 8 trials (confirmatory check)
# ---------------------------------------------------------------------------------------------
body = []
for s in range(10):
    body.append(f"{s + 1} & {TH['bias']['FC']['per_subject'][s]:+.3f} & {TH['bias']['q']['per_subject'][s]:+.3f} & "
                f"{AL['bias']['FC']['per_subject'][s]:+.3f} & {AL['bias']['q']['per_subject'][s]:+.3f} \\\\")
body.append("\\midrule")
body.append(f"Median & {TH['bias']['FC']['median_bias']:+.3f} & {TH['bias']['q']['median_bias']:+.3f} & "
            f"{AL['bias']['FC']['median_bias']:+.3f} & {AL['bias']['q']['median_bias']:+.3f} \\\\")
write("tab_s4_plasmode.tex", "\\begin{tabular}{rrrrr}\n\\toprule\n"
      " & \\multicolumn{2}{c}{Human EEG} & \\multicolumn{2}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n"
      "Subject & FC & $q$ & FC & $q$ \\\\\n\\midrule\n" + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

# ---------------------------------------------------------------------------------------------
# Numbers quoted in the text
# ---------------------------------------------------------------------------------------------
tf_max_mouse = max(max(v["med_sqrt"], v["med_inv"]) for k, v in AL["summary"].items() if k.endswith("|TF"))
M = {
    "EegFCeight": rng(TH, "primary", 8, "FC"), "EegLeight": L(TH, "primary", 8, "FC"),
    "EegFCsixteen": rng(TH, "primary", 16, "FC"), "EegLsixteen": L(TH, "primary", 16, "FC"),
    "EegFCthirtytwo": rng(TH, "primary", 32, "FC"), "EegLthirtytwo": L(TH, "primary", 32, "FC"),
    "EegNFeight": pct(TH["summary"]["primary|8|NF"]), "EegTFeight": rng(TH, "primary", 8, "TF"),
    "EegQL": "{}\\%, {}\\% and {}\\%".format(*(L(TH, "primary", R, "q") for R in (8, 16, 32))),
    "EegRsaFCeight": rng(TH, "rsa_cv_mds", 8, "FC"), "EegRsaLeight": L(TH, "rsa_cv_mds", 8, "FC"),
    "EegRsaNFeight": pct(TH["summary"]["rsa_cv_mds|8|NF"], 1),
    "EegCells": str(len(TH["cells_holding"])), "EegBias": f"${TH['bias']['FC']['median_bias']:+.3f}$",
    "EegKtwentyFC": rng(TH, "k20", 8, "FC"), "EegKtwentyL": L(TH, "k20", 8, "FC"),
    "EegPcaFC": rng(TH, "pca2", 8, "FC"), "EegPcaL": L(TH, "pca2", 8, "FC"),
    "MouseFCeight": rng(AL, "primary", 8, "FC"), "MouseLeight": L(AL, "primary", 8, "FC", 1),
    "MouseLsixteen": L(AL, "primary", 16, "FC", 1), "MouseNFeight": pct(AL["summary"]["primary|8|NF"]),
    "MouseQL": "{}\\% and {}\\%".format(*(L(AL, "primary", R, "q") for R in (8, 16))),
    "MouseRsaFCeight": rng(AL, "rsa_cv_mds", 8, "FC"), "MouseRsaLeight": L(AL, "rsa_cv_mds", 8, "FC", 1),
    "MouseRsaNFeight": pct(AL["summary"]["rsa_cv_mds|8|NF"], 1),
    "MouseCells": str(len(AL["cells_holding"])), "MouseBias": f"${AL['bias']['FC']['median_bias']:+.3f}$",
    "MouseTFmax": pct(tf_max_mouse, 1), "MouseUnwhitenedNF": pct(AL["summary"]["unwhitened|8|NF"]),
    "EegUnwhitenedFC": rng(TH, "unwhitened", 8, "FC"),
    "RatioEeg": f"{S2['things']['ratio']:.2f}", "RatioMouse": f"{S2['allen']['ratio']:.2f}",
    "RelEeg": f"{S2['things']['rel']:.2f}", "RelMouse": f"{S2['allen']['rel']:.2f}",
    "PosMouse": pct(S2["allen"]["pos"]), "MouseNonPos": str(118 - int(S2["allen"]["n_pos8"])),
    "MouseCellsMin": str(min(cells)), "MouseCellsMax": str(max(cells)),
    "MouseMale": str(10 - n_f), "MouseFemale": str(n_f),
}
assert [c for c in CELL_ORDER if c not in TH["cells_holding"]] == ["k20", "pca2"]
assert [c for c in CELL_ORDER if c not in AL["cells_holding"]] == []

# ---------------------------------------------------------------------------------------------
# The common checks (Table 3, Supplementary Tables S5 and S6, and their numbers)
# ---------------------------------------------------------------------------------------------
CK = json.load(open(ROOT / "results/checks/checks_summary.json"))
DS = {"things": ("Human EEG", (8, 16, 32)), "allen": ("Mouse V1", (8, 16))}
SHORT = {"leaves most false comfort": "leaves", "removes most false comfort": "removes",
         "undetermined": "undetermined", "not evaluable": "not evaluable"}


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


body = []
for ds, (name, Rs) in DS.items():
    body.append(f"\\multicolumn{{7}}{{l}}{{\\textit{{{name}}}}} \\\\")
    for X in CHECKS:
        cls = "--" if X == "N" else SHORT[CK["classification"][ds][X]]
        catch = "--" if X == "N" else ck_rng(ds, X, 8, "CATCH")
        withn = "--" if X == "N" else ck_rng(ds, X, 8, "FC_NX", 1 if ds == "allen" else 0)
        d = 1 if ds == "allen" else 0
        body.append(f"\\quad {CHECK_LABEL[X]} & {ck_pass(ds, X, 8)} & {ck_rng(ds, X, 8, 'FC', d)} & "
                    f"{ck_bound(ds, X, 8, 'L', d)}, {ck_bound(ds, X, 8, 'U', d)} & {catch} & {withn} & {cls} \\\\")
    body.append("\\addlinespace")
write("tab_checks.tex", "\\begin{tabular}{lrrrrrl}\n\\toprule\n"
      " & & False & Bounds & Plug-in false & With & \\\\\n"
      "Check & Passes & comfort & $L$, $U$ & comfort flagged & plug-in & Classification \\\\\n\\midrule\n"
      + "\n".join(body[:-1]) + "\n\\bottomrule\n\\end{tabular}\n")

for tab, ds in (("tab_s5_checks_eeg.tex", "things"), ("tab_s6_checks_mouse.tex", "allen")):
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
write("tab_s7_checks_plasmode.tex", "\\begin{tabular}{lrrrrrr}\n\\toprule\n"
      " & \\multicolumn{3}{c}{Human EEG} & \\multicolumn{3}{c}{Mouse V1} \\\\\n"
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
      "Check & Smaller & Larger & Subjects & Smaller & Larger & Subjects \\\\\n\\midrule\n"
      + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n")

biases = [CK["datasets"][ds]["bias"][X]["bias_lo"] for ds in DS for X in CHECKS
          if CK["datasets"][ds]["bias"][X]["n"] >= 5]                   # evaluable in the simulations
M["ChkBiasMin"], M["ChkBiasMax"] = f"${min(biases):+.3f}$", f"${max(biases):+.3f}$"
for ds, tag in (("things", "Eeg"), ("allen", "Mouse")):
    rep = CK["datasets"][ds]["reproduction"]
    M[f"Chk{tag}RepRows"], M[f"Chk{tag}RepDiff"] = str(rep["compared"]), str(rep["differing"])
    for X in CHECKS:
        M[f"Chk{tag}{X}Pass"] = ck_pass(ds, X, 8)
        d = 1 if ds == "allen" else 0
        M[f"Chk{tag}{X}FC"] = ck_rng(ds, X, 8, "FC", d)
        M[f"Chk{tag}{X}L"] = ck_bound(ds, X, 8, "L", d)
        M[f"Chk{tag}{X}U"] = ck_bound(ds, X, 8, "U", d)
        M[f"Chk{tag}{X}Catch"] = ck_rng(ds, X, 8, "CATCH")
        M[f"Chk{tag}{X}Yield"] = ck_rng(ds, X, 8, "YIELD", 1)
        M[f"Chk{tag}{X}FCNX"] = ck_rng(ds, X, 8, "FC_NX")
        M[f"Chk{tag}{X}PassNX"] = ck_pass(ds, X, 8, "PASS_NX")
M["ChkEegSHPassNXThirtyTwo"] = ck_pass("things", "SH", 32, "PASS_NX")
M["ChkEegSHFCNXThirtyTwo"] = ck_rng("things", "SH", 32, "FC_NX")
M["ChkEegSHFCThirtyTwo"] = ck_rng("things", "SH", 32, "FC")
M["ChkEegSHPassThirtyTwo"] = ck_pass("things", "SH", 32)
M["ChkEegBOOTPassThirtyTwo"] = ck_pass("things", "BOOT", 32)
M["ChkEegBOOTUThirtyTwo"] = ck_bound("things", "BOOT", 32, "U")
M["ChkEegBOOTPassSixteen"] = ck_pass("things", "BOOT", 16, d=2)
M["ChkMouseMinL"] = pct(min(ck("allen", X, 8, "FC")["L"] for X in ("TW", "SCD", "CV", "SEED")), 1)
M["ChkEegMinL"] = pct(min(ck("things", X, 8, "FC")["L"] for X in ("TW", "SCD", "CV", "SEED")))
M["ChkEegExtraSCD"] = pct(ck("things", "SCD", 8, "PASS")["med"] - ck("things", "SCD", 8, "PASS_NX")["med"])
M["ChkEegExtraSEED"] = pct(ck("things", "SEED", 8, "PASS")["med"] - ck("things", "SEED", 8, "PASS_NX")["med"])

(OUT / "numbers.tex").write_text("% generated by make_tables.py from committed results; do not edit\n"
                                 + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
print((OUT / "numbers.tex").read_text())
