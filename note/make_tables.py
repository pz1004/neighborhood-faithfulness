#!/usr/bin/env python3
"""Tables and in-text numbers for note/note.tex, drawn only from committed results.

  tables/tab_primary.tex  Table 1: called faithful, truly faithful, mislabel rate, false comfort
  tables/tab_snr.tex      Table 2: exploratory signal-to-noise (results/allen_confirmatory/explore_snr.txt)
  numbers.tex             every number the text quotes, as macros, so text and tables cannot drift

Sources: results/things_confirmatory/confirmatory_summary.json and
results/allen_confirmatory/allen_summary.json (pre-specified analyses), and
explore_snr.txt (exploratory). Nothing is re-measured here.
"""
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
TH = json.load(open(ROOT / "results/things_confirmatory/confirmatory_summary.json"))
AL = json.load(open(ROOT / "results/allen_confirmatory/allen_summary.json"))


def get(S, cell, R, q):
    return S["summary"][f"{cell}|{R}|{q}"]


def pct(x, d=0):
    return f"{100 * x:.{d}f}"


def rng(S, cell, R, q, d=0):
    """Range spanned by the two extrapolants' medians, in percent."""
    v = get(S, cell, R, q)
    lo, hi = sorted((v["med_sqrt"], v["med_inv"]))
    a, b = pct(lo, d), pct(hi, d)
    return a if a == b else f"{a}--{b}"


def L(S, cell, R, q, d=0):
    return pct(get(S, cell, R, q)["L"], d)


rows = []
for name, S, Rs, digits in (("Human EEG", TH, (8, 16, 32), 0), ("Mouse V1", AL, (8, 16), 1)):
    for pipe, cell in (("t-SNE", "primary"), ("cv-MDS", "rsa_cv_mds")):
        for R in Rs:
            nf = pct(S["summary"][f"{cell}|{R}|NF"], 1)
            rows.append(f"{name} & {pipe} & {R} & {nf} & {rng(S, cell, R, 'TF', 1)} & "
                        f"{rng(S, cell, R, 'q', 1)} & {rng(S, cell, R, 'FC', 1 if name == 'Mouse V1' else 0)} & "
                        f"{L(S, cell, R, 'FC', 1 if name == 'Mouse V1' else 0)} \\\\")
        rows.append("\\addlinespace")
rows = rows[:-1]
(OUT / "tables/tab_primary.tex").write_text(
    "\\begin{tabular}{llrrrrrr}\n\\toprule\n"
    " & & & Called & Truly & Mislabel & False & Bound \\\\\n"
    "Data & Map & Trials & faithful & faithful & rate & comfort & $L$ \\\\\n"
    "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")

snr = (ROOT / "results/allen_confirmatory/explore_snr.txt").read_text()
blocks = dict(zip(("allen", "things"), snr.split("THINGS-EEG2")))
def grab(block, pat):
    return float(re.search(pat, block).group(1))
S2 = {}
for key, b in blocks.items():
    S2[key] = dict(pos=grab(b, r"positive cross-validated d\^2: ([\d.]+)"),
                   rel=grab(b, r"distance matrix: ([\d.]+)"),
                   d8=grab(b, r"8th-neighbor distance ([\d.]+)"),
                   noise=grab(b, r"16-rep mean pair ([\d.]+)"),
                   ratio=grab(b, r"ratio ([\d.]+)"))
(OUT / "tables/tab_snr.tex").write_text(
    "\\begin{tabular}{lrr}\n\\toprule\n & Human EEG & Mouse V1 \\\\\n\\midrule\n"
    f"Image pairs with a positive cross-validated distance & {S2['things']['pos']:.3f} & {S2['allen']['pos']:.3f} \\\\\n"
    f"Split-half reliability of the distance matrix & {S2['things']['rel']:.3f} & {S2['allen']['rel']:.3f} \\\\\n"
    f"Distance to the 8th true neighbor, $d_{{(8)}}$ & {S2['things']['d8']:.2f} & {S2['allen']['d8']:.2f} \\\\\n"
    f"Noise between two 16-trial means, $\\sigma_{{16}}$ & {S2['things']['noise']:.2f} & {S2['allen']['noise']:.2f} \\\\\n"
    f"Ratio $d_{{(8)}}/\\sigma_{{16}}$ & {S2['things']['ratio']:.2f} & {S2['allen']['ratio']:.2f} \\\\\n"
    "\\bottomrule\n\\end{tabular}\n")

def fail_cells(S):
    names = {"k20": "$k=20$", "pca2": "PCA-2"}
    return [c for c in [n for n in ("primary", "rsa_cv_mds", "k5", "k20", "eps0.1", "eps0.3", "tsne_perp15",
                                    "tsne_perp50", "pca2", "mds", "early_window", "unwhitened")]
            if c not in S["cells_holding"]]

tf_max_mouse = max(max(v["med_sqrt"], v["med_inv"]) for k, v in AL["summary"].items()
                   if k.endswith("|TF"))
M = {
    "EegFCeight": rng(TH, "primary", 8, "FC"), "EegLeight": L(TH, "primary", 8, "FC"),
    "EegFCsixteen": rng(TH, "primary", 16, "FC"), "EegLsixteen": L(TH, "primary", 16, "FC"),
    "EegFCthirtytwo": rng(TH, "primary", 32, "FC"), "EegLthirtytwo": L(TH, "primary", 32, "FC"),
    "EegNFeight": pct(TH["summary"]["primary|8|NF"]), "EegTFeight": rng(TH, "primary", 8, "TF"),
    "EegQeight": rng(TH, "primary", 8, "q"), "EegQL": "{}\\%, {}\\% and {}\\%".format(*(L(TH, "primary", R, "q") for R in (8, 16, 32))),
    "EegRsaFCeight": rng(TH, "rsa_cv_mds", 8, "FC"), "EegRsaLeight": L(TH, "rsa_cv_mds", 8, "FC"),
    "EegRsaNFeight": pct(TH["summary"]["rsa_cv_mds|8|NF"], 1),
    "EegCells": str(len(TH["cells_holding"])), "EegBias": f"${TH['bias']['FC']['median_bias']:+.3f}$",
    "EegKtwentyFC": rng(TH, "k20", 8, "FC"), "EegKtwentyL": L(TH, "k20", 8, "FC"),
    "EegPcaFC": rng(TH, "pca2", 8, "FC"), "EegPcaL": L(TH, "pca2", 8, "FC"),
    "MouseFCeight": rng(AL, "primary", 8, "FC"), "MouseLeight": L(AL, "primary", 8, "FC", 1),
    "MouseLsixteen": L(AL, "primary", 16, "FC", 1), "MouseNFeight": pct(AL["summary"]["primary|8|NF"]),
    "MouseQeight": rng(AL, "primary", 8, "q"), "MouseQL": "{}\\% and {}\\%".format(*(L(AL, "primary", R, "q") for R in (8, 16))),
    "MouseRsaFCeight": rng(AL, "rsa_cv_mds", 8, "FC"), "MouseRsaLeight": L(AL, "rsa_cv_mds", 8, "FC", 1),
    "MouseRsaNFeight": pct(AL["summary"]["rsa_cv_mds|8|NF"], 1),
    "MouseCells": str(len(AL["cells_holding"])), "MouseBias": f"${AL['bias']['FC']['median_bias']:+.3f}$",
    "MouseTFmax": pct(tf_max_mouse, 1), "MouseUnwhitenedNF": pct(AL["summary"]["unwhitened|8|NF"]),
    "EegUnwhitenedFC": rng(TH, "unwhitened", 8, "FC"),
    "RatioEeg": f"{S2['things']['ratio']:.2f}", "RatioMouse": f"{S2['allen']['ratio']:.2f}",
    "RelEeg": f"{S2['things']['rel']:.2f}", "RelMouse": f"{S2['allen']['rel']:.2f}",
    "PosMouse": pct(S2["allen"]["pos"]),
}
assert fail_cells(TH) == ["k20", "pca2"] and fail_cells(AL) == [], (fail_cells(TH), fail_cells(AL))
(OUT / "numbers.tex").write_text("% generated by make_tables.py from committed results; do not edit\n"
                                 + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in M.items()))
print((OUT / "numbers.tex").read_text())
print((OUT / "tables/tab_primary.tex").read_text())
print((OUT / "tables/tab_snr.tex").read_text())
