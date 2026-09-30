#!/usr/bin/env python3
"""Figures for note/note.tex, drawn only from committed results. Nothing is re-measured here.

  Figure 1  (a) called faithful vs truly faithful at 8 trials; (b) false comfort vs trials
  Figure 2  false comfort at 8 trials in the 12 pre-specified variations, both datasets

Sources: results/things_confirmatory/confirmatory_summary.json and
results/allen_confirmatory/allen_summary.json.

**Estimates shown.** Each estimate is the midpoint of the two extrapolants' medians. Bars run
from the conservative bound L up to the upper end of the wider 95% interval.

**Palette** (as in the predecessor's report): blue #2a78d6 for EEG, orange #eb6834 for mouse
V1. Marker shape and line style carry the map type, so identity never rests on hue alone.
The PDFs carry no creation date, so unchanged results rebuild the same files.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "figures"
TH = json.load(open(ROOT / "results/things_confirmatory/confirmatory_summary.json"))["summary"]
AL = json.load(open(ROOT / "results/allen_confirmatory/allen_summary.json"))["summary"]

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE, REF = "#e1e0d9", "#c3c2b7", "#ffffff", "#a8a69e"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2,
    "axes.titlesize": 9, "axes.titlecolor": INK, "axes.titlelocation": "left", "axes.titlepad": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.size": 0, "ytick.major.size": 0,
    "legend.frameon": False, "legend.fontsize": 7.5,
    "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "pdf.fonttype": 42,
})
pct = FuncFormatter(lambda v, _: f"{v:.0%}")
NO_DATE = {"CreationDate": None}
MARK = dict(markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.0)


def est(S, cell, R, q):
    v = S[f"{cell}|{R}|{q}"]
    lo, hi = sorted((v["med_sqrt"], v["med_inv"]))
    return (lo + hi) / 2, v["L"], max(v["ci_sqrt"][1], v["ci_inv"][1]), lo, hi


def fig1():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 3.05), gridspec_kw=dict(width_ratios=[1.05, 1]))
    groups = [("EEG\nt-SNE", TH, "primary", BLUE), ("EEG\ncv-MDS", TH, "rsa_cv_mds", BLUE),
              ("mouse\nt-SNE", AL, "primary", ORANGE), ("mouse\ncv-MDS", AL, "rsa_cv_mds", ORANGE)]
    x = np.arange(len(groups))
    for j, (lab, S, cell, col) in enumerate(groups):
        nf = S[f"{cell}|8|NF"]
        tf, _, _, tlo, thi = est(S, cell, 8, "TF")
        a.bar(j - 0.19, nf, 0.36, color=col, alpha=0.9, label="called faithful" if j == 0 else None)
        a.bar(j + 0.19, tf, 0.36, color=SURFACE, edgecolor=col, linewidth=1.2,
              label="truly faithful" if j == 0 else None)
        a.plot([j + 0.19] * 2, [max(tlo, 0), thi], color=INK2, linewidth=1.0)
    a.set_xticks(x, [g[0] for g in groups])
    a.yaxis.set_major_formatter(pct)
    a.set_ylim(0, 0.22)
    a.set_title("a   Share of points, 8 trials")
    a.legend(loc="upper right")
    a.grid(axis="x", visible=False)

    series = [("EEG, t-SNE", TH, "primary", (8, 16, 32), BLUE, "o", "-"),
              ("EEG, cv-MDS", TH, "rsa_cv_mds", (8, 16, 32), BLUE, "s", "--"),
              ("mouse V1, t-SNE", AL, "primary", (8, 16), ORANGE, "o", "-"),
              ("mouse V1, cv-MDS", AL, "rsa_cv_mds", (8, 16), ORANGE, "s", "--")]
    for k, (lab, S, cell, Rs, col, mk, ls) in enumerate(series):
        off = (k - 1.5) * 0.9
        xs = np.array(Rs, float) * (1 + off / 40)
        pts = [est(S, cell, R, "FC") for R in Rs]
        mid = [p[0] for p in pts]
        b.errorbar(xs, mid, yerr=[[m - p[1] for m, p in zip(mid, pts)], [p[2] - m for m, p in zip(mid, pts)]],
                   color=col, marker=mk, linestyle=ls, linewidth=1.3, elinewidth=0.9, capsize=0,
                   label=lab, **MARK)
    b.axhline(0.5, color=REF, linewidth=0.9, linestyle=":")
    b.set_xscale("log", base=2)
    b.set_xticks([8, 16, 32], ["8", "16", "32"])
    b.set_xlim(6.5, 40)
    b.set_ylim(0, 1.04)
    b.yaxis.set_major_formatter(pct)
    b.set_xlabel("trials averaged")
    b.set_title("b   False comfort")
    b.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=7, columnspacing=1.2, handlelength=3.2)
    fig.tight_layout(w_pad=2.0)
    fig.savefig(OUT / "fig1_false_comfort.pdf", metadata=NO_DATE)
    plt.close(fig)


CELLS = [("primary", "primary (t-SNE, perplexity 30)"), ("rsa_cv_mds", "cv-MDS (MDS of cross-validated distances)"),
         ("k5", "k = 5"), ("k20", "k = 20"), ("eps0.1", "ε = 0.1"), ("eps0.3", "ε = 0.3"),
         ("tsne_perp15", "t-SNE perplexity 15"), ("tsne_perp50", "t-SNE perplexity 50"),
         ("pca2", "PCA-2"), ("mds", "MDS of Euclidean distances"),
         ("early_window", "alternative response window"), ("unwhitened", "unwhitened representation")]


def fig2():
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    y = np.arange(len(CELLS))[::-1]
    for S, col, mk, dy, lab in ((TH, BLUE, "o", 0.16, "human EEG"), (AL, ORANGE, "s", -0.16, "mouse V1")):
        for yi, (cell, _) in zip(y, CELLS):
            m, Lb, up, _, _ = est(S, cell, 8, "FC")
            ax.plot([Lb, up], [yi + dy] * 2, color=col, linewidth=1.1, alpha=0.8)
            ax.plot([m], [yi + dy], marker=mk, color=col, linestyle="none", **MARK,
                    label=lab if cell == "primary" else None)
            ax.plot([Lb], [yi + dy], marker="|", color=col, markersize=7)
    ax.axvline(0.5, color=REF, linewidth=0.9, linestyle=":")
    ax.set_yticks(y, [c[1] for c in CELLS])
    ax.set_xlim(0, 1.03)
    ax.xaxis.set_major_formatter(pct)
    ax.set_xlabel("false comfort at 8 trials (tick: conservative bound L)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_robustness.pdf", metadata=NO_DATE)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig1()
    fig2()
    print("figures written to", OUT)
