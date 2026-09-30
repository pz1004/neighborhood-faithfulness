#!/usr/bin/env python3
"""Figures for manuscript/manuscript.tex, drawn only from committed results. Nothing is re-measured.

  Figure 1  the estimator: design, the reach label, and the extrapolation in both datasets
  Figure 2  called faithful vs truly faithful at 8 trials, and false comfort vs trials
  Figure 3  false comfort at 8 trials in the 12 pre-specified variations
  Figure 4  common checks against the truth at 8 trials

Sources:
- results/things_confirmatory/{rows.jsonl, confirmatory_summary.json};
- results/allen_confirmatory/{rows.jsonl, allen_summary.json};
- results/checks/checks_summary.json.

**Estimates shown.** Each estimate is the midpoint of the two extrapolants' medians. Bars run
from the conservative bound L to the upper end of the wider 95% interval, or to U for the
checks.

**Palette.** Blue #2a78d6 for EEG, orange #eb6834 for mouse V1. Marker shape and line style
carry the map type or the check family, so identity never rests on hue alone. The PDFs carry
no creation date, so unchanged results rebuild the same files.
"""
import json, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import analyze_things_confirmatory as AT      # noqa: E402  (the registered estimators)
import analyze_allen_confirmatory as AA       # noqa: E402

OUT = Path(__file__).resolve().parent / "figures"
TH = json.load(open(ROOT / "results/things_confirmatory/confirmatory_summary.json"))["summary"]
AL = json.load(open(ROOT / "results/allen_confirmatory/allen_summary.json"))["summary"]
CHECKS = ROOT / "results/checks/checks_summary.json"

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE, REF = "#e1e0d9", "#c3c2b7", "#ffffff", "#a8a69e"
FILL = "#f3f2ee"
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


# ------------------------------------------------------------------------------------------
# Figure 1: the estimator
# ------------------------------------------------------------------------------------------

def _box(ax, x, y, w, h, text, fc=FILL, ec=AXIS, size=7.2, color=INK, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.003,rounding_size=0.015",
                                fc=fc, ec=ec, lw=lw, transform=ax.transAxes))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color=color,
            transform=ax.transAxes, linespacing=1.3)


def _arrow(ax, a, b, color=INK2, ls="-"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=8, lw=0.9, color=color,
                                 linestyle=ls, transform=ax.transAxes, shrinkA=0, shrinkB=0))


def _design(ax):
    ax.set_axis_off()
    ax.set_title("a   The estimator, for each subject and random split of the repetitions", pad=4)
    _box(ax, 0.0, 0.86, 0.79, 0.12,
         "block $F$ fits the representation used by both lanes: PCA, whitening, PCA to 20 dimensions",
         fc=SURFACE, size=7.0)
    w, h, xs = 0.17, 0.22, (0.0, 0.205, 0.41, 0.615)
    lane = [(0.52, ["block $U$: the\nanalyst's $R$ trials", "averages $c_i$", "t-SNE map;\ndisplayed\nneighbors $E_i$",
                    "naive label $N_i$:\nis $E_i$ faithful\nin the averages?"]),
            (0.10, ["block $H$: disjoint\nsub-blocks of\n$R_h$ trials", "sub-block\nmeans", "the same $E_i$",
                    "reference label\n$T_i(R_h)$: is $E_i$\nfaithful there?"])]
    for y, texts in lane:
        for j, (x, t) in enumerate(zip(xs, texts)):
            _box(ax, x, y, w, h, t, ec=BLUE if j == 3 else AXIS, lw=1.1 if j == 3 else 0.8,
                 fc=SURFACE if (y < 0.3 and j == 2) else FILL)
            if j:
                _arrow(ax, (xs[j - 1] + w, y + h / 2), (x, y + h / 2))
    _arrow(ax, (0.41 + w / 2, 0.52), (0.41 + w / 2, 0.10 + h), color=MUTED, ls="--")
    _box(ax, 0.83, 0.10, 0.17, 0.64,
         "disagreement\n$d_{FC}(R_h)$: share of\n$N_i = 1$ with\n$T_i(R_h) = 0$\n\nextrapolated to\n"
         "$R_h \\to \\infty$ gives\nthe false comfort\n(panels c, d)", fc=SURFACE, ec=INK2, size=7.0)
    _arrow(ax, (0.615 + w, 0.52 + h / 2), (0.83, 0.52 + h / 2))
    _arrow(ax, (0.615 + w, 0.10 + h / 2), (0.83, 0.10 + h / 2))


def _reach(ax):
    ax.set_title("b   The reach label (schematic)")
    k, m, kp = 10, 8, 20
    rows = [("ranks in the averages", np.array([1, 2, 3, 4, 5, 6, 8, 11, 14, 17])),
            ("ranks in the noise-free means", np.array([2, 5, 9, 13, 18, 24, 29, 33, 41, 52]))]
    for yi, (lab, rk) in zip((1.0, 0.0), rows):
        rk = np.sort(rk)
        reach = rk[m - 1]
        ok = reach <= kp
        col = BLUE if ok else ORANGE
        ax.scatter(rk, np.full(k, yi), s=15, color=INK2, zorder=3, lw=0)
        ax.scatter([reach], [yi], s=62, facecolor="none", edgecolor=col, lw=1.4, zorder=4)
        ax.text(1, yi + 0.3, lab, fontsize=6.9, color=MUTED, ha="left")
        ax.text(reach, yi - 0.36, f"reach {reach}: " + ("faithful" if ok else "not faithful"),
                ha="left" if ok else "center", fontsize=7, color=col, weight="bold")
    ax.axvline(kp, color=REF, lw=0.9, ls=":")
    ax.text(kp + 0.8, 1.55, "$k' = 20$", fontsize=7, color=MUTED)
    ax.set_xlim(0, 56)
    ax.set_ylim(-0.75, 1.75)
    ax.set_yticks([])
    ax.set_xticks([1, 10, 20, 30, 40, 50])
    ax.set_xlabel("rank of a displayed neighbor")
    ax.grid(False)
    ax.spines["left"].set_visible(False)


def _curves(ax, rows_path, mod, RH, color, title):
    rows = [json.loads(l) for l in open(rows_path)]
    rows = [r for r in rows if r["cell"] == "primary" and r["R"] == 8]
    dis = mod.disagreement(rows) if mod is AT else mod.disagreement(rows, RH)
    x = np.array([h ** -0.5 for h in RH])
    xf = np.linspace(0, x.max() * 1.03, 60)
    e1s, e2s = [], []
    for key in sorted(dis):
        y = np.array([dis[key][h][1] for h in RH])
        b1, a1 = np.polyfit(x, y, 1)
        b2, a2 = np.polyfit(x ** 2, y, 1)
        e1s.append(a1)
        e2s.append(a2)
        ax.plot(x, y, color=color, lw=0.8, alpha=0.6, marker="o", markersize=2.6, markeredgewidth=0)
        ax.plot(xf, a1 + b1 * xf, color=color, lw=0.6, alpha=0.35, ls="--")
        ax.plot(xf, a2 + b2 * xf ** 2, color=color, lw=0.6, alpha=0.35, ls=":")
    m1, m2 = np.median(e1s), np.median(e2s)
    ax.plot([0], [m1], marker="o", color=INK, markersize=5.5, markeredgecolor=SURFACE, zorder=5, clip_on=False)
    ax.plot([0], [m2], marker="^", color=INK, markersize=5.5, markeredgecolor=SURFACE, zorder=5, clip_on=False)
    ax.set_xticks([0] + list(x), ["∞"] + [str(h) for h in RH])
    ax.set_xlim(-0.015, x.max() * 1.06)
    ax.set_xlabel("reference trials $R_h$")
    ax.yaxis.set_major_formatter(pct)
    ax.set_ylim(0.6, 1.02)
    ax.set_title(title)
    return m1, m2


def fig_method():
    fig = plt.figure(figsize=(7.2, 5.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], width_ratios=[1.1, 1, 1], hspace=0.38, wspace=0.42,
                          left=0.03, right=0.99)
    _design(fig.add_subplot(gs[0, :]))
    _reach(fig.add_subplot(gs[1, 0]))
    c = fig.add_subplot(gs[1, 1])
    m1, m2 = _curves(c, ROOT / "results/things_confirmatory/rows.jsonl", AT, AT.RH, BLUE, "c   EEG, 8 trials")
    c.set_ylabel("disagreement $d_{FC}(R_h)$")
    c.annotate(f"{m1:.2f}", (0, m1), xytext=(6, -3), textcoords="offset points", fontsize=7, color=INK, va="top")
    c.annotate(f"{m2:.2f}", (0, m2), xytext=(6, 3), textcoords="offset points", fontsize=7, color=INK, va="bottom")
    d = fig.add_subplot(gs[1, 2])
    _curves(d, ROOT / "results/allen_confirmatory/rows.jsonl", AA, AA.RH, ORANGE, "d   Mouse V1, 8 trials")
    from matplotlib.lines import Line2D
    d.legend(handles=[Line2D([], [], color=INK2, lw=0.8, marker="o", markersize=2.6, label="one subject"),
                      Line2D([], [], color=INK2, lw=0.8, ls="--", label="fit linear in $R_h^{-1/2}$"),
                      Line2D([], [], color=INK2, lw=0.8, ls=":", label="fit linear in $1/R_h$"),
                      Line2D([], [], color=INK, lw=0, marker="o", markersize=5, label="median, $R_h^{-1/2}$"),
                      Line2D([], [], color=INK, lw=0, marker="^", markersize=5, label="median, $1/R_h$")],
             loc="lower right", fontsize=6.5, handlelength=2.2)
    fig.savefig(OUT / "fig1_method.pdf", metadata=NO_DATE, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------
# Figures 2 and 3: as in the note
# ------------------------------------------------------------------------------------------

def fig_false_comfort():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.15), gridspec_kw=dict(width_ratios=[1.05, 1]))
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
    fig.savefig(OUT / "fig2_false_comfort.pdf", metadata=NO_DATE)
    plt.close(fig)


CELLS = [("primary", "primary (t-SNE, perplexity 30)"), ("rsa_cv_mds", "cv-MDS (MDS of cross-validated distances)"),
         ("k5", "k = 5"), ("k20", "k = 20"), ("eps0.1", "ε = 0.1"), ("eps0.3", "ε = 0.3"),
         ("tsne_perp15", "t-SNE perplexity 15"), ("tsne_perp50", "t-SNE perplexity 50"),
         ("pca2", "PCA-2"), ("mds", "MDS of Euclidean distances"),
         ("early_window", "alternative response window"), ("unwhitened", "unwhitened representation")]


def fig_robustness():
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
    fig.savefig(OUT / "fig3_robustness.pdf", metadata=NO_DATE)
    plt.close(fig)


# ------------------------------------------------------------------------------------------
# Figure 4: common checks against the truth
# ------------------------------------------------------------------------------------------

CHECK_ROWS = [("N", "plug-in (N)", "o"), ("TW", "trustworthiness", "o"), ("SCD", "scDEED", "o"),
              ("CV", "cross-validated", "s"), ("SEED", "seed recurrence", "^"),
              ("BOOT", "bootstrap recurrence", "^"), ("SH", "split-half recurrence", "^")]


def fig_checks():
    CK = json.load(open(CHECKS))["datasets"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True, gridspec_kw=dict(wspace=0.32))
    y = np.arange(len(CHECK_ROWS))[::-1]
    for ax, ds, col, title in ((axes[0], "things", BLUE, "a   Human EEG, 8 trials"),
                               (axes[1], "allen", ORANGE, "b   Mouse V1, 8 trials")):
        S = CK[ds]["summary"]
        for yi, (X, lab, mk) in zip(y, CHECK_ROWS):
            f, p = S[f"{X}|8|FC"], S[f"{X}|8|PASS"]
            if f["n_defined"] < 5:
                ax.text(0.02, yi, "not evaluable", va="center", fontsize=7, color=MUTED)
            else:
                mid = (f["med_sqrt"] + f["med_inv"]) / 2
                ax.plot([f["L"], f["U"]], [yi, yi], color=col, lw=1.2, alpha=0.85)
                ax.plot([mid], [yi], marker=mk, color=col, linestyle="none", **MARK)
            ax.text(1.05, yi, f"{100 * p['med']:.1f}%", va="center", ha="left", fontsize=7, color=INK2,
                    transform=ax.get_yaxis_transform(), clip_on=False)
        ax.axvline(0.5, color=REF, lw=0.9, ls=":")
        ax.set_xlim(0, 1.02)
        ax.xaxis.set_major_formatter(pct)
        ax.set_xlabel("false comfort of the check (bar: $L$ to $U$)")
        ax.set_title(title)
        ax.grid(axis="y", visible=False)
        ax.text(1.05, len(CHECK_ROWS) - 0.35, "passes", ha="left", fontsize=7, color=MUTED,
                transform=ax.get_yaxis_transform(), clip_on=False)
    axes[0].set_yticks(y, [r[1] for r in CHECK_ROWS])
    fig.subplots_adjust(left=0.2, right=0.93, bottom=0.16, top=0.9)
    fig.savefig(OUT / "fig4_checks.pdf", metadata=NO_DATE, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig_method()
    fig_false_comfort()
    fig_robustness()
    fig_checks()
    print("figures written to", OUT)
