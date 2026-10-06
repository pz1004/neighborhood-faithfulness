#!/usr/bin/env python3
"""Figures for manuscript/manuscript.tex and supplementary.tex, drawn only from committed results.
Nothing is re-measured.

  Figure 1   the estimator: design, the reach label, and the disagreement curves in both datasets
  Figure 2   where the estimator can be trusted, in simulation (Layers A, B and C)
  Figure 3   the real data: disagreement with references of increasing depth, by budget
  Figure 4   common checks at 8 trials: false comfort against the share of points passed
  Figure S1  the 12 pre-specified variations at 8 trials

Sources:
- results/things_confirmatory/{rows.jsonl, confirmatory_summary.json};
- results/allen_confirmatory/{rows.jsonl, allen_summary.json};
- results/revision/validation/summary.json; results/revision/layer_c/calibration_*.json;
- results/gate_a/classification.json; results/benchmark/summary.json.

**What is drawn as an estimate.** A latent estimate (an extrapolated false comfort) is drawn only
where the result's status is estimable; everything else is drawn as the disagreement with the
deepest reference, which is reportable whatever the status (manuscript, Section 2.3.4).

**Palette.** Blue #2a78d6 for EEG, orange #eb6834 for mouse V1. Marker shape and line style
carry the status, the budget or the check, so identity never rests on hue alone. The PDFs carry
no creation date, so unchanged results rebuild the same files.
"""
import json, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import analyze_things_confirmatory as AT      # noqa: E402  (the registered estimators)
import analyze_allen_confirmatory as AA       # noqa: E402

OUT = Path(__file__).resolve().parent / "figures"
TH = json.load(open(ROOT / "results/things_confirmatory/confirmatory_summary.json"))["summary"]
AL = json.load(open(ROOT / "results/allen_confirmatory/allen_summary.json"))["summary"]
VAL = json.load(open(ROOT / "results/revision/validation/summary.json"))
GA = json.load(open(ROOT / "results/gate_a/classification.json"))
BM = json.load(open(ROOT / "results/benchmark/summary.json"))
LCDIR = ROOT / "results/revision/layer_c"

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, SURFACE, REF = "#e1e0d9", "#c3c2b7", "#ffffff", "#a8a69e"
FILL = "#f3f2ee"
COLOR = {"EEG": BLUE, "mouse": ORANGE}
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
STATUS_MARK = {"estimable": ("o", True), "extrapolation unresolved": ("^", False),
               "reference saturated": ("s", True), "insufficient acceptance": ("x", False)}
STATUS_LABEL = {"estimable": "estimable", "extrapolation unresolved": "unresolved",
                "reference saturated": "saturated", "insufficient acceptance": "too few passes"}


def rh_axis(ax, sizes, xmax=None):
    """An axis linear in R_h^(-1/2), with the reference sizes as ticks and infinity at 0."""
    x = [h ** -0.5 for h in sizes]
    ax.set_xticks([0] + x, ["∞"] + [str(h) for h in sizes])
    ax.set_xlim(-0.015, (xmax or max(x)) * 1.06)
    ax.set_xlabel("reference trials $R_h$")


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
    lane = [(0.55, ["block $U$: the\nanalyst's $R$ trials", "averages $c_i$", "t-SNE map;\ndisplayed\nneighbors $E_i$",
                    "naive label $N_i$:\nis $E_i$ faithful\nin the averages?"]),
            (0.17, ["block $H$: disjoint\nsub-blocks of\n$R_h$ trials", "sub-block\nmeans", "the same $E_i$",
                    "reference label\n$T_i(R_h)$: is $E_i$\nfaithful there?"])]
    for y, texts in lane:
        for j, (x, t) in enumerate(zip(xs, texts)):
            _box(ax, x, y, w, h, t, ec=BLUE if j == 3 else AXIS, lw=1.1 if j == 3 else 0.8,
                 fc=SURFACE if (y < 0.3 and j == 2) else FILL)
            if j:
                _arrow(ax, (xs[j - 1] + w, y + h / 2), (x, y + h / 2))
    _arrow(ax, (0.41 + w / 2, 0.55), (0.41 + w / 2, 0.17 + h), color=MUTED, ls="--")
    _box(ax, 0.83, 0.17, 0.17, 0.60,
         "disagreement\n$d_{FC}(R_h)$: share of\n$N_i = 1$ with\n$T_i(R_h) = 0$\n\nextrapolated to\n"
         "$R_h \\to \\infty$ where\nestimable (panel c)", fc=SURFACE, ec=INK2, size=7.0)
    _arrow(ax, (0.615 + w, 0.55 + h / 2), (0.83, 0.55 + h / 2))
    _arrow(ax, (0.615 + w, 0.17 + h / 2), (0.83, 0.17 + h / 2))
    ax.text(0.0, 0.03, "repetitions per image:  $F + R + H$ = 8 + $R$ + 40 of 80 (EEG);  10 + $R$ + 24 of 50 (mouse V1)",
            fontsize=7.0, color=INK2, transform=ax.transAxes, ha="left", va="center")


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


def _curves(ax, rows_path, mod, RH, color, title, extrapolate):
    rows = [json.loads(l) for l in open(rows_path)]
    rows = [r for r in rows if r["cell"] == "primary" and r["R"] == 8]
    dis = mod.disagreement(rows) if mod is AT else mod.disagreement(rows, RH)
    x = np.array([h ** -0.5 for h in RH])
    xf = np.linspace(0, x.max() * 1.03, 60)
    e1s, e2s = [], []
    for key in sorted(dis):
        y = np.array([dis[key][h][1] for h in RH])
        ax.plot(x, y, color=color, lw=0.8, alpha=0.6, marker="o", markersize=2.6, markeredgewidth=0)
        if extrapolate:
            b1, a1 = np.polyfit(x, y, 1)
            b2, a2 = np.polyfit(x ** 2, y, 1)
            e1s.append(a1)
            e2s.append(a2)
            ax.plot(xf, a1 + b1 * xf, color=color, lw=0.6, alpha=0.35, ls="--")
            ax.plot(xf, a2 + b2 * xf ** 2, color=color, lw=0.6, alpha=0.35, ls=":")
    rh_axis(ax, RH)
    ax.yaxis.set_major_formatter(pct)
    ax.set_ylim(0.6, 1.02)
    ax.set_title(title)
    if not extrapolate:
        ax.text(0.04, 0.63, "saturated: not extrapolated", fontsize=7, color=MUTED, transform=ax.transData)
        return None
    m1, m2 = np.median(e1s), np.median(e2s)
    ax.plot([0], [m1], marker="o", color=INK, markersize=5.5, markeredgecolor=SURFACE, zorder=5, clip_on=False)
    ax.plot([0], [m2], marker="^", color=INK, markersize=5.5, markeredgecolor=SURFACE, zorder=5, clip_on=False)
    return m1, m2


def fig_method():
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1.0], width_ratios=[1.1, 1, 1], hspace=0.38, wspace=0.42,
                          left=0.03, right=0.99)
    _design(fig.add_subplot(gs[0, :]))
    _reach(fig.add_subplot(gs[1, 0]))
    c = fig.add_subplot(gs[1, 1])
    m1, m2 = _curves(c, ROOT / "results/things_confirmatory/rows.jsonl", AT, AT.RH, BLUE, "c   EEG, 8 trials", True)
    c.set_ylabel("disagreement $d_{FC}(R_h)$")
    c.annotate(f"{m1:.2f}", (0, m1), xytext=(6, -3), textcoords="offset points", fontsize=7, color=INK, va="top")
    c.annotate(f"{m2:.2f}", (0, m2), xytext=(6, 3), textcoords="offset points", fontsize=7, color=INK, va="bottom")
    d = fig.add_subplot(gs[1, 2])
    _curves(d, ROOT / "results/allen_confirmatory/rows.jsonl", AA, AA.RH, ORANGE, "d   Mouse V1, 8 trials", False)
    d.legend(handles=[Line2D([], [], color=INK2, lw=0.8, marker="o", markersize=2.6, label="one subject"),
                      Line2D([], [], color=INK2, lw=0.8, ls="--", label="fit linear in $R_h^{-1/2}$"),
                      Line2D([], [], color=INK2, lw=0.8, ls=":", label="fit linear in $1/R_h$"),
                      Line2D([], [], color=INK, lw=0, marker="o", markersize=5, label="median, $R_h^{-1/2}$"),
                      Line2D([], [], color=INK, lw=0, marker="^", markersize=5, label="median, $1/R_h$")],
             loc="center left", bbox_to_anchor=(0.02, 0.55), fontsize=6.3, handlelength=2.2)
    fig.savefig(OUT / "fig1_method.pdf", metadata=NO_DATE, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------
# Figure 2: the estimator in simulation
# ------------------------------------------------------------------------------------------

def fig_validation():
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1.25, 1.0, 1.0], wspace=0.5))
    # (a) Layer B
    clip = 0.26
    off_scale = []
    a.axhspan(-0.03, 0.03, color=FILL, zorder=0)
    a.axhline(0, color=REF, lw=0.7)
    for ds in ("EEG", "mouse"):
        for key, cell in VAL[ds]["cells"].items():
            bias = cell["perf"]["smaller"]["median_bias"]
            if not np.isfinite(bias):
                continue
            y = float(np.clip(bias, -clip, clip))
            if cell["diagnostic"]:
                mk = "x" if "_S|" in key else "D"
                a.plot([cell["exact"]], [y], marker=mk, color=MUTED, linestyle="none", markersize=4.5,
                       markeredgewidth=1.0, fillstyle="none" if mk == "D" else "full")
                if abs(bias) > clip:
                    off_scale.append(bias)
                continue
            mk, filled = STATUS_MARK[cell["status"]]
            a.plot([cell["exact"]], [y], marker=mk, color=COLOR[ds], linestyle="none",
                   markersize=4.8, fillstyle="full" if filled else "none", markeredgewidth=1.0,
                   markeredgecolor=COLOR[ds] if not filled else SURFACE)
    a.set_xlim(-0.03, 1.03)
    a.set_ylim(-clip - 0.02, clip + 0.02)
    a.text(0.2, clip - 0.005, "off scale: " + ", ".join(f"{b:+.2f}" for b in sorted(off_scale)), fontsize=6.2,
           color=MUTED, va="top")
    a.xaxis.set_major_formatter(pct)
    a.set_xlabel("exact false comfort")
    a.set_ylabel("median bias, smaller extrapolant")
    a.set_title("a   Layer B: 100 datasets per cell")
    handles = [Line2D([], [], color=BLUE, marker="o", lw=0, markersize=4.8, markeredgecolor=SURFACE, label="EEG"),
               Line2D([], [], color=ORANGE, marker="o", lw=0, markersize=4.8, markeredgecolor=SURFACE, label="mouse V1"),
               Line2D([], [], color=INK2, marker="o", lw=0, markersize=4.8, markeredgecolor=SURFACE, label="estimable"),
               Line2D([], [], color=INK2, marker="^", lw=0, markersize=4.8, fillstyle="none", label="unresolved"),
               Line2D([], [], color=INK2, marker="s", lw=0, markersize=4.8, markeredgecolor=SURFACE, label="saturated"),
               Line2D([], [], color=MUTED, marker="x", lw=0, markersize=4.5, label="stress test"),
               Line2D([], [], color=MUTED, marker="D", lw=0, markersize=4.0, fillstyle="none", label="shared noise")]
    a.legend(handles=handles, loc="lower left", fontsize=6.0, ncol=2, handletextpad=0.3, columnspacing=0.8)
    # (b) Layer A
    series = [("EEG", "T1_l1", (8, 20, 40, 80, 160), "-", "EEG, data's noise"),
              ("EEG", "T1_S", (8, 20, 40, 80, 160), "--", "EEG, stress test"),
              ("mouse", "T1_S", (8, 12, 24, 48, 96), "--", "mouse, stress T1"),
              ("mouse", "T2_S", (8, 12, 24, 48, 96), ":", "mouse, stress T2")]
    for ds, key, grid, ls, lab in series:
        la = VAL[ds]["layer_a"][key]
        x = [h ** -0.5 for h in grid]
        y = [la["curve"][str(h)] for h in grid]
        b.plot(x, y, color=COLOR[ds], ls=ls, lw=1.2, marker="o", markersize=3, markeredgewidth=0, label=lab)
        b.plot([0], [la["exact"]], marker="*", color=COLOR[ds], markersize=7, clip_on=False, markeredgecolor=SURFACE)
    b.axvspan(40 ** -0.5 * 0.999, 8 ** -0.5 * 1.06, color=FILL, zorder=0)
    rh_axis(b, (8, 20, 40, 160))
    b.tick_params(axis="x", labelsize=7)
    b.set_ylim(0, 1.02)
    b.yaxis.set_major_formatter(pct)
    b.set_ylabel("disagreement $d_{FC}(R_h)$")
    b.set_title("b   Layer A: deeper references")
    b.text(0.255, 0.03, "real sizes", fontsize=6.3, color=MUTED, ha="center")
    b.legend(loc="lower left", bbox_to_anchor=(0.0, 0.08), fontsize=6.0, handlelength=2.4)
    # (c) Layer C
    runs = [("EEG_T1_tsne_R8", "FC", "EEG, T1"), ("EEG_T3_tsne_R8", "FC", "EEG, T3"),
            ("mouse_T1_tsne_R8", "q", "mouse 8, T1"), ("mouse_T3_tsne_R8", "q", "mouse 8, T3"),
            ("mouse_T1_tsne_R16", "q", "mouse 16, T1"), ("mouse_T3_tsne_R16", "q", "mouse 16, T3")]
    y = np.arange(len(runs))[::-1]
    for yi, (tag, q, lab) in zip(y, runs):
        o = json.load(open(LCDIR / f"calibration_{tag}.json"))
        r = o["published"]["L_above"] if q == "FC" else o["mislabel"]["L_above"]
        col = COLOR[tag.split("_")[0]]
        allowed = 0.025 + 2 * r["mcse"]
        c.plot([r["p"] - 2 * r["mcse"], r["p"] + 2 * r["mcse"]], [yi, yi], color=col, lw=1.2)
        c.plot([r["p"]], [yi], marker="o" if r["p"] <= allowed else "X", color=col, linestyle="none", **MARK)
        c.plot([allowed], [yi], marker="|", color=INK2, markersize=9, markeredgewidth=1.4)
    c.axvline(0.025, color=REF, lw=0.9, ls=":")
    c.set_yticks(y, [r[2] for r in runs], fontsize=7)
    c.set_xlim(0, 0.09)
    c.set_xticks([0, 0.025, 0.05, 0.075], ["0", "0.025", "0.05", "0.075"])
    c.set_xlabel("$\\Pr(L > \\Theta)$ over 1,000 studies")
    c.set_title("c   Layer C")
    c.grid(axis="y", visible=False)
    c.legend(handles=[Line2D([], [], color=INK2, marker="o", lw=0, markersize=5, markeredgecolor=SURFACE, label="met"),
                      Line2D([], [], color=INK2, marker="X", lw=0, markersize=5, markeredgecolor=SURFACE, label="missed"),
                      Line2D([], [], color=INK2, marker="|", lw=0, markersize=8, markeredgewidth=1.4, label="allowed")],
             loc="lower right", fontsize=6.3, handletextpad=0.3)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.9)
    fig.savefig(OUT / "fig2_validation.pdf", metadata=NO_DATE)
    plt.close(fig)


# ------------------------------------------------------------------------------------------
# Figure 3: the real data, by budget
# ------------------------------------------------------------------------------------------

def fig_depth():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(wspace=0.3))
    style = {8: ("o", "-"), 16: ("s", "--"), 32: ("^", ":")}
    for ax, ds, RH, Rs, title in ((axes[0], "EEG", (8, 20, 40), (8, 16, 32), "a   Human EEG"),
                                  (axes[1], "mouse", (8, 12, 24), (8, 16), "b   Mouse V1")):
        col = COLOR[ds]
        x = np.array([h ** -0.5 for h in RH])
        for R in Rs:
            mk, ls = style[R]
            e = BM[ds][f"N|{R}"]
            g = GA[ds][f"primary|{R}|FC"]
            fr = g["finite_reference"]
            st = STATUS_LABEL[g["status"]]
            ax.plot(x, e["d"], color=col, marker=mk, ls=ls, lw=1.2, label=f"{R} trials ({st})", **MARK)
            dx = {8: 0.007, 16: 0.0, 32: -0.007}[R]
            ax.plot([x[-1] + dx] * 2, fr["ci"], color=col, lw=1.0, alpha=0.7)
            ax.plot([x[-1] - 0.02], [e["rand_deep"]], marker=mk, color=col, fillstyle="none", linestyle="none",
                    markersize=5, markeredgewidth=1.0)
            if ds == "EEG" and R == 8:
                S = TH["primary|8|FC"]
                lo, hi = sorted((S["med_sqrt"], S["med_inv"]))
                ax.plot([0, x[-1]], [lo, e["d"][-1]], color=col, ls="--", lw=0.9, alpha=0.8)
                ax.plot([0, 0], [S["L"], hi], color=INK, lw=1.4, clip_on=False)
                ax.plot([0], [lo], marker="o", color=INK, markersize=5.5, markeredgecolor=SURFACE, zorder=5, clip_on=False)
                ax.annotate(f"estimate {lo:.2f}\nbound $L$ {S['L']:.2f}", (0, S["L"]), xytext=(7, -2),
                            textcoords="offset points", fontsize=6.6, color=INK, va="top")
        rh_axis(ax, RH)
        ax.set_ylim(0.45 if ds == "EEG" else 0.9, 1.01)
        ax.yaxis.set_major_formatter(pct)
        ax.set_title(title)
        ax.legend(loc="lower right", fontsize=6.6)
    axes[0].set_ylabel("disagreement of the passing\nneighborhoods with the reference")
    axes[1].text(0.03, 0.955, "open markers: random selection,\nat the deepest reference", fontsize=6.6, color=MUTED)
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.16, top=0.9)
    fig.savefig(OUT / "fig3_depth.pdf", metadata=NO_DATE)
    plt.close(fig)


# ------------------------------------------------------------------------------------------
# Figure 4: common checks at 8 trials
# ------------------------------------------------------------------------------------------

NATIVE = [("N", "N", "o"), ("TW", "TW", "o"), ("SCD", "SCD", "o"), ("SCD_alt", "SCD alt", "o"), ("CV", "CV", "s"),
          ("SEED", "SEED", "^"), ("SH", "SH", "^"), ("BOOT", "BOOT", "^"), ("cv-MDS N", "cv-MDS", "s")]
LABEL_AT = {"N": (-9, 4), "TW": (4, -7), "SCD": (-20, 2), "SCD_alt": (4, -6), "CV": (-12, 4), "SEED": (2, 4),
            "SH": (4, 2), "BOOT": (4, 2), "cv-MDS N": (4, 2)}
CURVES = [("N", "-"), ("TW", "--"), ("SCD", "-."), ("CV", ":"), ("BOOT", (0, (1, 2))), ("cv-MDS N", (0, (5, 2)))]
LEVELS = (5, 10, 20, 30, 50)


def fig_checks():
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), gridspec_kw=dict(wspace=0.36))
    panels = ((axes[0], "EEG", "chosen", "a   EEG: estimate", "random"),
              (axes[1], "EEG", "d_deep", "b   EEG: 40-trial reference", "rand_deep"),
              (axes[2], "mouse", "d_deep", "c   Mouse V1: 24-trial reference", "rand_deep"))
    for ax, ds, key, title, rkey in panels:
        col = COLOR[ds]
        latent = key == "chosen"
        for X, ls in CURVES:
            pts = [BM[ds][f"{X} @{lv}%|8"] for lv in LEVELS]
            ok = [p for p in pts if (p["status"] == "estimable" or not latent) and np.isfinite(p[key])]
            if ok:
                ax.plot([p["accept"] for p in ok], [p[key] for p in ok], color=col, ls=ls, lw=0.9, alpha=0.75)
        for X, lab, mk in NATIVE:
            e = BM[ds][f"{X}|8"]
            if e["accept"] <= 0 or not np.isfinite(e[key]):
                continue
            estimable = e["status"] == "estimable"
            if latent and not estimable:
                continue
            ax.plot([e["accept"]], [e[key]], marker=mk, color=col, linestyle="none",
                    fillstyle="full" if estimable else "none", markersize=5.5,
                    markeredgecolor=SURFACE if estimable else col, markeredgewidth=1.0)
            if ds == "EEG":
                ax.annotate(lab, (e["accept"], e[key]), xytext=LABEL_AT[X], textcoords="offset points", fontsize=6.0,
                            color=INK2)
        rand = np.nanmedian([BM[ds][f"{X}|8"][rkey] for X, _, _ in NATIVE if np.isfinite(BM[ds][f"{X}|8"][rkey])])
        ax.axhline(rand, color=MUTED, lw=1.0)
        ax.text(0.0055, rand + 0.006, "random selection", fontsize=6.2, color=MUTED, va="bottom")
        ax.set_xscale("log")
        ax.set_xlim(0.005, 0.7)
        ax.set_xticks([0.01, 0.05, 0.2, 0.5], ["1%", "5%", "20%", "50%"])
        ax.set_ylim((0.4, 1.0) if ds == "EEG" else (0.98, 1.001))
        ax.yaxis.set_major_formatter(pct if ds == "EEG" else FuncFormatter(lambda v, _: f"{v:.1%}"))
        ax.set_xlabel("share of points passed")
        ax.set_title(title, fontsize=8.3)
    axes[2].text(0.006, 0.9845, "every check is saturated\nat the 24-trial reference", fontsize=6.4, color=MUTED)
    axes[0].legend(handles=[Line2D([], [], color=INK2, lw=0.9, ls=ls, label=lab) for (X, ls), lab in
                            zip(CURVES, ("N", "TW", "SCD", "CV", "BOOT", "cv-MDS"))],
                   title="matched acceptance", title_fontsize=6.2, loc="lower right", fontsize=6.0, ncol=2,
                   handlelength=2.2, columnspacing=0.8)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.17, top=0.9)
    fig.savefig(OUT / "fig4_checks.pdf", metadata=NO_DATE)
    plt.close(fig)


# ------------------------------------------------------------------------------------------
# Figure S1: the pre-specified variations
# ------------------------------------------------------------------------------------------

CELLS = [("primary", "primary (t-SNE, perplexity 30)"), ("rsa_cv_mds", "cv-MDS (MDS of cross-validated distances)"),
         ("k5", "k = 5"), ("k20", "k = 20"), ("eps0.1", "ε = 0.1"), ("eps0.3", "ε = 0.3"),
         ("tsne_perp15", "t-SNE perplexity 15"), ("tsne_perp50", "t-SNE perplexity 50"),
         ("pca2", "PCA-2"), ("mds", "MDS of Euclidean distances"),
         ("early_window", "alternative response window"), ("unwhitened", "unwhitened representation")]


def fig_variations():
    fig, ax = plt.subplots(figsize=(5.8, 3.6))
    y = np.arange(len(CELLS))[::-1]
    for ds, S, mk, dy in (("EEG", TH, "o", 0.16), ("mouse", AL, "s", -0.16)):
        col = COLOR[ds]
        for yi, (cell, _) in zip(y, CELLS):
            g = GA[ds][f"{cell}|8|FC"]
            if g["status"] == "estimable":
                v = S[f"{cell}|8|FC"]
                lo, hi = sorted((v["med_sqrt"], v["med_inv"]))
                ax.plot([lo, hi], [yi + dy] * 2, color=col, lw=1.1, alpha=0.85)
                ax.plot([lo], [yi + dy], marker=mk, color=col, linestyle="none", **MARK)
            else:
                fr = g["finite_reference"]
                ax.plot(fr["ci"], [yi + dy] * 2, color=col, lw=1.0, alpha=0.6)
                ax.plot([fr["median"]], [yi + dy], marker=mk, color=col, linestyle="none", fillstyle="none",
                        markersize=5.5, markeredgewidth=1.0)
    ax.axvline(0.5, color=REF, lw=0.9, ls=":")
    ax.set_yticks(y, [c[1] for c in CELLS])
    ax.set_xlim(0.3, 1.03)
    ax.xaxis.set_major_formatter(pct)
    ax.set_xlabel("false comfort at 8 trials: estimate (filled)\nor disagreement with the deepest reference (open)")
    ax.grid(axis="y", visible=False)
    ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", lw=0, markersize=5, markeredgecolor=SURFACE, label="human EEG"),
                       Line2D([], [], color=ORANGE, marker="s", lw=0, markersize=5, markeredgecolor=SURFACE, label="mouse V1")],
              loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "figS1_variations.pdf", metadata=NO_DATE)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    fig_method()
    fig_validation()
    fig_depth()
    fig_checks()
    fig_variations()
    print("figures written to", OUT)
