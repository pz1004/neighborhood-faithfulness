#!/usr/bin/env python3
"""Apply the pre-registered decision of plan Amendment 2 to the THINGS-EEG2 pilot results.

Reads results/things_pilot/{records.jsonl, semi_synthetic.jsonl, diagnostics.json}; writes
results/things_pilot/pilot_analysis.txt and pilot_summary.json. Every statistic is a median
over the 10 subjects; bootstrap intervals resample subjects (2000 resamples, seed 3_700_000).
"""
import json, collections
import numpy as np

OUT = "results/things_pilot"
RU, RH = (8, 16, 32), (8, 20, 40)
KP = "20"
lines = []


def out(s=""):
    lines.append(s)
    print(s)


def line_fit(x, y):
    A = np.column_stack([np.ones_like(x), x])
    return np.linalg.lstsq(A, y, rcond=None)[0]              # intercept, slope


def boot_median(v, rng, B=2000):
    v = np.asarray(v, float)
    meds = [np.median(rng.choice(v, len(v))) for _ in range(B)]
    return float(np.median(v)), float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


recs = [json.loads(l) for l in open(f"{OUT}/records.jsonl")]
sim = [json.loads(l) for l in open(f"{OUT}/semi_synthetic.jsonl")]
diag = json.load(open(f"{OUT}/diagnostics.json"))
subs = sorted({r["subject"] for r in recs})
rng = np.random.default_rng(3_700_000)
tsne = lambda r: r["display"].startswith("tsne")

# ---------------------------------------------------------------- primary need estimator
xh = np.array([1 / np.sqrt(R) for R in RH])
q = {R: {} for R in RU}
d40 = {R: {} for R in RU}
for s in subs:
    for R in RU:
        rs = [r for r in recs if r["subject"] == s and r["R"] == R and tsne(r)]
        d = np.array([np.mean([r["disagree"][KP][str(Rh)] for r in rs]) for Rh in RH])
        q[R][s] = float(line_fit(xh, d)[0])
        d40[R][s] = float(d[-1])
q80 = {}
xr = np.array([1 / np.sqrt(R) for R in RU])
for s in subs:
    a, b = line_fit(xr, np.array([q[R][s] for R in RU]))
    q80[s] = float(a + b / np.sqrt(80))

out("THINGS-EEG2 NEED PILOT (plan Amendment 2): k=10, eps=0.2, k'=20, t-SNE displays")
out("Primary: reference-extrapolated naive mislabel rate q(R), median over subjects [95% bootstrap]")
summary = {}
for R in RU:
    m, lo, hi = boot_median(list(q[R].values()), rng)
    summary[f"q_{R}"] = (m, lo, hi)
    out(f"  q({R:2d}) = {m:.3f} [{lo:.3f}, {hi:.3f}]   raw disagreement with R_h=40: "
        f"{np.median(list(d40[R].values())):.3f}   per subject: "
        + " ".join(f"{q[R][s]:.3f}" for s in subs))
m80, lo80, hi80 = boot_median(list(q80.values()), rng)
summary["q_80"] = (m80, lo80, hi80)
out(f"  q(80), extrapolated with q = a + b/sqrt(R): {m80:.3f} [{lo80:.3f}, {hi80:.3f}]   per subject: "
    + " ".join(f"{q80[s]:.3f}" for s in subs))

q32 = summary["q_32"][0]
if q32 < 0.05:
    need = "FAILS (q(32) < 0.05, so q(80) < 0.05)"
elif m80 >= 0.05:
    need = "PASSES (q(80) >= 0.05)"
else:
    need = "INCONCLUSIVE (q(32) >= 0.05 > q(80))"
out(f"  need at native R: {need}")
out(f"  need at R=8: {'met' if summary['q_8'][0] >= 0.05 else 'not met'} (q(8) = {summary['q_8'][0]:.3f})")

# ---------------------------------------------------------------- secondary need estimators
out("\nSecondary: certified bounds from E-pair-t on H (R=40), medians over subjects")
for R in RU:
    rows = collections.defaultdict(list)
    for s in subs:
        rs = [r for r in recs if r["subject"] == s and r["R"] == R and tsne(r)]
        for f in ("lo", "hi", "resolved"):
            rows[f].append(np.mean([r["certified_bounds"][KP][f] for r in rs]))
    out(f"  R={R:2d}: naive mislabel in [{np.median(rows['lo']):.3f}, {np.median(rows['hi']):.3f}]"
        f"  (H resolves {np.median(rows['resolved']):.3f} of points)")

out("\nSecondary: semi-synthetic (truth = H means, noise = resampled U residuals), medians")
need_sim = [r for r in sim if r["kind"] == "need"]
for R in RU:
    per = [np.mean([r[f"mislabel_{KP}"] for r in need_sim if r["subject"] == s and r["R"] == R]) for s in subs]
    res = [np.mean([r[f"resolved_{KP}"] for r in need_sim if r["subject"] == s and r["R"] == R]) for s in subs]
    fl = sum(r[f"false_{KP}"] for r in need_sim if r["R"] == R)
    out(f"  R={R:2d}: exact naive mislabel {np.median(per):.3f}  (range {min(per):.3f}-{max(per):.3f});"
        f" E-pair-t resolved {np.median(res):.3f}; false certified labels {fl}")

# ---------------------------------------------------------------- use
out("\nUse: E-pair-t resolved share at k'=20 on the real user blocks, t-SNE, medians over subjects")
use = {}
for R in RU + (72,):
    per = []
    for s in subs:
        rs = [r for r in recs if r["subject"] == s and r["R"] == R and tsne(r)]
        per.append(np.mean([r["use"][KP]["faithful"] + r["use"][KP]["distorted"] for r in rs]))
    fa = np.median([np.mean([r["use"][KP]["faithful"] for r in recs if r["subject"] == s and r["R"] == R and tsne(r)]) for s in subs])
    nf = np.median([np.mean([r["naive_faithful"][KP] for r in recs if r["subject"] == s and r["R"] == R and tsne(r)]) for s in subs])
    use[R] = float(np.median(per))
    out(f"  R={R:2d}: resolved {use[R]:.3f} (faithful {fa:.3f}); naive map calls {nf:.3f} faithful")
use_ok = max(use[8], use[16]) >= 0.30
out(f"  use (resolved >= 0.30 with R <= 16): {'PASSES' if use_ok else 'FAILS'}")

# ---------------------------------------------------------------- validity and diagnostics
valid = [r for r in sim if r["kind"] == "valid"]
fw = {s: np.mean([(r["pair_viol"] + r["reach_viol"]) > 0 for r in valid if r["subject"] == s]) for s in subs}
fw_reach = {s: np.mean([r["reach_viol"] > 0 for r in valid if r["subject"] == s]) for s in subs}
worst = max(fw.values())
out(f"\nValidity under resampled real noise (300 draws per subject, R=8): E-pair-t failure rate "
    f"median {np.median(list(fw.values())):.3f}, max {worst:.3f}; reach claims only: max {max(fw_reach.values()):.3f}")
out(f"  validity condition (<= 0.10): {'holds' if worst <= 0.10 else 'FAILS'}")
out("\nDiagnostics (split 0, U+H residuals, 20-dim space), per subject:")
keys = ("mean_residual_variance", "eig_max_over_min", "eig_max_over_mean", "excess_kurtosis_median",
        "excess_kurtosis_max", "s2_cv_over_homoscedastic", "F_session_x_image", "F_session")
for k in keys:
    v = [d[k] for d in diag]
    out(f"  {k:28s} median {np.median(v):7.3f}  range {min(v):7.3f} - {max(v):7.3f}")

if need.startswith("FAILS"):
    rec_ = "stop (need fails at native R)"
elif use_ok:
    rec_ = "proceed to a second dataset"
else:
    rec_ = "need-power tension holds on real data: reframe before any further data"
out(f"\nPre-registered recommendation: {rec_}")
open(f"{OUT}/pilot_analysis.txt", "w").write("\n".join(lines) + "\n")
json.dump(dict(summary=summary, need=need, use=use, use_ok=use_ok, validity_max=worst,
               recommendation=rec_, q_per_subject={str(R): q[R] for R in RU}, q80=q80),
          open(f"{OUT}/pilot_summary.json", "w"), indent=1)
