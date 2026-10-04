#!/usr/bin/env python3
"""Analysis of the simulation study (docs/revision_protocol.md, Sections 3.5 to 4).

    python3 analyze_validation.py --phase dev           # the development pilot and its decisions
    python3 analyze_validation.py --phase validation    # performance, by status, and Layer A

Each repetition's pooled counts (run_validation.py) give its disagreement curves, the estimates
of every method and the observables of Section 4. A cell's status is computed from the medians
over its repetitions, as a real cell's is from the medians over its subjects.

The development phase applies the rules of Section 3.7, in this order:
1. Saturation, under the default thresholds, decides which cells are unsaturated.
2. The candidate estimator is chosen there by the selection rule of Section 3.5.
3. The status thresholds are tuned for that estimator.

It writes results/revision/dev/decisions.json, which the validation run reads once it is
committed with the addendum.
"""
import argparse, collections, functools, gzip, itertools, json, os
import numpy as np

import analyze_allen_confirmatory as AA
from reach import sim

OUT = "results/revision"
G = {"sqrt": lambda h: h ** -0.5, "inv": lambda h: 1.0 / h}
DEFAULT = dict(A_min=20, c_sat=0.95, f_sat=0.02, c_fd=0.03, c_div=0.10)
GRID = dict(A_min=(10, 20, 50), c_sat=(0.90, 0.95, 0.98), f_sat=(0.01, 0.02, 0.05),
            c_fd=(0.02, 0.03, 0.05), c_div=(0.05, 0.10, 0.15))
CANDIDATES = ("a", "b", "c")          # Section 3.5: R_h^-1/2 alone; more reference sizes; cv references
ESTIMATORS = ("sqrt", "inv", "smaller", "larger") + CANDIDATES
FAMILY_OF = {"sqrt": "plug", "inv": "plug", "smaller": "plug", "larger": "plug", "a": "plug", "b": "plug", "c": "cv"}
BIAS_OK, RMSE_OK = 0.03, 0.08         # estimable cells must meet these (Sections 3.5, 4)
STATUSES = ("estimable", "extrapolation unresolved", "reference saturated", "insufficient acceptance")
DIAGNOSTIC = ("_D", "_S")             # dependent-noise and stress cells: diagnostics, never matching cells
# How Section 5 is applied, fixed by the addendum before any validation seed is used.
GATE_A = dict(
    accuracy_of="the smaller of the two published extrapolants, on which the bound L rests; the published range "
                "[smaller, larger] is reported with its bracket coverage",
    matching="validation cells of the same dataset and map type, with the same status, whose median d(R_h,max) and "
             "fall are both within 0.05 of the real cell's; the dependent-noise (D) and stress (S) cells are "
             "diagnostics and never match",
    criterion="every matching cell must meet the accuracy criterion (the worst matching cell decides)",
    readings="two classifications are reported: strict (Section 3.5: no candidate qualified, so every latent "
             "estimate is model-dependent) and local (Section 5, cell by cell); the author chooses at Gate A",
    layer_c="under the local reading, for real cells classified estimable with their matching cells accurate; the "
            "guard in the first 200 studies only, as the timing check projects over 6 h per scenario with it",
)


@functools.lru_cache(maxsize=None)
def load(phase, ds):
    """{(cell, base, rep, layer, regen, R): {(ref, Rh): counts}} from a run's rows."""
    reps = collections.defaultdict(dict)
    for line in gzip.open(f"{OUT}/{phase}/{ds}.jsonl.gz", "rt"):
        r = json.loads(line)
        key = (r["cell"], r["base"], r["rep"], r["layer"], r["regen"], r["R"])
        reps[key][(r["ref"], r["Rh"])] = np.array([r[f] for f in sim.FIELDS], float)
    return reps


def fc(c, ref, h):
    """Disagreement d_FC = n10 / (n1 * nb) against reference (ref, h); NaN where nothing is called faithful."""
    n1nb, n10 = c[(ref, h)][0], c[(ref, h)][1]
    return n10 / n1nb if n1nb else np.nan


def intercept(sizes, y, g):
    if not np.isfinite(y).all():
        return np.nan
    return AA.intercept(np.array([G[g](h) for h in sizes]), np.asarray(y, float))


def record(c, L, exact_key=("exact", 0)):
    """Estimates and observables of one repetition (or one real subject) from its pooled counts."""
    RH, ex = L["RH"], L["extra"]
    plug = [fc(c, "plug", h) for h in RH]
    rec = dict(exact=fc(c, *exact_key) if exact_key in c else np.nan,
               acc=float(c[("plug", RH[0])][5]) if ("plug", RH[0]) in c else np.nan,
               sqrt=intercept(RH, plug, "sqrt"), inv=intercept(RH, plug, "inv"))
    rec["smaller"], rec["larger"] = min(rec["sqrt"], rec["inv"]), max(rec["sqrt"], rec["inv"])
    rec["a"] = rec["sqrt"]
    if ("plug", ex) in c:
        rec["b"] = intercept(RH + (ex,), plug + [fc(c, "plug", ex)], "sqrt")
    if ("cv", RH[0]) in c:
        cv = [fc(c, "cv", h) for h in RH]
        rec["c"] = intercept(RH, cv, "sqrt")
        rec["cv_inv"] = intercept(RH, cv, "inv")
        rec["obs_cv"] = observables(cv, RH, rec["c"], rec["cv_inv"])
    rec["obs_plug"] = observables(plug, RH, rec["sqrt"], rec["inv"])
    rec["curve"] = {f"plug|{h}": fc(c, "plug", h) for h in RH + (ex,) if ("plug", h) in c}
    rec["curve"].update({f"deep|{h}": fc(c, "deep", h) for h in L["deep"] if ("deep", h) in c})
    return rec


def observables(d, RH, e_sqrt, e_inv):
    """Section 4's observables from a disagreement curve at the published reference sizes."""
    x = [G["sqrt"](h) for h in RH]
    pred = d[0] + (d[1] - d[0]) / (x[1] - x[0]) * (x[2] - x[0])
    return dict(dmax=d[-1], fall=d[0] - d[-1], fd=pred - d[-1], div=abs(e_sqrt - e_inv))


def status(recs, th, family="plug"):
    """The status of a cell (Section 4) from its subjects' or repetitions' records."""
    acc = np.array([r["acc"] for r in recs], float)
    if (acc >= th["A_min"]).sum() < len(recs) / 2:
        return "insufficient acceptance"
    obs = [r[f"obs_{family}"] for r in recs]
    med = lambda k: float(np.nanmedian([o[k] for o in obs]))
    if med("dmax") >= th["c_sat"] and med("fall") <= th["f_sat"]:
        return "reference saturated"
    if float(np.nanmedian([abs(o["fd"]) for o in obs])) > th["c_fd"] or med("div") > th["c_div"]:
        return "extrapolation unresolved"
    return "estimable"


def performance(recs, est, rng=None):
    """Section 3.6's point-estimate measures of one method over repetitions, with Monte Carlo errors."""
    pairs = [(r[est], r["exact"]) for r in recs if est in r and np.isfinite(r[est]) and np.isfinite(r["exact"])]
    if len(pairs) < 2:
        return dict(n=len(pairs))
    e = np.array([a - b for a, b in pairs])
    est_v = np.array([a for a, _ in pairs])
    n = len(e)
    mse = float((e ** 2).mean())
    rng = rng or np.random.default_rng(0)
    boot = np.median(e[rng.integers(0, n, size=(1000, n))], axis=1)
    return dict(n=n, median_bias=float(np.median(e)), mcse_median=float(boot.std(ddof=1)),
                mean_bias=float(e.mean()), mcse_mean=float(e.std(ddof=1) / np.sqrt(n)),
                mae=float(np.abs(e).mean()), rmse=float(np.sqrt(mse)),
                mcse_rmse=float((e ** 2).std(ddof=1) / (2 * np.sqrt(mse) * np.sqrt(n))) if mse > 0 else 0.0,
                q90=float(np.percentile(np.abs(e), 90)),
                out_of_range=float(((est_v < 0) | (est_v > 1)).mean()))


def layer_b(phase, ds):
    """{(cell, R): [record per repetition]} for Layer B."""
    L = sim.LAYOUT[ds]
    out = collections.defaultdict(list)
    for (cell, b, rep, layer, regen, R), c in sorted(load(phase, ds).items()):
        if layer == "B":
            out[(cell, R)].append(dict(record(c, L), base=b, rep=rep))
    return out


def layer_a(phase, ds):
    """{cell: [per repetition: mean error of each estimator over regenerations, and the mean curve]}."""
    L = sim.LAYOUT[ds]
    by = collections.defaultdict(dict)
    for (cell, b, rep, layer, regen, R), c in load(phase, ds).items():
        if layer == "A":
            by[(cell, b, rep, R)][regen] = c
    out = collections.defaultdict(list)
    for (cell, b, rep, R), regs in sorted(by.items()):
        exact = fc(regs[-1], "exact", 0)
        recs = [record({**c, ("exact", 0): regs[-1][("exact", 0)]}, L) for g, c in regs.items() if g >= 0]
        sizes = L["RH"] + L["deep"]
        curve = np.nanmean([[r["curve"][f"plug|{h}"] if h in L["RH"] else r["curve"][f"deep|{h}"] for h in sizes]
                            for r in recs], axis=0)
        out[cell].append(dict(exact=exact, R=R,
                              err={e: float(np.nanmean([r[e] - exact for r in recs if e in r])) for e in ESTIMATORS},
                              curve=dict(zip(map(str, sizes), curve.tolist())),
                              deep_sqrt=intercept(sizes, curve, "sqrt") - exact,
                              grid_sqrt=intercept(L["RH"], curve[:len(L["RH"])], "sqrt") - exact))
    return out


def bracket(recs):
    """Share of repetitions whose exact value lies within the published range [smaller, larger]."""
    ok = [r for r in recs if np.isfinite(r["smaller"]) and np.isfinite(r["exact"])]
    return float(np.mean([r["smaller"] <= r["exact"] <= r["larger"] for r in ok])) if ok else np.nan


def summarize_cells(B, th, est_family):
    rows = {}
    for (cell, R), recs in sorted(B.items()):
        rows[f"{cell}|{R}"] = dict(bracket=bracket(recs), diagnostic=any(t in cell for t in DIAGNOSTIC),
            status=status(recs, th, est_family), exact=float(np.nanmedian([r["exact"] for r in recs])),
            obs={k: float(np.nanmedian([r["obs_plug"][k] for r in recs])) for k in ("dmax", "fall", "fd", "div")},
            perf={e: performance(recs, e) for e in ESTIMATORS})
    return rows


def choose_estimator(B, narrow=False):
    """Section 3.5: among candidates with |median bias| <= 0.03 in every unsaturated cell (under the
    default thresholds, excluding cells with too few faithful calls), the smallest RMSE; ties go to (a).
    narrow=True is a sensitivity, not the rule: screening cells with default status 'estimable' only."""
    usable = {k: v for k, v in B.items() if status(v, DEFAULT) in ("estimable", "extrapolation unresolved")}
    if narrow:
        usable = {k: v for k, v in usable.items() if status(v, DEFAULT) == "estimable"
                  and k[0].split(":")[1].split("_")[1].startswith("l")}
    table = {}
    for e in CANDIDATES:
        biases = {f"{k[0]}|{k[1]}": performance(v, e).get("median_bias", np.nan) for k, v in usable.items()}
        errs = [r[e] - r["exact"] for v in usable.values() for r in v
                if e in r and np.isfinite(r[e]) and np.isfinite(r["exact"])]
        table[e] = dict(qualifies=bool(biases) and all(abs(b) <= BIAS_OK for b in biases.values()),
                        rmse=float(np.sqrt(np.mean(np.square(errs)))) if errs else np.nan, median_bias=biases)
    ok = [e for e in CANDIDATES if table[e]["qualifies"]]
    if not ok:
        return "published", table, sorted(f"{k[0]}|{k[1]}" for k in usable)
    best = min(ok, key=lambda e: (round(table[e]["rmse"], 3), CANDIDATES.index(e)))
    return best, table, sorted(f"{k[0]}|{k[1]}" for k in usable)


def tune_thresholds(B, est):
    """Section 4: the most permissive thresholds for which every estimable development cell meets
    |median bias| <= 0.03 and RMSE <= 0.08 with the chosen estimator; ties go to the defaults."""
    fam = FAMILY_OF.get(est, "plug")
    perf = {k: performance(v, est if est != "published" else "smaller") for k, v in B.items()}
    best, best_key = None, None
    for vals in itertools.product(*GRID.values()):
        th = dict(zip(GRID, vals))
        est_cells = [k for k, v in B.items() if status(v, th, fam) == "estimable"]
        if not all(abs(perf[k].get("median_bias", np.inf)) <= BIAS_OK and perf[k].get("rmse", np.inf) <= RMSE_OK
                   for k in est_cells):
            continue
        key = (-len(est_cells), sum(th[p] != DEFAULT[p] for p in th))
        if best_key is None or key < best_key:
            best, best_key = th, key
    return best, (-best_key[0] if best_key else 0)


def coverage(rows):
    """Section 3.7: does the exact false comfort across cells span the four ranges, with cells near
    0.25 and 0.5?"""
    ex = [v["exact"] for v in rows.values() if np.isfinite(v["exact"])]
    bins = dict(below_025=sum(x < 0.25 for x in ex), from_025_to_05=sum(0.25 <= x < 0.5 for x in ex),
                from_05_to_09=sum(0.5 <= x < 0.9 for x in ex), above_09=sum(x >= 0.9 for x in ex),
                near_025=sum(abs(x - 0.25) <= 0.1 for x in ex), near_05=sum(abs(x - 0.5) <= 0.1 for x in ex))
    return bins, all(v > 0 for v in bins.values())


def print_cells(name, rows, est_show=("sqrt", "inv", "a", "b", "c")):
    print(f"\n{name}: per cell and budget, medians over repetitions; bias = median of estimate minus exact")
    print(f"  {'cell':12s} {'R':>3s}  {'exact':>6s}  {'d_max':>6s} {'fall':>6s} {'fd':>7s} {'div':>6s}  status"
          "                     " + "  ".join(f"{e:>14s}" for e in est_show) + "  bracket")
    for key, v in rows.items():
        cell, R = key.split("|")
        o = v["obs"]
        perf = "  ".join((f"{v['perf'][e]['median_bias']:+.3f} ({v['perf'][e]['rmse']:.3f})"
                          if v["perf"][e].get("n", 0) >= 2 else f"{'n/a':>14s}") for e in est_show)
        print(f"  {cell:12s} {R:>3s}  {v['exact']:6.3f}  {o['dmax']:6.3f} {o['fall']:6.3f} {o['fd']:+7.3f} {o['div']:6.3f}"
              f"  {v['status']:26s} {perf}  {v['bracket']:5.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", choices=("dev", "validation"), required=True)
    ap.add_argument("--lams", default="0.1,0.5,1,2", help="development: the noise grid the validation will use")
    a = ap.parse_args()
    out = {}
    if a.phase == "dev":
        for ds in ("EEG", "mouse"):
            B = layer_b("dev", ds)
            rows = summarize_cells(B, DEFAULT, "plug")
            print_cells(f"{ds} (development, default thresholds)", rows)
            cov, ok = coverage(rows)
            print(f"  exact false comfort over cells: {cov} -> {'spans the ranges' if ok else 'DOES NOT span the ranges'}")
            out[ds] = dict(B=B, rows=rows, coverage=cov, coverage_ok=ok)
        # the estimator and thresholds are chosen on both datasets' development cells together
        Ball = {(f"{ds}:{k[0]}", k[1]): v for ds in out for k, v in out[ds]["B"].items()}
        est, table, usable = choose_estimator(Ball)
        print(f"\nEstimator selection over {len(usable)} usable development cells:")
        for e, t in table.items():
            worst = max((abs(b) for b in t["median_bias"].values()), default=np.nan)
            print(f"  ({e}) qualifies {t['qualifies']}, pooled RMSE {t['rmse']:.3f}, largest |median bias| {worst:.3f}")
        print(f"  -> chosen: {est}")
        est_n, table_n, usable_n = choose_estimator(Ball, narrow=True)
        print(f"Sensitivity, not the rule (screening cells estimable under the defaults, {len(usable_n)} cells): "
              + ", ".join(f"({e}) {'qualifies' if t['qualifies'] else 'fails'}" for e, t in table_n.items())
              + f" -> {est_n}")
        th, n_est = tune_thresholds(Ball, est)
        tuned = th is not None
        if not tuned:
            th = dict(DEFAULT)
            n_est = sum(status(v, th, FAMILY_OF.get(est, "plug")) == "estimable" for v in Ball.values())
        print(f"Thresholds: {th} ({'tuned' if tuned else 'no set met the criterion, so the defaults stand'};"
              f" {n_est} estimable development cells)")
        al = json.load(open(f"{OUT}/dev/alpha.json"))
        dec = dict(r=al["r"], alpha={ds: al[ds]["alpha"] for ds in ("EEG", "mouse")},
                   lams=[float(x) for x in a.lams.split(",")],
                   estimator=est, thresholds=th, thresholds_tuned=tuned, n_estimable_dev_cells=n_est, gate_a=GATE_A,
                   sensitivity_narrow=dict(estimator=est_n, usable_cells=usable_n,
                                           qualifies={e: t["qualifies"] for e, t in table_n.items()}),
                   selection={e: {k: v for k, v in t.items()} for e, t in table.items()}, usable_cells=usable,
                   coverage={ds: dict(bins=out[ds]["coverage"], ok=out[ds]["coverage_ok"]) for ds in out},
                   cells={ds: out[ds]["rows"] for ds in out})
        json.dump(dec, open(f"{OUT}/dev/decisions.json", "w"), indent=1, default=float)
        print(f"-> {OUT}/dev/decisions.json")
        return
    dec = json.load(open(f"{OUT}/dev/decisions.json"))
    th, est = dec["thresholds"], dec["estimator"]
    fam = FAMILY_OF.get(est, "plug")
    summary = dict(estimator=est, thresholds=th)
    for ds in ("EEG", "mouse"):
        B = layer_b("validation", ds)
        rows = summarize_cells(B, th, fam)
        print_cells(f"{ds} (validation, frozen thresholds)", rows)
        A = layer_a("validation", ds)
        la = {}
        print(f"\n{ds} Layer A (first split's maps fixed; references regenerated): mean error over"
              " regenerations, medians over repetitions")
        for cell, recs in A.items():
            la[cell] = dict(exact=float(np.nanmedian([r["exact"] for r in recs])),
                            err={e: float(np.nanmedian([r["err"][e] for r in recs])) for e in ESTIMATORS},
                            deep_sqrt=float(np.nanmedian([r["deep_sqrt"] for r in recs])),
                            grid_sqrt=float(np.nanmedian([r["grid_sqrt"] for r in recs])),
                            curve={h: float(np.nanmedian([r["curve"][h] for r in recs])) for h in recs[0]["curve"]})
            v = la[cell]
            print(f"  {cell:12s} exact {v['exact']:.3f}  sqrt {v['err']['sqrt']:+.3f}  inv {v['err']['inv']:+.3f}"
                  f"  (b) {v['err']['b']:+.3f}  (c) {v['err']['c']:+.3f}   sqrt over the real grid {v['grid_sqrt']:+.3f},"
                  f" with deep sizes {v['deep_sqrt']:+.3f}   curve " + " ".join(f"{h}:{x:.3f}" for h, x in v["curve"].items()))
        by_status = {}
        for s in STATUSES:
            recs = [r for (cell, R), v in B.items() for r in v if rows[f"{cell}|{R}"]["status"] == s]
            by_status[s] = {e: performance(recs, e) for e in ESTIMATORS}
            by_status[s]["bracket"] = bracket(recs)
        summary[ds] = dict(cells=rows, layer_a=la, by_status=by_status)
    json.dump(summary, open(f"{OUT}/validation/summary.json", "w"), indent=1, default=float)
    print(f"-> {OUT}/validation/summary.json")


if __name__ == "__main__":
    main()
