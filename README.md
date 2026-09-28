# Certified neighborhood reach

Simultaneous, finite-sample confidence intervals on how far each displayed neighborhood of a
2-D display (t-SNE, UMAP, PCA) reaches into the true neighbor ranking, for data measured
with replicates. From the intervals comes a faithful / distorted / unresolved map that is
valid at every tolerance.

**Status: Gates 1 and 2 decided (2026-09-28); Gate 3 needs real replicate data, which needs
downloads.** The protocol is [docs/research_plan.md](docs/research_plan.md), version 1 of
2026-09-28, written before any Gate 2 measurement, with Amendment 1 at its end. The gate
reports are [docs/gate1_prior_art.md](docs/gate1_prior_art.md) and
[docs/gate2_report.md](docs/gate2_report.md).

## Provenance

This project succeeds `setwise-certificates-neighborhood-fidelity` (dropped on 2026-09-23;
evidence base frozen at commit `9bc9f03`). Section 2 of the plan reviews the proposal that
started it. `pilot/2026-09-28/` archives the scratch pilot the proposal rested on, unchanged.
The pilot imports the predecessor repository by absolute path and is not part of the
reproduction. `reach/predecessor/` is the predecessor's code, copied verbatim; see its
`PROVENANCE.md`.

## Layout

| path | contents |
|---|---|
| `docs/research_plan.md` | the protocol: claim, review of the proposal, engines, gate criteria, seeds, amendments |
| `docs/gate1_prior_art.md` | Gate 1: prior-art audit of the combination (four literatures) |
| `docs/gate2_report.md` | Gate 2: validity under seven noise models, robustness, power, need, planning rule |
| `reach/engines.py` | E-ball, E-pair-t, E-pair-pool and the reference E-cmp (plan Section 4) |
| `reach/target.py` | reach, its closed-form box bounds, displays, naive maps (plan Section 3) |
| `reach/noise.py` | the seven noise models (plan Section 5.2) |
| `reach/quantiles.py` | Student and chi-square quantiles verified in arbitrary precision |
| `reach/data.py` | the development geometries: the predecessor's digits panel and synthetic families |
| `reach/predecessor/` | the predecessor's `setwise/` package, verbatim |
| `tests/test_reach.py` | structural tests and the pilot reproduction (plan Section 5.7) |
| `run_gate2_validity.py` | family-wise failure of every engine, 24 configurations, 7900 draws |
| `run_gate2_power.py` | power, need and planning-rule measurements, 20 configurations |
| `analyze_gate2.py` | applies the pre-registered Gate 2 criteria to the saved runs |
| `run_gate2_synthetic_ceiling.py` | same-panel comparison with the predecessor's ball ceiling (Gate 2 report, Section 5) |
| `run_gate2_need_crosstab.py` | naive-map errors against the certificate, point by point (Gate 2 report, Section 6) |
| `results/gate2/` | raw results (JSONL), logs, manifests and the analysis printout |
| `pilot/2026-09-28/` | the scratch pilot and the tail check behind review finding 1 |

## Reproducing

Run from the repository root:

```bash
python3 tests/test_reach.py > results/gate2/test_reach.log
```

```bash
python3 run_gate2_validity.py --workers 30 > results/gate2/validity.log
```

```bash
python3 run_gate2_power.py --workers 32 > results/gate2/power.log && python3 analyze_gate2.py
```

```bash
python3 run_gate2_synthetic_ceiling.py > results/gate2/synthetic_ceiling.txt && python3 run_gate2_need_crosstab.py > results/gate2/need_crosstab.txt
```

On a 36-core workstation the validity runs take about 9 minutes, the power runs under 2,
and each tabulation about a minute. The only data is scikit-learn's bundled digits; nothing
is downloaded. Seeds are fixed per configuration and draw (plan Section 7). The t-SNE
displays of the power runs use `n_jobs=1`; the pilot reproduction uses `n_jobs=8`, as the
pilot did.
