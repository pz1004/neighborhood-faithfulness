# Certified neighborhood reach

Simultaneous, finite-sample confidence intervals on how far each displayed neighborhood of a
2-D display (t-SNE, UMAP, PCA) reaches into the true neighbor ranking, for data measured
with replicates. From the intervals comes a faithful / distorted / unresolved map that is
valid at every tolerance.

**Status (2026-09-29): Gates 1 and 2 passed; the Gate 3 need pilot on THINGS-EEG2 ended
in the protocol's "reframe" branch.** Real EEG replicates showed strong need: naive
neighborhood maps mislabel about 20% of points, and most of their "faithful" calls are
wrong. They also showed no use: the certificate resolves no point even at 72 repetitions,
because of signal-to-noise rather than the error control. The per-point certificate is
therefore informative only where naive maps are already right.
[docs/gate3_pilot_things.md](docs/gate3_pilot_things.md) has the details. The protocol is
[docs/research_plan.md](docs/research_plan.md), with Amendments 1 and 2 at its end. The
earlier reports are [docs/gate1_prior_art.md](docs/gate1_prior_art.md) and
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
| `docs/gate3_pilot_things.md` | Gate 3 need pilot on THINGS-EEG2: need, use, validity under real noise, diagnostics |
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
| `fetch_things_eeg2_test.py` | fetches only the preprocessed test file of each THINGS-EEG2 subject (range requests, CRC-32 checked) |
| `reach/things.py` | THINGS-EEG2 loading (restricted unpickler), session reconstruction, blocks, representation |
| `tests/test_things.py` | loader, session and block checks that touch no data value |
| `run_things_need_pilot.py` | the pilot of plan Amendment 2 |
| `analyze_things_need_pilot.py` | applies Amendment 2's decision rule |
| `explore_things_snr.py` | exploratory, post hoc: why nothing resolves, and false comfort |
| `data/things_eeg2/` | the downloaded test files (not tracked) and their SHA-256 manifest (tracked) |
| `results/gate2/` | raw results (JSONL), logs, manifests and the analysis printout |
| `results/things_pilot/` | the pilot's records, semi-synthetic runs, diagnostics and printouts |
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

The THINGS-EEG2 pilot downloads 2.1 GB (the ten subjects' test files only):

```bash
python3 fetch_things_eeg2_test.py > data/fetch.log && python3 tests/test_things.py > results/things_pilot/test_things.log
```

```bash
python3 run_things_need_pilot.py > results/things_pilot/run.log && python3 analyze_things_need_pilot.py && python3 explore_things_snr.py > results/things_pilot/explore_snr.txt
```

On a 36-core workstation the validity runs take about 9 minutes, the power runs under 2, and
each tabulation about a minute. Gates 1 and 2 use only scikit-learn's bundled digits. The
THINGS-EEG2 pilot takes about 3 minutes on 10 workers, after the download. Seeds are fixed
per configuration and draw (plan Section 7). The t-SNE displays of the power runs use
`n_jobs=1`; the pilot reproduction uses `n_jobs=8`, as the pilot did.
