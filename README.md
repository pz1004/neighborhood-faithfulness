# Certified neighborhood reach

Simultaneous, finite-sample confidence intervals on how far each displayed neighborhood of a
2-D display (t-SNE, UMAP, PCA) reaches into the true neighbor ranking, for data measured
with replicates. From the intervals comes a faithful / distorted / unresolved map that is
valid at every tolerance.

**Status (2026-09-29): Gates 1 and 2 passed. The Gate 3 need pilot on THINGS-EEG2 ended
in the protocol's "reframe" branch, and the reframe's prior-art check ended in its
"short empirical note" branch. The note's pre-registered headline held on THINGS-EEG2
and on an independent second dataset; the note is drafted in [note/](note/).**
- Real EEG replicates showed strong need. Naive neighborhood maps mislabel about 20% of
  points, and most of their "faithful" calls are wrong.
- They also showed no use. The certificate resolves no point even at 72 repetitions,
  because of signal-to-noise rather than the error control.
- The reframed "neighborhood noise ceiling" is routine as a method. Only the empirical
  finding is new, and [docs/note_protocol.md](docs/note_protocol.md) tests it with fixed
  claims and thresholds.
- The finding held ([docs/note_confirmatory_results.md](docs/note_confirmatory_results.md)).
  At 8 trials, 81–87% of the faithful calls of t-SNE maps are wrong (conservative bound
  64%). For MDS maps of cross-validated distances it is 62–72% (bound 55%). The headline
  holds in 10 of 12 pre-specified variations. The claim at 32 trials fell just short.
- It was then confirmed independently on a second species and modality
  ([docs/note_allen_results.md](docs/note_allen_results.md)): Allen two-photon, mouse V1,
  10 mice. Essentially every faithful call is wrong at 8 and 16 trials (bound 0.997), in all
  12 variations, and for RSA-style maps too.

Details: [docs/gate3_pilot_things.md](docs/gate3_pilot_things.md) and
[docs/reframe_prior_art.md](docs/reframe_prior_art.md). The protocol is
[docs/research_plan.md](docs/research_plan.md), with Amendments 1 to 3 at its end. The
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
| `docs/reframe_prior_art.md` | prior-art check of the reframed neighborhood noise ceiling: protocol, verdicts, decision |
| `docs/note_protocol.md` | pre-registered confirmatory analysis for the empirical note (plan Amendment 3) |
| `docs/note_confirmatory_results.md` | its results: verdict, primary cell, RSA pipeline, robustness |
| `docs/note_protocol_allen.md` | pre-registered second-dataset protocol: Allen two-photon, mouse V1 (plan Amendment 4) |
| `docs/note_allen_results.md` | its results and the joint conclusion with THINGS-EEG2 |
| `note/` | the empirical note, draft: `note.tex`, with `make_figures.py` and `make_tables.py` drawing every figure, table and quoted number from the committed results; the PDF is built, not tracked |
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
| `run_things_confirmatory.py` | the confirmatory analysis of `docs/note_protocol.md`: 12 cells, both pipelines, estimator check |
| `analyze_things_confirmatory.py` | its estimators, subject bootstrap and pre-registered decision rules |
| `tests/test_things_confirmatory.py` | synthetic known-truth checks of that code; touches no data |
| `fetch_allen_2p.py` | applies the Allen selection rule and fetches recordings, SHA-256 recorded |
| `reach/allen.py` | Allen recordings: per-presentation responses, time strata, blocks, halves |
| `run_allen_confirmatory.py` | the second-dataset analysis of `docs/note_protocol_allen.md` |
| `analyze_allen_confirmatory.py` | its estimators, decision rules and joint conclusion with THINGS-EEG2 |
| `tests/test_allen_confirmatory.py` | synthetic checks with the recordings' file layout; opens no recording |
| `explore_allen_snr.py` | exploratory, post hoc: image structure and neighbor distance against noise, both datasets |
| `data/things_eeg2/` | the downloaded test files (not tracked) and their SHA-256 manifest (tracked) |
| `results/gate2/` | raw results (JSONL), logs, manifests and the analysis printout |
| `results/things_pilot/` | the pilot's records, semi-synthetic runs, diagnostics and printouts |
| `results/things_confirmatory/` | the confirmatory analysis's test log, rows, semi-synthetic check and printout |
| `data/allen_2p/` | Allen metadata and recordings (not tracked) and their SHA-256 manifest (tracked) |
| `results/allen_confirmatory/` | the second-dataset test log, rows, semi-synthetic check and printouts |
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

The confirmatory analysis:

```bash
python3 tests/test_things_confirmatory.py > results/things_confirmatory/test_things_confirmatory.log && python3 run_things_confirmatory.py > results/things_confirmatory/run.log && python3 analyze_things_confirmatory.py
```

The second dataset downloads 4.4 GB (10 Allen recordings, selected by the protocol's rule):

```bash
curl -o data/allen_2p/ophys_experiments.json https://allen-brain-observatory.s3.us-west-2.amazonaws.com/visual-coding-2p/ophys_experiments.json && python3 fetch_allen_2p.py > data/allen_2p/fetch.log
```

```bash
python3 tests/test_allen_confirmatory.py > results/allen_confirmatory/test_allen_confirmatory.log && python3 run_allen_confirmatory.py > results/allen_confirmatory/run.log && python3 analyze_allen_confirmatory.py && python3 explore_allen_snr.py > results/allen_confirmatory/explore_snr.txt
```

The note (needs matplotlib and a LaTeX installation with latexmk):

```bash
python3 note/make_figures.py && python3 note/make_tables.py && latexmk -pdf -cd note/note.tex
```

On a 36-core workstation the validity runs take about 9 minutes, the power runs under 2, and
each tabulation about a minute. Gates 1 and 2 use only scikit-learn's bundled digits. The
THINGS-EEG2 pilot takes about 3 minutes on 10 workers, after the download. Seeds are fixed
per configuration and draw (plan Section 7). The t-SNE displays of the power runs use
`n_jobs=1`; the pilot reproduction uses `n_jobs=8`, as the pilot did.
