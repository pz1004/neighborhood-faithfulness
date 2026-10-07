# Neighborhood Faithfulness in Maps of Trial-Averaged Neural Responses

### Estimated from Held-Out Repetitions

**Sooyoung Jang**\
Department of Computer Engineering, Hanbat National University, Daejeon 34158, Republic of Korea\
Correspondence: `syjang@hanbat.ac.kr`

> **Status:** manuscript submitted for publication, not yet peer reviewed. This work was
> supported by the National Research Foundation of Korea (NRF) grant funded by the Korea
> government (MSIT) (RS-2026-25602117); the author declares no competing interests.

## Abstract

Two-dimensional maps of trial-averaged responses, made with t-SNE or multidimensional
scaling, are read neighbor by neighbor and checked against the averaged data (the plug-in
check). Because the map is fitted to the same noisy averages, such checks share its errors.
We define the *false comfort* of a check as the share of the neighborhoods it passes that are
not faithful to the noise-free responses, and estimate it from held-out repetitions,
extrapolating away the reference's own noise. A simulation study defines observable statuses
that say when the estimate can be trusted, and tests whether its lower bound is calibrated.
In human EEG (THINGS-EEG2; 10 participants, 200 images) with 8 trials per image, an
estimated 81% of the t-SNE neighborhoods that passed the plug-in check were not faithful
(calibrated lower bound 64%), against 93% for a random selection. At 16 and 32 trials the
estimate was unresolved; disagreement with a 40-trial held-out reference fell from 88% to
67%. In mouse primary visual cortex (Allen Brain Observatory; 10 mice, 118 images), a
24-trial held-out reference confirmed under 0.1% of the passing neighborhoods. Of six common
checks computed from the same trials, four selected faithful neighborhoods 2.0–4.2 times more
often than a random selection in EEG, and two passed almost nothing. Even the best-scoring 5%
of neighborhoods had an estimated false comfort of 68–73%. Neighborhoods in maps of few-trial
averages should be read against held-out repetitions deep enough for false comfort to be
estimable.

## Repository layout

```
analyses/             one package per stage, in the order of the research plan
  certificates/       per-point certificates on synthetic data (Gates 1 and 2)
  fetch/              fetching the recordings
  eeg_pilot/          the exploratory EEG pilot (Amendment 2)
  eeg_confirmatory/   the confirmatory EEG analysis (Amendment 3)
  mouse_v1/           mouse V1 (Amendment 4)
  checks/             common checks against the truth (Amendment 5)
  validation/         the estimator's simulation study (Amendment 6)
  benchmark/          sensitivity analyses and the benchmark of the checks (Amendment 6)
  calibration/        statuses, reporting classes and Layer C (Amendment 6)
scripts/              one runner per stage, in order
reach/                the library: both datasets' loaders, the maps, the six checks, the
                      simulation study, and the per-point certificates of the pilot
tests/                synthetic known-truth tests; none opens a recording
docs/                 the research plan and each analysis protocol, with its results report
results/              every derived result: rows, logs, run manifests and printouts
data/                 SHA-256 manifests and fetch logs of the recordings
pilot/                the archived scratch pilot the research plan started from
```

Within a stage, `run_*.py` drives an analysis and writes its rows and log to `results/`,
`analyze_*.py` applies the pre-specified estimators and decision rules to those rows, and
`explore_*.py` is exploratory and post hoc.

`reach/predecessor/` is copied verbatim from the author's earlier project (see its
`PROVENANCE.md`). The archived pilot imports that project by absolute path and is not part
of the reproduction; `tests/test_reach.py` checks the copied code against the pilot's saved
output.

Not in the repository: the THINGS-EEG2 and Allen Brain Observatory recordings, which their
providers distribute and stage 2 below fetches, and the manuscript. Every number the
manuscript reports is computed from `results/`; see [Results](#results).

## Installation

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
```

Results were produced with Python 3.13.5, NumPy 2.1.3, SciPy 1.15.3 and scikit-learn 1.6.1,
as each run's `manifest.json` records, with h5py 3.12.1 and mpmath 1.3.0, on a 36-thread
Intel Core i9-10980XE. The t-SNE maps depend on the scikit-learn version. The run times
below are for that machine.

## Data

Two public datasets. Each fetcher records the SHA-256 digest of every file it writes in
`data/*/manifest.json`, and the confirmatory runs refuse a file whose digest differs.

| Dataset | Used here | Download | Fetched by |
|---|---|---:|---|
| THINGS-EEG2, preprocessed ([OSF 3JK45](https://osf.io/3jk45/), component [anp5v](https://osf.io/anp5v/)) | the test set of 10 participants: 200 images × 80 repetitions over 4 sessions | 2.1 GB | `analyses/fetch/fetch_things_eeg2_test.py`, which reads only the test file of each participant's archive |
| Allen Brain Observatory, Visual Coding two-photon, natural scenes | 10 experiments in mouse V1, one per mouse, selected by the rule in `docs/note_protocol_allen.md`: 118 images × 50 repetitions | 4.4 GB | `analyses/fetch/fetch_allen_2p.py` |

## Reproducing the results

Run the stages in order; each depends on the ones before it. Each runner can be started from
any directory, runs that stage's tests and scripts from the repository root, and writes to the
paths of the committed outputs it reproduces. Seeds are fixed per configuration and draw
(`docs/research_plan.md`, Section 7).

```bash
./scripts/01_certificates.sh       # Gates 1-2: per-point certificates, synthetic data    ~15 min
./scripts/02_fetch.sh              # the recordings, 6.5 GB
./scripts/03_eeg_pilot.sh          # Amendment 2: the exploratory EEG pilot              ~3 min
./scripts/04_eeg_confirmatory.sh   # Amendment 3: the confirmatory EEG analysis          ~8 min
./scripts/05_mouse_v1.sh           # Amendment 4: mouse V1                               ~3 min
./scripts/06_checks.sh             # Amendment 5: common checks against the truth        ~75 min
./scripts/07_validation.sh         # Amendment 6: the estimator's simulation study       ~2.5 h
./scripts/08_benchmark.sh          # Amendment 6: sensitivity analyses and the benchmark ~1 h
./scripts/09_calibration.sh        # Amendment 6: statuses, classes and Layer C          ~19 h
```

**Stage 9 takes most of a day.** Nearly all of it is Layer C, 1,000 simulated studies in each of
six scenarios (`analyses/calibration/run_layer_c.sh`); the runner stops if any of them fails.

`analyses/benchmark/run_points.py` re-creates the primary and cv-MDS maps and every check of
stages 4 to 6, and reports how many committed counts it did not reproduce; the committed logs
report 0. Stage 1's validity run resumes from `results/gate2/validity.jsonl`; delete that file
to recompute every draw. Stage 7 writes the development decisions to
`results/revision/dev/decisions_rerun.json` and leaves the committed ones alone.

Any script can also be run on its own from the repository root, for example
`python3 analyses/checks/run_checks.py --dataset things`.

| Stage | Folder | Scripts (`.py`) | Writes to `results/` |
|---|---|---|---|
| 1 | `analyses/certificates/` | `run_gate2_validity`, `run_gate2_power`, `run_gate2_synthetic_ceiling`, `run_gate2_need_crosstab`, `analyze_gate2` | `gate2/` |
| 2 | `analyses/fetch/` | `fetch_things_eeg2_test`, `fetch_allen_2p` | (`data/`) |
| 3 | `analyses/eeg_pilot/` | `run_things_need_pilot`, `analyze_things_need_pilot`, `explore_things_snr` | `things_pilot/` |
| 4 | `analyses/eeg_confirmatory/` | `run_things_confirmatory`, `analyze_things_confirmatory` | `things_confirmatory/` |
| 5 | `analyses/mouse_v1/` | `run_allen_confirmatory`, `analyze_allen_confirmatory`, `explore_allen_snr`, `describe_allen_subjects` | `allen_confirmatory/` |
| 6 | `analyses/checks/` | `run_checks`, `analyze_checks` | `checks/` |
| 7 | `analyses/validation/` | `analyze_review_audit`, `run_validation`, `analyze_validation` | `review_audit/`, `revision/` |
| 8 | `analyses/benchmark/` | `make_trial_manifest`, `run_session_disjoint`, `analyze_session_disjoint`, `run_points`, `analyze_benchmark` | `manifest/`, `session_disjoint/`, `points/`, `benchmark/` |
| 9 | `analyses/calibration/` | `analyze_gate_a`, `run_study_calibration`, `analyze_study_calibration`, `run_layer_c.sh`, `analyze_layer_c_rule` | `gate_a/`, `revision/layer_c/` |

The protocols in `docs/` name the scripts by their earlier paths in the repository root. They
moved into `analyses/` in commit `5198734`, after every committed result had been produced,
with their names unchanged; `git log --follow analyses/<stage>/<script>` traces each one.

## Results

| | |
|---|---|
| EEG, 8 trials per image — *calibrated* | An estimated 81% of the t-SNE neighborhoods that pass the plug-in check are not faithful (lower bound 64%), against 93% for a random selection of as many neighborhoods. In 1,000 simulated studies the bound exceeded the true value in 2.0% and 1.6% of them (± 0.4%) under two truths, against a nominal 2.5% |
| EEG, 16 and 32 trials — *finite-reference* | Unresolved; disagreement with a 40-trial held-out reference falls from 88% to 67% |
| Mouse V1, 8 trials — *finite-reference* | A 24-trial held-out reference confirmed 5 of 5,919 passing neighborhoods (under 0.1%) |
| Six common checks, EEG — *model-dependent* | Four select faithful neighborhoods 2.0–4.2 times more often than a random selection, and bootstrap and split-half recurrence pass almost nothing; even the best-scoring 5% have an estimated false comfort of 68–73% |
| MDS of cross-validated distances, EEG, 8 trials — *model-dependent* | An estimated 62–72% of the passing neighborhoods are not faithful |

Each result is reported in the class its observable status allows (`docs/revision_protocol.md`,
Section 5, as applied in `results/gate_a/classification.json`). Only a *calibrated* result
carries a threshold claim; a *model-dependent* one is an estimate without one; a
*finite-reference* one is only the disagreement with the deepest held-out reference. The
numbers are in `results/things_confirmatory/confirmatory_summary.json`,
`results/benchmark/summary.json`, `results/gate_a/classification.json` and
`results/revision/layer_c/`.

## Protocols and the commit trail

Each protocol and the code implementing it were committed before the analysis it governs was
run. The EEG recordings had first been analyzed in the exploratory pilot (Amendment 2); the
confirmatory analysis used new splits and seeds, but not new data. Amendment 6 was written
after the first results and an internal review; its simulation study, statuses and reporting
rules were committed before they were run. `git show <commit>` shows each step:

| Commit | Date | What it fixed or recorded |
|---|---|---|
| `ce69f4e` | 2026-09-28 | Research plan, before any analysis |
| `5024fe6` | 2026-09-29 | Amendment 2: the exploratory THINGS-EEG2 pilot, before its data were loaded |
| `c7e5a93` | 2026-09-29 | Pilot results |
| `598562e` | 2026-09-29 | Amendment 3: the confirmatory EEG protocol, its code and synthetic tests |
| `5c2bd7f` | 2026-09-29 | Confirmatory EEG results |
| `985098d` | 2026-09-29 | Amendment 4: the mouse V1 protocol, before any recording was downloaded |
| `ea36f1d` | 2026-09-29 | Mouse V1 code and synthetic tests, before any response value was read |
| `0782aa2` | 2026-09-29 | Mouse V1 recordings downloaded, with SHA-256 digests |
| `a2f39ca` | 2026-09-29 | Mouse V1 results |
| `766e4b7` | 2026-10-01 | Amendment 5: the comparison of common checks, its code and synthetic tests |
| `3e3204e` | 2026-10-01 | Results of the comparison |
| `c5df05b` | 2026-10-04 | Amendment 6, written after the first results and an internal review: the simulation study, statuses, reporting classes and sensitivity analyses |
| `817f138` | 2026-10-04 | Simulation code and synthetic tests, before the development pilot |
| `cfa8ba2` | 2026-10-04 | The development pilot's decisions, before any validation seed |
| `200cf63` | 2026-10-04 | Code for the classes and the study-level calibration, and the trial manifest, before they were used |
| `5e16193` | 2026-10-04 | The mislabel rate's status from its own curve, before the validation was analyzed |
| `a3c282d` | 2026-10-04 | Validation results and the class of every real result |
| `a1047b0` | 2026-10-04 | Results of the sensitivity analyses and the benchmark |
| `1587f15` | 2026-10-04 | How the study-level calibration is judged, before its output was read |
| `62df069` | 2026-10-05 | Study-level calibration results |
| `262d8c3` | 2026-10-05 | The author's choice of how the results are reported |

The research plan and its six amendments are in `docs/research_plan.md`; the protocols, the
addendum and the results reports are beside it in `docs/`. Internal working documents that
some of them mention (prior-art audits, a gate report, the internal review and its audit) are
not part of this release; the evidence the audit computed is reproduced by
`analyses/validation/analyze_review_audit.py` into `results/review_audit/`.

## Citation

```bibtex
@unpublished{jang2026neighborhood,
  title  = {Neighborhood faithfulness in maps of trial-averaged neural responses,
            estimated from held-out repetitions},
  author = {Jang, Sooyoung},
  year   = {2026},
  note   = {Manuscript submitted for publication}
}
```

## Acknowledgement

The EEG recordings are THINGS-EEG2, and the mouse recordings come from the Allen Brain
Observatory. Please cite both datasets and their papers when you use the derived results:

- A. T. Gifford, K. Dwivedi, G. Roig and R. M. Cichy, "A large and rich EEG dataset for
  modeling human visual object recognition," *NeuroImage* 264 (2022) 119754,
  [doi:10.1016/j.neuroimage.2022.119754](https://doi.org/10.1016/j.neuroimage.2022.119754);
  data on OSF, [doi:10.17605/OSF.IO/3JK45](https://doi.org/10.17605/OSF.IO/3JK45).
- S. E. J. de Vries, J. A. Lecoq, M. A. Buice et al., "A large-scale standardized
  physiological survey reveals functional organization of the mouse visual cortex,"
  *Nature Neuroscience* 23 (2020) 138–151,
  [doi:10.1038/s41593-019-0550-9](https://doi.org/10.1038/s41593-019-0550-9); data from the
  Allen Institute MindScope Program, *Allen Brain Observatory – 2-photon Visual Coding*
  (2016), <https://observatory.brain-map.org/visualcoding>.

The scDEED check in `reach/checks.py` reimplements L. Xia, C. Lee and J. J. Li, "Statistical
method scDEED for detecting dubious 2D single-cell embeddings and optimizing t-SNE and UMAP
hyperparameters," *Nature Communications* 15 (2024) 1753.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

The THINGS-EEG2 and Allen Brain Observatory data are not covered by it and are not
redistributed here; they carry their providers' terms.
