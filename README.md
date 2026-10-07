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
reach/          the library: both datasets' loaders, the maps, the six checks, the
                simulation study, and the per-point certificates of the pilot
tests/          synthetic known-truth tests; none opens a recording
docs/           the research plan and each analysis protocol, with its results report
results/        every derived result: rows, logs, run manifests and printouts
data/           SHA-256 manifests and fetch logs of the recordings
fetch_*.py      data acquisition
run_*.py        one driver per analysis; each writes its rows and log to results/
analyze_*.py    the pre-specified estimators and decision rules, applied to those rows
explore_*.py    exploratory, post hoc analyses
pilot/          the archived scratch pilot the research plan started from
```

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
| THINGS-EEG2, preprocessed ([OSF 3JK45](https://osf.io/3jk45/), component [anp5v](https://osf.io/anp5v/)) | the test set of 10 participants: 200 images × 80 repetitions over 4 sessions | 2.1 GB | `fetch_things_eeg2_test.py`, which reads only the test file of each participant's archive |
| Allen Brain Observatory, Visual Coding two-photon, natural scenes | 10 experiments in mouse V1, one per mouse, selected by the rule in `docs/note_protocol_allen.md`: 118 images × 50 repetitions | 4.4 GB | `fetch_allen_2p.py` |

## Reproducing the results

Run the stages in order from the repository root; each depends on the ones before it. They
follow the research plan and its amendments, and every command writes to the path of the
committed output it reproduces. Seeds are fixed per configuration and draw
(`docs/research_plan.md`, Section 7).

**1. Per-point certificates on synthetic data** (research plan, Gates 1 and 2; scikit-learn's
bundled digits; about 15 minutes)

```bash
python3 tests/test_reach.py > results/gate2/test_reach.log
python3 run_gate2_validity.py --workers 30 > results/gate2/validity.log
python3 run_gate2_power.py --workers 32 > results/gate2/power.log
python3 analyze_gate2.py
python3 run_gate2_synthetic_ceiling.py > results/gate2/synthetic_ceiling.txt
python3 run_gate2_need_crosstab.py > results/gate2/need_crosstab.txt
```

**2. Data** (6.5 GB)

```bash
python3 fetch_things_eeg2_test.py > data/fetch.log
curl -o data/allen_2p/ophys_experiments.json https://allen-brain-observatory.s3.us-west-2.amazonaws.com/visual-coding-2p/ophys_experiments.json
python3 fetch_allen_2p.py > data/allen_2p/fetch.log
```

**3. The exploratory EEG pilot** (Amendment 2; about 3 minutes)

```bash
python3 tests/test_things.py > results/things_pilot/test_things.log
python3 run_things_need_pilot.py > results/things_pilot/run.log
python3 analyze_things_need_pilot.py
python3 explore_things_snr.py > results/things_pilot/explore_snr.txt
```

**4. The confirmatory EEG analysis** (Amendment 3, `docs/note_protocol.md`; about 8 minutes)

```bash
python3 tests/test_things_confirmatory.py > results/things_confirmatory/test_things_confirmatory.log
python3 run_things_confirmatory.py > results/things_confirmatory/run.log
python3 analyze_things_confirmatory.py
```

**5. Mouse V1** (Amendment 4, `docs/note_protocol_allen.md`; about 3 minutes)

```bash
python3 tests/test_allen_confirmatory.py > results/allen_confirmatory/test_allen_confirmatory.log
python3 run_allen_confirmatory.py > results/allen_confirmatory/run.log
python3 analyze_allen_confirmatory.py
python3 explore_allen_snr.py > results/allen_confirmatory/explore_snr.txt
python3 describe_allen_subjects.py
```

**6. Common checks against the truth** (Amendment 5, `docs/checks_protocol.md`; about 75
minutes)

```bash
python3 tests/test_checks.py > results/checks/test_checks.log
python3 run_checks.py --dataset things > results/checks/things/run.log
python3 run_checks.py --dataset allen > results/checks/allen/run.log
python3 analyze_checks.py
```

**7. Validating the estimator** (Amendment 6, `docs/revision_protocol.md`; about 2.5 hours).
The development runs fixed the noise grid, the estimator and the thresholds in
`docs/revision_addendum.md` before any validation seed was used.

```bash
python3 analyze_review_audit.py
python3 tests/test_estimand.py > results/revision/test_estimand.log
python3 tests/test_validation.py > results/revision/test_validation.log
python3 run_validation.py --phase dev --alpha
python3 run_validation.py --phase dev > results/revision/dev/run.log
python3 run_validation.py --phase dev --lams 0.25,0.35 --append --cells T1_l0.25,T2_l0.25,T3_l0.25,T1_l0.35,T2_l0.35,T3_l0.35 > results/revision/dev/run_extension.log
python3 analyze_validation.py --phase dev --lams 0.1,0.25,0.35,0.5,1,2 > results/revision/dev/analysis.txt
python3 run_validation.py --phase validation > results/revision/validation/run.log
python3 analyze_validation.py --phase validation > results/revision/validation/analysis.txt
```

**8. Sensitivity analyses and the benchmark** (Amendment 6; about 1 hour). `run_points.py`
re-creates the primary and cv-MDS maps and every check, and compares each count with the
committed rows; its committed logs report 0 counts not reproduced.

```bash
python3 tests/test_stage2.py > results/revision/test_stage2.log
python3 make_trial_manifest.py
python3 run_session_disjoint.py > results/session_disjoint/run.log
python3 analyze_session_disjoint.py
python3 run_points.py --dataset things > results/points/things.log
python3 run_points.py --dataset allen > results/points/allen.log
python3 analyze_benchmark.py
```

**9. Statuses, reporting classes and the study-level calibration** (Amendment 6; about 19
hours, nearly all of it the 1,000 simulated studies of `run_layer_c.sh`)

```bash
python3 analyze_gate_a.py
python3 tests/test_layer_c.py > results/revision/test_layer_c.log
bash run_layer_c.sh
python3 analyze_layer_c_rule.py
```

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
`analyze_review_audit.py` into `results/review_audit/`.

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
