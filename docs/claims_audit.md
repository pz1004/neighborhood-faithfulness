# Claims audit for Gate B (plan Amendment 6, Stage 3)

**What this is.** Every sentence of the abstract, the highlights and the conclusions of the rewritten manuscript (`manuscript/manuscript.tex`, `manuscript/highlights.txt`), mapped to:
- the committed result it rests on;
- that result's status and reporting class, under the reading the author chose at Gate A (`docs/revision_results.md`, Section 7.3);
- the table or figure that shows it.

Every number below is written into the `.tex` files by `manuscript/make_tables.py` from the files named; `python3 manuscript/make_tables.py --check` verifies them. The highlights are plain text and are checked here by hand.

**The classes** (manuscript Section 2.3.4):
- **calibrated:** threshold claims allowed;
- **model-dependent:** an estimate, without a threshold claim;
- **finite-reference:** only the disagreement with the deepest reference.

## 1. Abstract

| # | Sentence (abridged) | Rests on | Status and class | Shown in |
|---|---|---|---|---|
| A1 | Maps of trial-averaged responses are read neighbor by neighbor and checked against the averaged data (the plug-in check). | Background; Introduction and its citations | — | — |
| A2 | Because the map is fitted to the same noisy averages, such checks share its errors. | The argument of Section 1; Supplementary Section S1.3 | Argument, not a measurement | — |
| A3 | False comfort is defined, and estimated from held-out repetitions with an extrapolation for the reference's noise. | Methods 2.1–2.2 | Method | Figure 1 |
| A4 | A simulation study defines observable statuses and tests whether the lower bound is calibrated. | `results/revision/validation/summary.json`; `results/revision/layer_c/calibration_*.json` | Method | Table 2, Figure 2 |
| A5 | EEG at 8 trials: an estimated 81% of the passing t-SNE neighborhoods were not faithful (calibrated lower bound 64%), against 93% for a random selection. | 81 and 64: `results/things_confirmatory/confirmatory_summary.json` (`primary\|8\|FC`, smaller extrapolant 0.811, L 0.642). Status estimable with 3 matching cells (worst −0.026 / 0.049): `results/gate_a/classification.json`. Calibrated under T1 (0.020 ± 0.004) and T3 (0.016 ± 0.004): `results/revision/layer_c/`. 93: `results/benchmark/summary.json` (`N\|8`, random 0.934) | **calibrated**: the only threshold claim in the paper | Table 3, Table 4, Figure 2c, Figure 3a |
| A6 | At 16 and 32 trials the estimate was unresolved; disagreement with a 40-trial held-out reference fell from 88% to 67%. | Statuses `primary\|16\|FC` and `primary\|32\|FC`, unresolved; finite-reference medians 0.877 and 0.673: `classification.json` | **finite-reference** | Table 4, Figure 3a |
| A7 | Mouse V1: a 24-trial held-out reference confirmed under 0.1% of the passing neighborhoods. | 5 of 5,919 calls at 8 trials: `classification.json` (`primary\|8\|FC`, `pooled`), a count added after the classification and marked so in the text | **finite-reference** (saturated) | Table 4, Figure 3b |
| A8 | Of six common checks, four selected faithful neighborhoods 2.0–4.2 times more often than a random selection in EEG, and two passed almost nothing. | Enrichment: TW 3.06, SCD 2.16, CV 4.17, SEED 2.01. Passes: SH 1.2%, BOOT 0.0%. Both from `results/benchmark/summary.json` | **model-dependent**: estimable, not calibrated; descriptive by the author's choice | Table 5, Figure 4a |
| A9 | Even the best-scoring 5% of neighborhoods had an estimated false comfort of 68–73%. | `X @5%\|8` for the four estimable scores: N 0.715, TW 0.681, CV 0.680, cv-MDS 0.734. SCD and BOOT are unresolved at 5% and excluded | **model-dependent** | Table 5, Supplementary Table S14, Figure 4a |
| A10 | Neighborhoods should be read against held-out repetitions deep enough for false comfort to be estimable. | Follows from A5–A9 | Recommendation | — |

The abstract has 249 words with its headings (limit 250).

## 2. Highlights

| # | Highlight (characters) | Rests on | Class |
|---|---|---|---|
| H1 | Checks of 2-D maps share the averages' noise; held-out repetitions do not (73) | Section 1 and Supplementary S1.3. The reference's noise is independent of the analyst's by design (Section 2.2) | Argument and design |
| H2 | A simulation study flags where false comfort can be estimated, and where not (76) | Section 3.1, Table 2, Figure 2 | Method result |
| H3 | EEG, 8 trials: an estimated 81% of passing neighborhoods were not faithful (74) | As A5 | **calibrated** |
| H4 | Mouse V1: a 24-trial reference confirmed under 0.1% of passing neighborhoods (76) | As A7 | **finite-reference** |
| H5 | Common checks beat chance 2-4 fold, yet their best 5% had 68-73% false comfort (78) | As A8 and A9; "2-4" rounds 2.0–4.2 | **model-dependent** |

The old highlight "~100% (mouse V1)" was removed (audit item N8): the mouse false comfort is saturated and is not reported as an estimate.

## 3. Conclusions (Section 4.4)

| # | Sentence | Rests on | Class |
|---|---|---|---|
| C1 | Neighborhoods that look faithful are often not faithful: in human EEG at 8 trials, an estimated four in five were not, with a calibrated bound. | As A5 (81%) | **calibrated**. "Often" is supported by EEG at 8 trials only; mouse V1 is finite-reference |
| C2 | Checks computed from the same trials help, but they do not repair this. | "Help": A8 (enrichment). "Do not repair": A9, estimates of 68–73% | **model-dependent**; no threshold is claimed for any check |
| C3 | Whether false comfort can be estimated depends on the reference's depth; the statuses say when; otherwise report the disagreement. | Section 3.1 (Layer A, the stress test, the statuses) and Section 3.4 | Method result |
| C4 | Neighborhoods should be read against held-out repetitions. | Follows from C1–C3 | Recommendation |
| C5 | The estimator and its statuses are released with the code, protocols and results. | The repository | — |

## 4. Other claims checked in the text

These sentences carry a threshold or a comparison and were checked against their class.
- **"Most neighborhoods that passed the plug-in check were not faithful"** (Discussion, first paragraph; heading of Section 3.3). It applies only to EEG t-SNE at 8 trials, the calibrated result.
- **Mouse V1.**
  - It is described only by its disagreement with the reference (≥ 99.7% in every mouse) and its confirmed count.
  - The simulation statement (saturated curves came only with exact false comfort ≥ 0.974 under equal noise) is labeled as context, not an estimate.
- **The checks.**
  - Every check estimate is labeled model-dependent.
  - The pre-specified classification ("leaves most false comfort") appears only in Supplementary Tables S7 and S8, marked as not used for inference.
- **At 16 and 32 trials** only finite-reference values are quoted. That includes the checks at 32 trials, where Section 3.7 quotes disagreement with the 40-trial reference rather than extrapolated values.
- **C1-RSA (cv-MDS)** is model-dependent and unvalidated: its estimate (62–72%) is labeled so, with no threshold claim.
- **The mislabel rate (C3)** is model-dependent in both datasets: estimates only.
- **Generality** (12 variations): the pre-specified verdict is reported, together with the statement that the claim cannot be supported under the classes.

## 5. The audit's required fixes (plan Stage 3, item 1)

| Fix | Where it now stands |
|---|---|
| Figure 2a caption (cv-MDS "calls far more than are faithful") | The figure was replaced. The plug-in check's labels are in Table 4, by class |
| "Four checks passed as many neighborhoods as the plug-in check or more" | Replaced (Section 3.7). TW is rank-matched, SCD and SEED pass more, CV passes fewer; in mouse V1, CV and scDEED pass fewer |
| "Hence" (false comfort from the marginals) | Removed. Section 3.3 states the estimate and the enrichment directly |
| "Reliably distinguishable" (Discussion) | Replaced by "had a positive cross-validated distance" |
| Plasmode truth "itself a finite average" (Limitations) | Now "exact by construction; what it cannot show is how closely real data resemble it" |
| "Fixed representation" (Section 2.1) | The truth is in the representation fitted on block F of each split |
| "Conservative" bound | Now "the pre-specified lower bound L", with its target (the population median), its calibration (Section 3.1) and the sensitivity bound x(2) |
| MVNN used all trials (Section 2.4.1) | Disclosed, as a dependence shared by the analyst's trials and the reference |
| "General" (12 variations) | The pre-specified verdict is reported, with the statement that it cannot be supported (Section 3.6) |
| scDEED null justified by its pass rate | Justified by fidelity to scDEED's procedure; both nulls are reported (Sections 2.8, 2.11, 3.7; S1.4) |
| "Recovered closely" (plasmode) | The plasmode result is now scoped to 8 trials and one map (Sections 2.3.1, 3.1) |
| The guard's scope | Disclosed: computed from the t-SNE cell and applied to every claim (Section 2.3.1) |
| Extrapolant ranges at 16–32 trials and in mouse V1 | Table 4 shows estimates only where the class allows. Supplementary Tables S11, S12, S15 and S16 keep the pre-specified outputs, with statuses and notes |
| "Which the estimator makes practical" | Removed. Practice item 2 states the cost: 48 and 34 extra repetitions per image |
| Section 3.6 title "Why the check fails" | Now "Neighbor distances inside the trial noise (exploratory)", with the causal sentence reworded |
| S1.2 premise "trace terms nearly equal" (N7) | Replaced by a limit statement; the spread (1.63×) is quoted |
| S1.5 certificates (E8) | Reports that no certificate was issued in 900 draws, and drops the multiplicity claim |
| Same-trial argument (R30a) | Limited to what sufficiency implies (Section 4.2.4, S1.3); the protocol's stronger prediction is reported as wrong |
| Highlights and cover letter (N8) | Rewritten to the classes; the cover letter discloses that the reporting rules were set after the first results |
| Total repetitions (R29) | Table 1, Section 2.4.3, Figure 1a |
| Random comparator, enrichment, matched acceptance, continuous overlap (R24, R26, R40) | Sections 2.2.3, 2.8, 2.9, 3.3, 3.7; Table 5; Figure 4; S14 |
| Estimand definitions (R6–R9) | Section 2.1 (undefined at zero passes; protocol-level pooling; truth per representation); Section 2.2.2 (population median) |
| Scope-bounding results in their own section (N9) | Section 3.1, "Where the estimator can be trusted" |
| Title with "faithfulness" (R2) | "Neighborhood faithfulness in maps of trial-averaged neural responses, estimated from held-out repetitions" |
| Morris et al. 2019 (R43) | Cited (Sections 1, 2.3.2; S1.6); the record was checked against Crossref |

## 6. Build checks

- `python3 manuscript/make_tables.py --check` passes.
- `manuscript.tex`, `supplementary.tex` and `cover_letter.tex` build with no LaTeX warnings and no overfull boxes; BibTeX reports 0 warnings.
- The manuscript has 32 pages, the supplement 17 and the cover letter 2.
