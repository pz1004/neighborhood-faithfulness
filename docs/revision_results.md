# Results of the revision analyses (plan Amendment 6)

**What this reports.** The protocol is `docs/revision_protocol.md`; how it was applied is fixed by `docs/revision_addendum.md`. This document reports the analyses of Stages 1 and 2 and gives the classifications that Gate A needs. Nothing here changes the original claims' recorded verdicts. It decides how each result may be reported.

**Status (2026-10-04).** Every analysis is complete except the study-level calibration, Layer C (Section 4), which is still running and marked *pending*.

## 1. What ran

| Analysis | Seeds | Size | Output |
|---|---|---|---|
| Development pilot, with its extension | 6,000,000 | 20 repetitions per cell, 26 cells per dataset | `results/revision/dev/`, Addendum 1 |
| Validation, Layers A and B | 6,500,000 | 100 repetitions per cell, 26 cells per dataset; 2 h | `results/revision/validation/` |
| Gate A classification | — | every real result | `results/gate_a/` |
| Trial manifest and mouse overlap accounting | — | — | `results/manifest/` |
| EEG blocks from separate sessions | 8,100,000 | 12 designs per subject and budget | `results/session_disjoint/` |
| Per-point rerun and benchmark | the committed seeds; 8,200,000 for the alternative scDEED null | every map | `results/points/`, `results/benchmark/` |
| Layer C | 7,000,000 and 7,500,000 | 1,000 virtual studies per run | `results/revision/layer_c/` (*pending*) |

**Deviations and clarifications.** Each was made before the analysis it affects.
- **The noise grid** gained λ = 0.25 and 0.35, which were added rather than replaced (addendum, Section 1).
- **The status thresholds** are the defaults, because no threshold set met the tuning rule (addendum, Section 1).
- **Gate A's open points** were fixed by the addendum (Section 2):
  - accuracy is judged on the smaller extrapolant;
  - matching is within dataset and map;
  - the diagnostic cells never match;
  - the worst matching cell decides;
  - two readings are reported, strict and local.
- **The mislabel rate gets its own status.** Its status and matching use the observables of its own disagreement curve (commit `5e16193`, before the validation was analyzed). Inheriting false comfort's status would have marked an accurately estimated quantity "unresolved".
- **The audit's overlap figure is corrected.** It estimated 0.07 adjacent U–H presentation pairs per image pair in one order; the manifest counts 0.24 over all orders.

## 2. The estimator in simulation (Layer B)

Each cell's median bias over its 100 repetitions (estimate minus exact false comfort), with the RMSE in brackets (`results/revision/validation/analysis.txt`).

**EEG:**

| Cell | Exact | Status | Smaller extrapolant | 1/R_h | Bracket coverage |
|---|---|---|---|---|---|
| T1, λ = 1, 8 trials | 0.826 | estimable | −0.006 (0.033) | +0.046 (0.049) | 0.58 |
| T2, λ = 1, 8 trials | 0.783 | estimable | −0.026 (0.049) | +0.044 (0.047) | 0.70 |
| T3, λ = 1, 8 trials | 0.784 | estimable | −0.012 (0.045) | +0.051 (0.053) | 0.58 |
| T1, λ = 0.5, 8 trials | 0.445 | estimable | −0.081 (0.084) | +0.018 (0.027) | 0.86 |
| T1, λ = 1, 16 trials | 0.646 | unresolved | +0.003 (0.048) | +0.108 (0.109) | 0.49 |
| T3, λ = 1, 32 trials | 0.404 | unresolved | +0.082 (0.120) | +0.240 (0.233) | 0.30 |
| T1, λ = 2 | 0.997 | saturated | −0.000 (0.006) | +0.002 (0.007) | 0.08 |
| Stress T1 | 0.107 | unresolved | +0.162 (0.210) | +0.389 (0.394) | 0.17 |
| Dependent noise, ρ = 0.3 | 0.873 | estimable | −0.123 (0.132) | −0.059 (0.065) | 0.00 |

**Mouse V1:**

| Cell | Exact | Status | Smaller extrapolant | 1/R_h | Bracket coverage |
|---|---|---|---|---|---|
| T1, λ = 1, 8 trials | 0.974 | saturated | +0.008 (0.016) | +0.014 (0.020) | 0.13 |
| T1, λ = 1, 16 trials | 0.905 | estimable | +0.041 (0.046) | +0.055 (0.066) | 0.07 |
| T1, λ = 0.25, 8 trials | 0.333 | estimable | −0.073 (0.085) | +0.012 (0.027) | 0.65 |
| Stress T2 | 0.105 | unresolved | +0.752 (0.743) | +0.803 (0.798) | 0.00 |
| Dependent noise, ρ = 0.3 | 0.994 | estimable | −0.231 (0.254) | −0.170 (0.184) | 0.00 |

**What the simulation shows:**
1. **In the regime of the real EEG data** (λ = 1, 8 trials), the published conservative estimate is nearly unbiased: −0.026 to −0.006, with RMSE ≤ 0.049. 1/R_h overstates by 0.04–0.05.
2. **No extrapolant is right across regimes.**
   - At moderate false comfort, the smaller extrapolant underestimates by 0.04–0.14, and 1/R_h is nearly unbiased.
   - At 16 trials, 1/R_h overstates by 0.11–0.14.
   - At 32 trials, both overstate.
   - In EEG the exact value usually lies between the two, but the published range covers it in only 30–94% of repetitions, depending on the cell.
   - In mouse V1 at 16 trials, both extrapolants overstate, and the range covers the exact value in only 3–7% of repetitions.
3. **Saturated cells were accurate** (RMSE ≤ 0.011 in both datasets).
   - In every mechanism where the analyst's and the reference's trials are equally noisy, a saturated curve came with an exact false comfort of at least 0.974.
   - Only the stress mechanism, whose analyst is far less noisy than its reference, produced high but wrong curves. Even those curves were not saturated.
4. **The stress test confirms the review's concern.**
   - A shallow, noisy reference against an accurate analyst inflates the estimate by +0.16 to +0.75.
   - The statuses flag every such cell as unresolved.
5. **Noise shared by blocks U and H biases the estimate downward,** by −0.12 in EEG and −0.23 in mouse at ρ = 0.3, as audit item N4 expected. It cannot produce a high false comfort.
6. **The mislabel rate is recovered well in the regime of the real data**, with bias within ±0.008. At moderate false comfort it is underestimated by up to 0.07.
7. **By status, pooled over cells:**

   | Status | Dataset | Smaller extrapolant: median bias (RMSE) | 1/R_h: median bias (RMSE) |
   |---|---|---|---|
   | Estimable | EEG | −0.032 (0.058) | +0.005 (0.033) |
   | Estimable | mouse V1 | −0.028 (0.095) | +0.000 (0.071) |
   | Unresolved | both | — (0.15–0.32) | — (0.21–0.34) |
   | Saturated | both | — (≤ 0.011) | — |

## 3. Deeper references (Layer A)

Layer A keeps the first split's maps fixed and regenerates only the references, 20 times, at the real grid and at 2× and 4× its largest size.

**The disagreement keeps falling beyond the real grid.** In EEG (T1, λ = 1, 8 trials), it falls from 0.968 (8 trials) through 0.922, 0.883 and 0.849 to 0.832 (160 trials), against an exact value of 0.801.
- **Within the real grid.** The R_h^−1/2 extrapolant fitted to {8, 20, 40} lands within 0.004 of the exact value.
- **Adding the deep sizes** moves it to −0.020.

**The stress failure is a shallow-reference failure in EEG.**
- In EEG, references of 80 and 160 trials bring the stress cell's bias from +0.149 to +0.017.
- In mouse V1 they do not: the bias falls from +0.553 to +0.276, and from +0.757 to +0.495.

## 4. Study-level calibration (Layer C) — *pending*

**What runs, under the local reading:**
- **EEG C1:** false comfort at 8 trials, with the guard, under T1 and T3.
- **Mouse V1 C3:** the mislabel rate at 8 and 16 trials, under T1 and T3. In mouse V1, T3 is T2.

Each run has 1,000 virtual studies of 10 subjects, drawn from 30 virtual subjects whose truths come from 200 oracle datasets each.

**How Layer C is judged.** These rules were fixed and committed before any Layer C output was analyzed. Protocol Section 5 leaves item 2 open; the other items restate the protocol and the addendum.
1. **The criterion** is the protocol's: Pr(L > Θ) ≤ 0.025 + 2 Monte Carlo errors over the 1,000 studies.
   - **For EEG C1,** it applies to L of false comfort at 8 trials.
   - **For mouse V1 C3,** it applies to the mislabel rate's bound L_q at 8 trials and at 16, since C3 claims both budgets. As published, L_q has no guard.
2. **Both scenarios must meet it.** The worst scenario decides, just as the worst matching cell decides in Gate A (addendum, Section 2, item 3). A claim calibrated under T1 only remains model-dependent.
3. **The guard ran in the first 200 studies only** (addendum, Section 2, item 5). The other 800 studies use L without it. The guard can only lower L, so this makes the criterion harder to meet, not easier. If the guard triggers, the rate over the first 200 studies is reported too.
4. **The order-statistic bound** (the 2nd-smallest subject estimate) is evaluated beside the published bound (protocol Section 3.4). The author may choose it at Gate A. It would be reported as a revision analysis, not as the pre-specified bound.
5. **Reported but not decisive:**
   - Pr(U < Θ);
   - the coverage and width of [L, U];
   - the rate of declaring L > 0.5, or L_q > 0.10;
   - the same rates against each study's own sample median.

## 5. Revision analyses on the real data

### 5.1 Trial manifest and overlap (`results/manifest/`)

**EEG: chronological order within each session is recoverable.**
- **Source.** The dataset authors' preprocessing code (github.com/gifale95/eeg_encoding, `02_eeg_preprocessing`) selects each test image's 20 presentations per session with `shuffle(..., random_state=20200220)`. That is the seed of their cross-session shuffle, which reach/things.py already replays.
- **One permutation for every image.** With exactly 20 presentations per image and session, the permutation (11, 1, 17, 0, …) is the same for every image. The manifest therefore gives every stored repetition its chronological rank within its session.
- **Onset times** are not recoverable from the preprocessed files.
- **The review's assumption.** It took only session membership to be known.

**Mouse V1:**
- **Overlap of response windows.** Consecutive presentations start every 8 frames (30.1 Hz). Consecutive 0.5 s windows overlap by 7 of 15 frames; 0.25 s windows do not overlap.
- **Pairs straddling U and H.** About 28% of consecutive presentation pairs fall one in U and one in H. That is 0.24 per image pair over all orders, against 16 × 24 = 384 pairs of U and H trials per image pair.

### 5.2 EEG blocks from separate sessions (`results/session_disjoint/`)

**Design.**
- F is drawn from one session and U from another; H is the other two sessions, 40 repetitions.
- The control is a balanced random split with the same block sizes.
- The primary map is used throughout.

**Results** (medians over 10 subjects):

| Budget | Design | d at R_h = 8, 20, 40 | Smaller extrapolant | Random selection | Status |
|---|---|---|---|---|---|
| 8 trials | control | 0.970, 0.913, 0.881 | 0.810 | 0.937 | estimable |
| 8 trials | separate sessions | 0.972, 0.930, 0.900 | 0.848 | 0.945 | estimable |
| 16 trials | control | — | 0.632 | — | unresolved |
| 16 trials | separate sessions | — | 0.696 | — | estimable |

**Paired difference per subject** (separate sessions minus control): +0.052 at 8 trials (range −0.019 to +0.103) and +0.076 at 16.

**Reading.**
- Separating the sessions *raises* the estimate. The session-balanced design shares session-specific components between the analyst's and the reference's trials, and shared noise makes them agree (Section 2, item 5).
- The published estimates are therefore, if anything, on the low side of a cross-session estimate.
- Across sessions the responses themselves can also change, so this is a sensitivity analysis, not a measurement of leakage.
- The control reproduces the published estimate (0.810 against 0.811).

### 5.3 Per-point rerun and benchmark (`results/points/`, `results/benchmark/`)

**Reproduction.** The rerun re-created every map with the committed seeds. It reproduced every committed count of the comparison of checks and of the cv-MDS cell, in both datasets (`reproduction.json`: none differs).

**EEG at 8 trials.** These are medians over subjects; the *total* column is the repetitions an assessment uses (F + R + H).

| Check | Accepts | Total | Estimated false comfort (smaller) | Random selection | Enrichment | Finite-reference d_40 | Status |
|---|---|---|---|---|---|---|---|
| plug-in | 18.8% | 56 | 0.811 | 0.934 | 2.8 | 0.877 | estimable |
| trustworthiness | 19.0% | 56 | 0.792 | 0.934 | 3.1 | 0.865 | estimable |
| scDEED | 32.5% | 56 | 0.863 | 0.934 | 2.2 | 0.912 | estimable |
| scDEED, alternative null | 8.2% | 56 | 0.702 | 0.934 | 4.5 | 0.804 | estimable |
| cross-validated distances | 7.9% | 56 | 0.728 | 0.935 | 4.2 | 0.824 | estimable |
| seed recurrence | 34.6% | 56 | 0.870 | 0.934 | 2.0 | 0.916 | estimable |
| split-half recurrence | 1.2% | 56 | 0.467 | 0.934 | 5.7 | 0.654 | unresolved |
| bootstrap recurrence | 0.0% | 56 | — | — | — | — | too few acceptances |
| cv-MDS plug-in label | 2.0% | 56 | 0.623 | 0.964 | 10.2 | 0.754 | estimable |

**Matched acceptance.** Each map accepts the 5% of points with the best score, ties included. The estimated false comfort is then:

| Check | False comfort | Note |
|---|---|---|
| plug-in (smallest naive reach) | 0.715 | |
| trustworthiness | 0.681 | |
| scDEED score | 0.635 | unresolved |
| cross-validated reach | 0.680 | |
| bootstrap support | 0.560 | unresolved |
| cv-MDS naive reach | 0.734 | |

At 50% acceptance every check is near the random selection, 0.88–0.93.

**What the benchmark shows:**
1. **Same-trial checks carry information.** Every check selects faithful neighborhoods 2 to 10 times more often than a random selection of as many neighborhoods, so "cannot reveal where the averages are wrong" (abstract) is too strong.
2. **They are still mostly wrong.** Faithful neighborhoods are rare (about 6% in EEG at 8 trials), so even the strictest 5% by any score is mostly not faithful.
3. **The scDEED comparison depends on the null.** The registered null (measured features permuted) is permissive. The alternative null (representation coordinates permuted) passes a quarter as many neighborhoods, and its precision is twice the registered null's. The comparison's verdict for scDEED therefore depends on that choice.
4. **The cost of the assessment.** It needs F + H = 48 repetitions in EEG beside the analyst's R, and 34 in mouse V1.
5. **Continuous overlap with the deepest reference** (the share of displayed neighbors among the reference's k′ nearest):
   - in EEG, 0.38 for all neighborhoods and 0.51 for the plug-in check's passes;
   - in mouse V1, 0.22 and 0.26, against about 0.17 by chance.

   A faithful neighborhood needs 0.8.
6. **Mouse V1.** Every check is saturated, or has too few acceptances (bootstrap and split-half recurrence). No enrichment can be determined.

### 5.4 cv-MDS diagnostics

| | EEG, 8 trials | EEG, 16 | EEG, 32 | Mouse V1, 8 trials | Mouse V1, 16 |
|---|---|---|---|---|---|
| Negative cross-validated squared distances | 4.9% | 1.1% | 0.1% | 25.1% | 16.7% |
| Kruskal stress-1 | 0.30 | 0.27 | 0.26 | 0.37 | 0.31 |

- **Zeros in the MDS input.** Every negative squared distance becomes a zero in √max(0, d̂²).
- **Ranking convention.** The naive labels rank the signed d̂², so they have no ties at zero. The references use plug-in distances.

## 6. Every real result's status and class (Gate A)

`results/gate_a/classification.txt` lists every cell. The claims:

**EEG:**

| Claim | Published range | Status | Matching cells (worst bias / RMSE) | Strict | Local |
|---|---|---|---|---|---|
| C1, false comfort, 8 trials | 0.811–0.865 | estimable | 3 (−0.026 / 0.049) | model-dependent | Layer C |
| C1-RSA, cv-MDS, 8 trials | 0.623–0.720 | estimable | none | model-dependent (unvalidated) | model-dependent (unvalidated) |
| C2, false comfort, 32 trials | 0.476–0.635 | unresolved | — | finite-reference | finite-reference |
| C3, mislabel rate, 8 trials | 0.169–0.172 | estimable | 8 (−0.035 / 0.039) | model-dependent | model-dependent |
| C3, mislabel rate, 16 trials | 0.208–0.219 | estimable | 6 (−0.038 / 0.038) | model-dependent | model-dependent |
| C3, mislabel rate, 32 trials | 0.210–0.243 | estimable | 6 (−0.059 / 0.060) | finite-reference | finite-reference |

**Mouse V1:**

| Claim | Published range | Status | Matching cells (worst bias / RMSE) | Strict | Local |
|---|---|---|---|---|---|
| C1, false comfort, 8 trials | 1.000 | saturated | — | finite-reference | finite-reference |
| C1-RSA, cv-MDS, 8 trials | 1.000 | saturated | — | finite-reference | finite-reference |
| C3, mislabel rate, 8 trials | 0.161 | estimable | 10 (−0.012 / 0.014) | model-dependent | Layer C |
| C3, mislabel rate, 16 trials | 0.172 | estimable | 9 (−0.005 / 0.010) | model-dependent | Layer C |

**Finite-reference values** (median disagreement with the deepest reference, 95% subject-bootstrap interval):

| Result | Value |
|---|---|
| EEG C1 | 0.877 [0.768, 0.928] |
| EEG C2 | 0.673 [0.500, 0.739] |
| Mouse V1 C1 | 1.000 [0.998, 1.000] |

**The 12 variations.**
- **EEG, at 8 trials:**
  - estimable with accurate matching cells, as C1: early window, k = 5, perplexities 15 and 50;
  - estimable but unmatched: ε = 0.1 and 0.3, MDS of Euclidean distances, and cv-MDS;
  - unresolved: k = 20 and PCA-2;
  - saturated: unwhitened.
- **Mouse V1:** saturated in every variation.
- **The generality claim.** C1 was to count as general if it held in at least 10 of 12 variations. At most 5 EEG variations can carry a calibrated threshold claim, so under the revision's rules that claim cannot be supported. Robustness becomes descriptive.

**The checks at 8 trials.**
- **EEG:**
  - estimable with accurate matching cells: trustworthiness, scDEED, cross-validated distances and seed recurrence;
  - unresolved: split-half recurrence;
  - too few acceptances: bootstrap recurrence.
- **Mouse V1:** every check is saturated or has too few acceptances.
- **Under the local reading,** the checks' pre-registered classification would need Layer C for the checks. The protocol makes that conditional on the author keeping the classification as an inferential claim.

## 7. Decisions for the author (Gate A)

*To be completed when the pending analyses are in.*
