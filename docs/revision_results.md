# Results of the revision analyses (plan Amendment 6)

**What this reports.** The protocol is `docs/revision_protocol.md`; how it was applied is fixed by `docs/revision_addendum.md`. This document reports the analyses of Stages 1 and 2 and gives the classifications that Gate A needs. Nothing here changes the original claims' recorded verdicts. It decides how each result may be reported.

**Status (2026-10-05).** Every analysis is complete, including the study-level calibration (Layer C, Section 4). Section 7 sets out the decisions for the author at Gate A.

## 1. What ran

| Analysis | Seeds | Size | Output |
|---|---|---|---|
| Development pilot, with its extension | 6,000,000 | 20 repetitions per cell, 26 cells per dataset | `results/revision/dev/`, Addendum 1 |
| Validation, Layers A and B | 6,500,000 | 100 repetitions per cell, 26 cells per dataset; 2 h | `results/revision/validation/` |
| Gate A classification | — | every real result | `results/gate_a/` |
| Trial manifest and mouse overlap accounting | — | — | `results/manifest/` |
| EEG blocks from separate sessions | 8,100,000 | 12 designs per subject and budget | `results/session_disjoint/` |
| Per-point rerun and benchmark | the committed seeds; 8,200,000 for the alternative scDEED null | every map | `results/points/`, `results/benchmark/` |
| Layer C | 7,000,000 and 7,500,000 | 1,000 virtual studies per run, 6 runs; about 5 h per EEG run | `results/revision/layer_c/` |
| Layer C, the bound's rule (descriptive, added after the runs were read) | the Layer C resamples | every run | `results/revision/layer_c/rule.txt` |

**Deviations, clarifications and corrections.** Each deviation and clarification was made before the analysis it affects. The corrections fix statements in the audit.
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
- **The audit's item E5 is corrected** (Section 4). The published bound is not, as a rule, the 3rd-smallest subject estimate. It is slightly liberal: on continuous populations it errs in 3.0–3.4% of studies, against a nominal 2.5%, not 5.5%.

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

## 4. Study-level calibration (Layer C)

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

**Results** (`results/revision/layer_c/calibration_*.txt`; ± is one Monte Carlo error).

**EEG C1: false comfort at 8 trials, with the guard.**

| Truth | Θ (range of θ_v) | Study medians: smaller / 1/R_h | Pr(L > Θ) | Criterion | Coverage of [L, U] | Mean width | Declares L > 0.5 |
|---|---|---|---|---|---|---|---|
| T1 | 0.821 (0.550–0.988) | 0.811 / 0.867 | 0.020 ± 0.004 | met (≤ 0.034) | 0.976 | 0.27 | 97.6% |
| T3 | 0.794 (0.453–0.985) | 0.777 / 0.844 | 0.016 ± 0.004 | met (≤ 0.033) | 0.982 | 0.34 | 84.9% |

- **EEG C1 is calibrated under both truths.** Under the local reading it is therefore "estimable, calibrated" (Section 6).
- **The guard never triggered** in the first 200 studies (median plasmode bias −0.008 under T1 and −0.004 under T3), so L was the same with and without it.
- **The order-statistic bound** errs less often (Pr(L > Θ) 0.011 and 0.006) but declares L > 0.5 less often (84.8% and 68.0%). The published bound needs no replacement.
- **Against each study's own sample median,** L exceeds it in 0.1% of studies.
- **The mislabel rate,** for completeness: Pr(L_q > Θ_q) is 0.017 and 0.010. EEG C3 nevertheless stays model-dependent, because its matching cells' worst bias (−0.035) exceeds 0.03.

**Why the bound met or missed its level, and a correction to the audit.** This comes from `results/revision/layer_c/rule.txt` (`analyze_layer_c_rule.py`, descriptive, added after these runs were read).
- **Audit item E5 is wrong as a general rule.** It read L as the 3rd-smallest of the 10 subject estimates, which would err with probability 0.055 even if every estimate were exact.
  - That held for the real data's spacing (C1: L = x(3) = 0.642).
  - Across the virtual studies, L equals x(3) in 41–53% of them and lies below it otherwise.
- **The bootstrap rule itself is slightly liberal.** Its nominal error rate is 2.5%.
  - On continuous populations of known median, it errs in 3.0–3.4% of studies (normal, uniform and exponential; 20,000 studies each).
  - On the exact values of the virtual populations, which are discrete, it errs in 2.1–4.3%.
- **In EEG the estimator's slight downward bias compensates.** Its median within-study bias is −0.006 under T1 and −0.014 under T3, which brings the error rate to 2.0% and 1.6%.
- **What this settles.**
  - Audit item N3 said "conservative" presumed a coverage never shown. For EEG C1 the coverage has now been shown.
  - The margin, however, rests on the estimator erring slightly low in this regime, as it does in the matching cells (Section 2).
  - E5's concern was right in direction but wrong in size and mechanism.

**Mouse V1 C3: the mislabel rate at 8 and 16 trials, without the guard.**

| Truth | Budget | Θ_q (range of θ_q,v) | Median within-study bias | Pr(L_q > Θ_q) | Criterion | Declares L_q > 0.10 |
|---|---|---|---|---|---|---|
| T1 | 8 trials | 0.138 (0.129–0.155) | −0.000 | 0.040 ± 0.006 | **not met** (≤ 0.037) | 100% |
| T1 | 16 trials | 0.139 (0.120–0.149) | +0.002 | 0.066 ± 0.008 | **not met** (≤ 0.041) | 100% |
| T3 (= T2) | 8 trials | 0.165 (0.131–0.266) | −0.003 | 0.024 ± 0.005 | met (≤ 0.035) | 100% |
| T3 (= T2) | 16 trials | 0.210 (0.142–0.307) | −0.006 | 0.013 ± 0.004 | met (≤ 0.032) | 100% |

- **Mouse C3 stays model-dependent.** It met the criterion under T3 at both budgets and missed it under T1 at both. The worst run decides (item 2 above).
- **Why it missed under T1.** The virtual subjects' mislabel rates are tightly clustered under T1, with an interquartile range about 0.008 wide, so estimation error dominates the spread of the subjects' estimates.
  - Nothing offsets the rule's slight liberality at 8 trials, where the bias is −0.0002, and the error rate was 4.0%.
  - A small upward bias at 16 trials (+0.002) took it to 6.6%.
  - Under T3 the rates spread over 0.13–0.31, the bias is slightly negative, and the bound met the criterion.
- **The threshold itself was never wrongly declared.** Every virtual subject's mislabel rate exceeds 0.10 in all four runs, so no declaration that L_q > 0.10 was wrong. The miss concerns the bound's coverage of Θ_q, which is what the criterion tests.
- **False comfort in mouse V1 is not a claim here,** since its real curves are saturated, but the same runs show the following.
  - At 8 trials, Θ is 0.974 under T1 and 0.978 under T3. The smaller extrapolant overstates by +0.004 within studies, and Pr(L > Θ) is 0.049 and 0.059.
  - At 16 trials, Θ is 0.906 and 0.889. It overstates by +0.026 and +0.036, and Pr(L > Θ) is 0.127 and 0.213, in line with Layer B's +0.041 at 16 trials (Section 2).
  - So the published bound would not be calibrated for mouse false comfort at either budget.

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
| C1, false comfort, 8 trials | 0.811–0.865 | estimable | 3 (−0.026 / 0.049) | model-dependent | estimable, calibrated (Layer C, both truths) |
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
| C3, mislabel rate, 8 trials | 0.161 | estimable | 10 (−0.012 / 0.014) | model-dependent | model-dependent (Layer C missed under T1) |
| C3, mislabel rate, 16 trials | 0.172 | estimable | 9 (−0.005 / 0.010) | model-dependent | model-dependent (C3 claims both budgets; Section 4) |

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

### 7.1 Where each claim stands

Classes under the local reading, after Layer C. Under the strict reading the one calibrated entry, EEG C1, becomes model-dependent; nothing else changes.

| Dataset | Claim | Class | What the paper can say |
|---|---|---|---|
| EEG | C1: false comfort of t-SNE maps at 8 trials | estimable, calibrated | an estimate of 0.81 (published range 0.81–0.87), and that it exceeds 0.5 (L = 0.642) |
| EEG | C1-RSA: cv-MDS at 8 trials | model-dependent (unvalidated) | an estimate of 0.62–0.72, labeled model-dependent, with no threshold claim |
| EEG | C2: 32 trials | finite-reference only | disagreement with a 40-trial reference: 0.673 [0.500, 0.739] |
| EEG | C3: the mislabel rate | model-dependent at 8 and 16 trials; finite-reference only at 32 | estimates of 0.17 (8 trials) and 0.21–0.22 (16), with no threshold claim; the finite-reference disagreement at 32 |
| Mouse V1 | C1 and C1-RSA | finite-reference only | disagreement with a 24-trial reference: 1.000 [0.998, 1.000]; it confirmed 0.08% of the displayed faithful neighborhoods (5 of 5,919 calls, pooled over subjects, maps and reference sub-blocks) |
| Mouse V1 | C3 | model-dependent (Layer C met under T3, missed under T1) | estimates of 0.16–0.17, with no threshold claim |
| Both | the checks | EEG: model-dependent for four checks (Layer C not run), finite-reference for split-half recurrence, too few acceptances for bootstrap recurrence; mouse V1: finite-reference or too few acceptances | enrichment over a random selection and acceptance against error (Section 5.3), descriptively |
| Both | generality over the 12 variations | cannot be supported | robustness, descriptively |

**The stop rule does not apply.** In EEG the checks select faithful neighborhoods 2–10 times more often than a random selection. The budgets also differ: the disagreement with the same 40-trial reference is 0.877 at 8 trials and 0.673 at 32.

### 7.2 Decisions

1. **Strict or local reading.** Only EEG C1 depends on the choice.
   - **Recommendation: local.**
     - The threshold claim stays within the regime its matching cells cover: EEG, t-SNE, 8 trials.
     - There it is accurate under three truth families, and its bound is calibrated under both population scenarios.
     - Shared noise would bias it downward (Section 2, item 5), and separating the sessions raised it (Section 5.2).
     - The strict reading would drop a validated result because the estimator fails in other regimes, which the statuses already detect.
2. **Direction for each dataset.**
   - **EEG, recommended: a scoped latent claim.**
     - C1 is the calibrated result at 8 trials.
     - C1-RSA and C3 are reported as model-dependent estimates.
     - At 16 and 32 trials, the result is the finite-reference disagreement curve, extrapolated nowhere.
   - **Mouse V1, recommended: finite-reference.**
     - The statement is "a held-out reference of 24 trials confirmed fewer than 0.1% of the neighborhoods the maps showed as faithful". The plan's wording, "no displayed neighborhood was confirmed", is too strong: 5 of 5,919 were (`results/gate_a/classification.txt`).
     - It comes with the simulation's context: wherever the analyst's and the reference's trials were equally noisy, a saturated curve came with an exact false comfort of at least 0.974 (Section 2, item 3). That context is evidence about the method, not an estimate for these mice.
     - C3 is reported as a model-dependent estimate.
     - The alternative is to move mouse V1 to the supplement. I advise against it: this is the case that shows why the statuses are needed.
3. **The checks' pre-registered classification** ("leaves most of the false comfort").
   - Keeping it as an inferential claim needs Layer C for each check: 200 studies × 2 truths × 5 checks, about a day of compute.
   - **Recommendation: make it descriptive,** and lead with the benchmark. It answers the practical question directly. In EEG, every check that accepts anything beats a random selection. Yet even its strictest 5% has an estimated false comfort of 0.68–0.73 where the status is estimable; two scores are unresolved there.
4. **The bound.** The published bound was calibrated where it is claimed, so it stays.
   - Its margin rests on the estimator's slight downward bias in this regime, since the rule alone errs in about 3% of studies (Section 4).
   - The order-statistic bound gives the same conclusion for C1 (x(2) = 0.605 > 0.5) and errs in about 1% of studies even with exact values.
   - **Recommendation:** keep L, and report x(2) beside it as a sensitivity bound.

**These follow under any direction,** and are done in Stage 3:
- the corrections in Section 4 (audit item E5) and Section 1;
- "cannot reveal" in the abstract (Section 5.3, item 1);
- generality as descriptive (Section 6);
- the mouse "≈100%" as a finite-reference statement (audit items N2 and N8);
- the scDEED comparison reported under both nulls (Section 5.3, item 3);
- the title, chosen at Gate B, using "faithfulness" (audit item R2).

### 7.3 The author's decisions (2026-10-05)

The author chose every recommendation above:
1. **The local reading.**
2. **EEG:** a scoped latent claim. C1 is the calibrated result at 8 trials; C1-RSA and C3 are model-dependent estimates; at 16 and 32 trials only the finite-reference disagreement is reported.
3. **Mouse V1:** in the main text, as a finite-reference result with the simulation's context; C3 as a model-dependent estimate.
4. **The checks:** descriptive. The benchmark leads, and the inferential classification is dropped.

**The bound** stays the published L, with the 2nd-smallest estimate reported beside it as a sensitivity bound. That was the stated default, and the author did not object.

The rewrite (Stage 3) follows these decisions and ends at Gate B.
