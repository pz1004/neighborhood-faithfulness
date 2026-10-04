# Protocol: validating the false-comfort estimator, and the revision analyses

**Version.** Version 1 (draft for the author's review), 2026-10-04. Plan Amendment 6 points here.

**Pre-registration.** As with every earlier protocol, this one, the code that implements it and its synthetic tests are committed before any of its analyses run. Two further commits follow before the validation runs (Section 8):
- an addendum, fixing what the development pilot decides;
- the validation code itself.

**Scope.** This protocol answers the colleague's revision plan (`docs/review_comments.md`). It follows the audit of that plan (`docs/review_audit.md`). It also absorbs the review's proposed "estimand specification" and "unit-test report": Section 2 is the first, and the test log of Section 7 is the second.

## 1. Purpose, and what is known before this protocol

**Purpose.** The manuscript at tag `jnm-baseline` reports latent false comfort, extrapolated from held-out references, with subject-bootstrap bounds. The audit found four gaps (its items E1–E8):
- the estimator was validated at one budget and for one map only;
- its larger extrapolant overstates at 16 and 32 trials;
- the mouse references are saturated;
- the bound L is under-covered for a population median.

This protocol does three things:
1. **Validate** the estimator end to end, and find where it fails (Section 3).
2. **Define** observable statuses that say, for any cell, whether a latent estimate may be reported (Section 4).
3. **Fix in advance** the rules by which the revised manuscript reports each result (Section 5).

It also specifies the revision's descriptive analyses on the real data (Section 6).

**What is known before this protocol:**
- Every committed result:
  - the confirmatory runs;
  - the comparison of checks;
  - the pilot.
- Items E1–E8 of the audit, computed on 2026-10-04 by `analyze_review_audit.py` from the committed rows, before this protocol was written.
- Nothing from the simulations of Section 3, or from the real-data analyses of Section 6, has been computed.

**What kind of analysis this is.** Everything here is a *revision analysis*: designed after the results and E1–E8 were seen, but specified before any of its own runs. The paper labels each analysis as one of three kinds:
- **original pre-specified:** Amendments 3–5, unchanged;
- **revision:** this protocol;
- **exploratory.**

The original claims (C1, C1-RSA, C2, C3 and the checks' classification) keep their recorded verdicts. This protocol does not re-decide them. It decides how they may be reported.

## 2. Estimands

**Notation.**
- Objects (stimuli) are i = 1, …, n.
- A *split* z is a random draw of blocks F, U, H. The representation W_z is fitted on block F of that split, and the truth is μ_i projected by W_z.
- The truth label T_iz of a displayed neighborhood depends on z. "Faithful" is always relative to the representation in which it is judged.
- For a check X, A_i = 1 if X accepts (passes) the neighborhood of i.
- T^(r)_i is the label from a reference of r held-out trials.

**Quantities:**
- **Map level.** For one map M:
  - FC_X(M) = Σ_i A_i (1 − T_i) / Σ_i A_i;
  - it is *undefined* (not zero) when Σ_i A_i = 0.
- **Protocol level.** For subject s and budget R:
  - θ_sX(R) = E_Z[Σ_i A_i (1 − T_i)] / E_Z[Σ_i A_i];
  - Z covers the draw of the blocks, the representation fitted on F, the map seed and the check's own randomness;
  - θ weights maps by how many neighborhoods they accept;
  - the published per-subject estimate is the ratio of sums over 10 splits × 3 seeds, extrapolated over reference sizes. It estimates θ, not the false comfort of any one map. The difference is small here (audit E6) but it is defined.
- **Acceptance and yield.**
  - acceptance PASS_sX = E_Z[Σ A_i] / n;
  - faithful yield Y_sX = E_Z[Σ A_i T_i] / n;
  - unfaithful yield U_sX = PASS_sX − Y_sX;
  - θ_sX = U_sX / PASS_sX when PASS_sX > 0;
  - all of these are computed within subject before any median is taken.
- **Population and sample targets.**
  - The inferential target of every interval is Θ_X(R), the median of θ_sX(R) over the population of subjects that the 10 represent.
  - The median over the 10 observed subjects is descriptive.
- **Finite-reference target.**
  - d^(r)_sX = E_Z[Σ A_i (1 − T^(r)_i)] / E_Z[Σ A_i] for a reference of r held-out trials.
  - It is reported as "disagreement with an r-trial reference" and never called false comfort.
- **Random comparator and enrichment.**
  - The comparator selects, in each map, as many neighborhoods at random as X accepts. Its expected false comfort is FC_rand,X.
  - Enrichment is (1 − FC_X) / (1 − FC_rand,X), not determined when 1 − FC_rand,X < 0.005.
  - Both are extrapolated exactly like FC_X.
- **Continuous overlap.**
  - o^(r)_i = |E_i ∩ N_k′^(r)(i)| / k, where N_k′^(r)(i) is the set of the k′ objects nearest to i in the reference.
  - A neighborhood is faithful at the reference exactly when o ≥ m/k.
  - It is a finite-reference measure.

## 3. Simulation study (ADEMP; Morris, White & Crowther 2019)

### 3.1 Aims
1. **Accuracy.** How accurately the published extrapolants, and any revised estimator (Section 3.5), recover θ:
   - at R ∈ {8, 16, 32} (EEG) and {8, 16} (mouse V1);
   - for the t-SNE and cv-MDS maps.
2. **Failure regimes.** Where the estimator fails, and whether the observable statuses of Section 4 flag the failure.
3. **Coverage.** The coverage of the published bounds L and U, and of an order-statistic alternative, for the population target Θ.

### 3.2 Data-generating mechanisms

Each dataset's mechanisms are built from each real subject's split-0 blocks, the published plasmode's split. All of them work in the *measured* features: EEG channel × time; mouse cells. This lets every simulated dataset refit the representation on its own block F.

**Truth families:**
- **T1, the published plasmode truth:** μ = the mean of block H.
- **T2, low-rank:** μ = m + P_r(μ_T1 − m).
  - m is the mean over images.
  - P_r projects onto the first r = 5 principal components of the image means μ_T1 − m.
  - The result has smooth, better separated structure.
- **T3, gap-matched:** μ = μ_T2 + α(μ_T1 − μ_T2), with α ∈ [0, 1].
  - **How α is fixed:** by the development pilot. The median relative gap between each image's 8th and 20th nearest true distance, in the base split's representation, must match the cross-validated estimate from the real data: 0.159 in EEG and 0.298 in mouse V1 (`results/*/explore_snr.txt`).
  - **If no α matches,** the nearest is used and recorded.

**Noise.**
- A trial is μ + λ·e.
- e is drawn with replacement from block U's within-image residuals, rescaled by √(n_U/(n_U − 1)) as in the published plasmode, independently for every trial.
- **Screening levels:** λ ∈ {0.1, 0.5, 1, 2}.
- **Zero-noise control:** blocks U and H are noise-free (λ = 0), and block F keeps λ = 1 so that the whitening can be fitted.

**Design.**
- Each dataset keeps its real design: 80 trials per image in 4 sessions (EEG), or 50 in 2 time strata (mouse).
- The published splitter, halves, budgets and reference grids are used.

**Dependent noise (mechanism D).**
- **The offset:** for every image i and session (EEG) or stratum (mouse) g, add s_ig = √ρ·λ·ẽ_ig to every trial of image i in g.
  - ẽ_ig is one residual drawn from image i's pool.
  - The trial-level noise is scaled by √(1 − ρ), so the total noise per trial is unchanged.
- **Levels:** ρ ∈ {0.1, 0.3}, at T1 and λ = 1.
- **The truth** is μ, without the offsets. Because blocks are balanced over sessions or strata, U and H share the mean offset.
- **The audit's expectation, under test (N4):** the estimate falls below the exact value.

**Stress test (mechanism S).** This is a conditional-reference diagnostic, not an estimate of how often the situation arises in the real data.
- **Setup:**
  - the analyst's block U is drawn at λ_U = 0.1, so most displayed neighborhoods are truly faithful;
  - blocks F and H are drawn at λ = 1;
  - run at T1 and T2.
- **Question:** does a shallow, noisy reference make the estimator report a high false comfort where the exact value is low?

### 3.3 Estimands in the simulation
- **Per repetition:** the exact false comfort is computed from exact labels. The labels are taken against μ, projected through each split's fitted representation, and pooled over the repetition's maps, as θ is.
- **Layer C:** a virtual subject's θ_v is computed by oracle Monte Carlo:
  - 200 independent datasets × 1 random split × 3 seeds, with exact labels;
  - its Monte Carlo error is reported.

### 3.4 Validation layers

**Layer A: references only.**
- **Fixed:** the analyst's data, the fitted representation, the maps and the acceptance decisions, from Layer B's repetitions.
- **Regenerated, K = 20 times:** only the references.
- **Reference sizes:**
  - the real grid: {8, 20, 40} EEG, {8, 12, 24} mouse;
  - plus 2× and 4× its largest size: 80 and 160 EEG, 48 and 96 mouse;
  - the published numbers of sub-blocks, and one sub-block at the deeper sizes.
- **Outputs:**
  - the extrapolation error against the exact value of the fixed maps;
  - whether the disagreement curve keeps the extrapolant's form beyond the real grid;
  - how the finite-depth error relates to the extrapolation error.

**Layer B: end to end.**
- **One repetition is one virtual subject:**
  - a new dataset from the mechanism;
  - the published procedure on it, which is:
    - 10 random splits;
    - a representation refitted on F in each split;
    - 3 map seeds;
    - naive labels;
    - references at the real grid;
  - exact labels;
  - per-subject estimates from both extrapolants and any candidate of Section 3.5;
  - the observable statuses of Section 4.
- **Reference trials** beyond the design, for the deeper reference sizes, are drawn from the same mechanism and used only for the depth diagnostics.
- **Size:** 100 repetitions per cell, 10 for each base subject.

**Cells of Layers A and B, per dataset:**

| Cells | Map | R | Count |
|---|---|---|---|
| T1–T3 × λ ∈ {0.1, 0.5, 1, 2} (screening) | t-SNE | 8 | 12 |
| zero-noise control, T1 | t-SNE | 8 | 1 |
| D at ρ ∈ {0.1, 0.3}, T1, λ = 1 | t-SNE | 8 | 2 |
| S at T1 and T2 | t-SNE | 8 | 2 |
| T1–T3 at λ = 1 | cv-MDS | 8 | 3 |
| T1 and T3 at λ = 1 | t-SNE | 16 and 32 (EEG); 16 (mouse) | 4 (EEG); 2 (mouse) |

That is 24 cells in EEG and 22 in mouse V1.

**Layer C: study-level calibration.**
- **When it runs:** only for claims whose cells are "estimable" after Layers A and B and the real-data statuses (Section 5). The likely candidates are EEG C1 (t-SNE) and C1-RSA (cv-MDS) at R = 8.
- **Population:**
  - 30 virtual subjects: the 10 base subjects × λ ∈ {0.8, 1.0, 1.25};
  - built under T1 (scenario 1) and under T3 (scenario 2);
  - Θ is the median of their θ_v.
- **Virtual studies:** 1,000 per scenario and map. Each study:
  - samples 10 virtual subjects with replacement;
  - draws a new dataset for each;
  - reproduces the published procedure exactly:
    - 10 splits × 3 seeds;
    - references at the real grid;
    - both extrapolants;
    - per-subject estimates and their median;
    - 10,000 bootstrap resamples of subjects, giving L and U;
    - the not-evaluable rule;
    - the bias guard, from a 20-draw plasmode per subject, as published.
  - **If the guard makes a study too slow** (over 6 h of wall time per scenario at the timing check of Section 9), it is computed in the first 200 studies only, and the rate at which it triggers is reported.
- **Also evaluated:**
  - the order-statistic bounds, the 2nd-smallest and 2nd-largest subject estimates;
  - any revised estimator from Section 3.5.
- **The checks** enter Layer C only if their classification is kept as an inferential claim:
  - 200 studies, with a Monte Carlo error of 0.011 at a 2.5% tail;
  - bootstrap recurrence excluded.

### 3.5 Methods evaluated
- **The published methods:**
  - the R_h^−1/2 and 1/R_h extrapolants;
  - the "smaller" and "larger" combinations;
  - L, U and the guard.
- **Candidates for revision,** chosen on development simulations only:
  - (a) the R_h^−1/2 extrapolant alone, with 1/R_h kept as a sensitivity;
  - (b) more reference sizes from the same block H: R_h ∈ {8, 10, 20, 40} in EEG (adding 4 sub-blocks of 10), {6, 8, 12, 24} in mouse;
  - (c) cross-validated reference distances: inner products between the two halves of each sub-block, which are unbiased for the squared distance.
- **Selection rule.** Among the candidates with |median bias| ≤ 0.03 in every unsaturated development cell, choose the smallest RMSE; on a tie, the simplest (a). If none qualifies, the published pair is kept and reported as model-dependent.

### 3.6 Performance measures

Each measure is reported with its Monte Carlo standard error, overall and by status (Section 4).
- **Point estimates:**
  - bias (mean and median);
  - mean absolute error;
  - RMSE;
  - 90th percentile of |error|;
  - the share of estimates outside [0, 1];
  - the share not evaluable.
- **Bounds:**
  - Pr(L > Θ) and Pr(U < Θ), against 0.025;
  - coverage of [L, U], against 0.95;
  - width;
  - the rate of informative declarations (L above, or U below, the thresholds 0.5 and 0.25);
  - the rate of erroneous declarations;
  - all of these also against the sample-median target.
- **Screening triggers** (not calibration): |bias| > 0.05 or RMSE > 0.10.

### 3.7 Development pilot

- **Seeds:** separate seeds, at 20 repetitions per cell.
- **What it does:**
  - times the first 20 runs (Section 9);
  - sets r (T2) and α (T3);
  - checks that the exact θ across cells spans below 0.25, 0.25–0.5, 0.5–0.9 and above 0.9, with cells near 0.25 and 0.5. If not, it replaces λ values (for example 0.25 or 4) and records the change;
  - fixes the status thresholds (Section 4);
  - applies the selection rule of Section 3.5.
- **What it produces:** its choices and results go into a dated addendum to this protocol, committed before any validation seed is used.
- **Isolation:** development data are never reused for validation.

## 4. Estimation statuses (observable)

**Per subject and cell, from quantities an analyst can compute:**
- the pooled number of accepted neighborhoods;
- the deepest reference's disagreement d(R_h,max);
- the fall d(R_h,min) − d(R_h,max);
- the finite-depth error: the line through the two shallower sizes predicts the deepest, and the error is predicted minus observed;
- the divergence |FC_sqrt − FC_inv|.

**A cell's status:**
1. **Insufficient acceptance:** fewer than 5 subjects have at least A_min pooled accepted neighborhoods.
2. **Reference saturated:** over subjects, the median d(R_h,max) is ≥ c_sat and the median fall is ≤ f_sat.
3. **Extrapolation unresolved:** over subjects, the median |finite-depth error| is > c_fd, or the median divergence is > c_div.
4. **Estimable:** otherwise.

**Thresholds.**

| Threshold | Default |
|---|---|
| A_min | 20 |
| c_sat | 0.95 |
| f_sat | 0.02 |
| c_fd | 0.03 |
| c_div | 0.10 |

The development pilot replaces these defaults with the most permissive values for which estimable development cells meet |median bias| ≤ 0.03 and RMSE ≤ 0.08. The validation then reports performance overall and by status.

No status may depend on quantities an analyst cannot observe, such as true rank margins.

## 5. Reporting rules, decided before any result (Gate A)

For each real cell (dataset × map × budget × quantity), Stage 1 gives:
- **its status;**
- **the accuracy of the chosen estimator in the matching validation cells.** Those are cells with the same status whose median d(R_h,max) and fall are both within ±0.05 of the real cell's.
  - If no validation cell matches, the real cell is *unvalidated* and can be model-dependent at best.

**The four classes:**
- **Estimable, calibrated.** Threshold claims are allowed ("exceeds 50%"). All of these must hold:
  - the status is estimable;
  - in the matching cells, |median bias| ≤ 0.03 and RMSE ≤ 0.08;
  - for a threshold claim, Layer C gives Pr(L > Θ) ≤ 0.025 + 2 Monte Carlo errors, or the corresponding rule for U.
- **Model-dependent.** The estimate is reported and labeled as model-dependent; there are no threshold claims. All of these must hold:
  - the status is estimable;
  - |median bias| ≤ 0.05 and RMSE ≤ 0.10 in the matching cells;
  - the calibration criterion fails or was not run.
- **Finite-reference only.** Applies in every other case. The following are reported, never called false comfort:
  - d^(r) at the deepest reference, with its subject-bootstrap interval;
  - the disagreement curve over reference sizes;
  - the enrichment over random at that depth.
- **Stop (recommendation).** Neither latent nor finite-reference quantities distinguish the checks from a random selection, or the budgets from each other, in either dataset.

The author chooses the paper's direction at Gate A after seeing these classifications. A different direction for each dataset is allowed.

## 6. Revision analyses on the real data (all descriptive)

1. **Trial manifest** (`make_trial_manifest.py` → `results/manifest/`).
   - **EEG:**
     - subject, image and stored position;
     - the session from the replayed permutation;
     - the merged index;
     - the block in each confirmatory split.
     - Whether the merged index gives the chronological order within a session is established from the authors' public preprocessing code and stated either way. No data are downloaded.
   - **Mouse V1:**
     - experiment, image and presentation rank;
     - onset frame, epoch (1–3) and stratum;
     - the block in each split.
2. **EEG blocks from separate sessions** (`run_session_disjoint.py`). This is a sensitivity analysis, not a measurement of leakage.
   - **Design.** For each of the 12 ordered pairs of sessions (a, b):
     - F = 8 random repetitions from session a;
     - U = 16 from session b, with R ∈ {8, 16} its first R;
     - H = the 40 repetitions of the other two sessions, balanced over them.
   - **Analysis:**
     - the primary map with three seeds;
     - references at {8, 20, 40};
     - statuses, estimates, finite-reference disagreement and enrichment.
   - **Control.** 12 random splits by the published splitter with the same block sizes (2, 4 and 10 per session).
   - **Report.** The paired difference per subject, and the session mixture behind each target.
   - **What it cannot show.** It does not bear on R = 32.
3. **Mouse overlap accounting.**
   - **Source:** the downloaded NWB files, checked against their SHA-256 digests.
   - **What is counted:**
     - per experiment, the overlap in frames between consecutive response windows (0.5 s and 0.25 s);
     - per confirmatory split, the share of consecutive presentation pairs that fall in U and H;
     - the expected count per image pair.
4. **Per-point rerun** (`run_points.py`).
   - **What it re-creates,** with the committed seeds, for every budget:
     - the primary and cv-MDS maps;
     - every check's labels.
   - **Reproduction.** It asserts that every committed count is reproduced exactly.
   - **What it stores, per point:**
     - the naive reach;
     - each check's score (per-point trustworthiness; the scDEED score and its null's 95th percentile; the cross-validated reach; the number of seed recurrences; the bootstrap support; the split-half recurrences);
     - the reference reach and the continuous overlap per sub-block.
5. **Benchmark** (`analyze_benchmark.py`). For every check and budget, it reports:
   - the native operating point;
   - the random comparator matched per map, and enrichment;
   - yield, accepted counts and contributing subjects;
   - the total repetitions used (F + U + H), beside R.
   - **Matched-acceptance curves** at 5, 10, 20, 30 and 50%, for each check with a continuous score:

     | Check | Accepts first |
     |---|---|
     | plug-in | the smallest naive reach |
     | trustworthiness | the highest τ_i |
     | scDEED | the highest score |
     | cross-validated distances | the smallest cross-validated reach |
     | bootstrap recurrence | the highest support |

     Ties are included. Seed and split-half recurrence keep their native points only.
   - **Status.** Each quantity carries its Section 4 status. Latent values are reported only where the Section 5 rules allow; finite-reference values are reported everywhere.
6. **scDEED.**
   - The registered null is documented: one permutation per map, with the null pooled over that map's n scores.
   - As a sensitivity, the null that permutes the 20 representation coordinates, with the same seeds and one permutation per map.
7. **cv-MDS diagnostics.** Per subject and budget:
   - the share of negative cross-validated squared distances;
   - the share of zeros in the MDS input √max(0, d̂²);
   - the Kruskal stress-1 of each MDS fit;
   - the ranking convention: signed d̂² for the naive labels, plug-in distances for the references.
8. **Statuses** for every real cell and check, by the rules of Section 4.

Stages 1 and 2 report in `docs/revision_results.md`.

## 7. Correctness tests

`tests/test_estimand.py`, with its log in `results/revision/test_estimand.log`, runs before any analysis of this protocol.
- **Zero noise:**
  - naive labels equal the truth labels;
  - every disagreement is 0;
  - both extrapolants return 0.
- **Exact references** (noise-free block H, noisy U): the disagreement equals the exact false comfort at every reference size, and both extrapolants return it.
- **Zero acceptance:** false comfort is NaN, not 0. Fewer than 5 subjects with a pass gives "not evaluable".
- **Full acceptance:** a check that passes everything has FC = 1 − TF.
- **Ties:** ties are broken by object index, deterministically.
- **Count identities:**
  - within each map: the plug-in check's joint counts equal the confirmatory counts; X_b ≤ n10; X_a ≤ X_pass·nb;
  - within each subject: acceptance = faithful yield + unfaithful yield;
  - these hold in the synthetic cases and in every committed row.
- **Weighting:**
  - the pooled disagreement does not change when every map has proportionally more sub-blocks;
  - nb is constant across maps at each reference size in the committed rows;
  - more splits with identical counts change no estimate.

## 8. Integrity

**Order of commits:**
1. This protocol, `analyze_review_audit.py` and the tests. This is Stop 1, on the author's approval.
2. The development pilot runs, on development seeds.
3. The addendum is committed, fixing r, α, the λ grid, the status thresholds and the estimator.
4. The validation runs and the Stage 2 analyses.
5. The results are committed.
6. Gate A.

**Rules:**
- **No tuning on validation truths.** Development and validation seeds are disjoint, and the code asserts it. No validation seed is used before the addendum is committed.
- **Deviations** from this protocol are recorded in `docs/revision_results.md`.
- **The audit's numbers** motivate this protocol and are reported as computed before it.

## 9. Seeds, compute and stop conditions

**Seed ranges.** Seeds are tuples, as in `run_checks.py`. Every range is disjoint from those already used (3,090,000 to 5,700,000, including the display seeds 4,000,101–4,000,203):

| Range | Use |
|---|---|
| 6,000,000 | development pilot |
| 6,500,000 | Layers A and B |
| 7,000,000 | Layer C studies |
| 7,500,000 | the oracle θ_v |
| 8,100,000 | the session-disjoint analysis |
| 8,200,000 | the alternative scDEED null |

**Compute estimates.** A t-SNE map takes about 0.8 s on 200 objects and 0.5 s on 118, on one thread of the 36.

| Stage | Size | Time |
|---|---|---|
| Layer B | ≈140,000 maps | ≈1–1.5 h wall |
| Layer A | — | under an hour |
| Layer C | — | ≈5 h per scenario and map, ≈20 h for EEG C1 and C1-RSA under two scenarios |
| Stage 2 | — | ≈2 h, like the comparison of checks (53 and 21 min) |

**Timing check.** It is made on the first 20 development runs. If Layer C would exceed 8 h per scenario and map, it runs 500 studies instead of 1,000 (a Monte Carlo error of 0.007 at a 2.5% tail), and the change is recorded.

**No downloads.**

## 10. Out of scope

- New datasets.
- Other embedding methods (for example UMAP).
- Re-processing the raw EEG before the dataset's noise normalization (audit R20).
- Re-splitting the mouse timeline (audit R22).
- Per-point certification.
- Any change to the original claims' decision rules.
