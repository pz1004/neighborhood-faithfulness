# Protocol: do common checks catch false comfort?

Version 1, 2026-10-01. Pre-registered: this protocol, the code that implements it and its
synthetic tests are committed together, before the analysis reads any data file. Plan
Amendment 5 points here.

## 1. Purpose, and what is known before the run

The confirmatory analyses (`docs/note_confirmatory_results.md`, `docs/note_allen_results.md`)
found that most neighborhoods that pass the plug-in check are not faithful to the
noise-free responses: 81–86% in THINGS-EEG2 and essentially all in mouse V1, at 8 trials.
An analyst has other ways to decide which neighborhoods of a map to trust. This protocol
asks whether any of the common ones would have caught the false comfort. Every check is
computed from the analyst's own trials, as an analyst would compute it. This is the
comparison with existing methods that a methods journal expects.

**What is known before the run.**
- For the plug-in check (N) in the primary cell: everything in the two results files.
- For the six other checks: nothing. None has been computed on either dataset.
- The maps are the confirmatory maps, re-created with the same splits, seeds and code
  (Section 3), so N's results on exactly these maps are known.
- The synthetic tests of the code (`results/checks/test_checks.log`) ran on synthetic data
  only.

## 2. The checks

Notation:
- E_i is the set of the k = 10 displayed neighbors of object i, with m = ⌈(1 − ε)k⌉ = 8
  and k′ = 20.
- *E_i recurs in a map M* means the reach of E_i computed from M's 2-D distances is at most
  k′: at least 8 of the 10 displayed neighbors are among i's 20 nearest objects in M.

Each check gives every object a pass or a fail. A pass means the check would let the
analyst read that neighborhood.

| check | family | the neighborhood passes if |
|---|---|---|
| **N**, plug-in | score against the averages | the reach of E_i computed from the analyst's averages is at most k′ (the confirmatory check) |
| **TW**, trustworthiness | score against the averages | its per-point trustworthiness is among the n_N highest in the map, ties included, where n_N is the number of points N passes in that map |
| **SCD**, scDEED | score against the averages, with a null | its scDEED reliability score reaches the 95th percentile of the permutation null ("trustworthy") |
| **CV**, cross-validated | distances | the reach of E_i computed from cross-validated squared distances is at most k′ |
| **SEED**, seed recurrence | stability | E_i recurs in the maps of the same averages made with the other two display seeds |
| **BOOT**, bootstrap recurrence | stability | E_i recurs in at least 95% of 50 maps of bootstrap resamples of the trials |
| **SH**, split-half recurrence | stability | E_i recurs in the map of each half of the trials |

**Details.**
- **TW.** T_i = 1 − 2 / (k(2n − 3k − 1)) · Σ_{j ∈ E_i} max(0, r_ij − k), where r_ij is the
  rank of j among i's neighbors in the averages (Venna & Kaski). The mean of T_i over
  points is the usual trustworthiness score. Matching TW's pass count to N's compares the
  two as rankings at the same number of neighborhoods read.
- **SCD** follows Xia, Lee & Li (2024, *Nat. Commun.* 15:1753), Methods, steps 1–6, with
  its default neighborhood x = 50%:
  - *score:* the Pearson correlation of two vectors of 2-D distances from i. One runs to
    i's round(n/2) nearest objects in the averages, ordered by their distance in the
    averages. The other runs to i's round(n/2) nearest objects in the map, ordered by their
    2-D distance.
  - *null:* scDEED permutes the measured features and then applies the pipeline's PCA step.
    Here, every measured feature of the analyst's averages (EEG channel-time samples; mouse
    cells) is permuted independently across objects. The permuted averages then go through
    the same fixed representation fitted on block F. They are mapped with the display's
    t-SNE settings and seed, and the n scores of that map form the null.
  - *why at the level of measured features:* permuting the 20 coordinates of the
    representation instead would keep the signal's unequal spread across coordinates. In
    the synthetic tests that null resembled the data so closely that scDEED called 11% of
    points trustworthy on strongly structured data, against 74% with the measured-feature
    null.
- **CV** uses the cross-validated squared distances of the RSA pipeline
  (`docs/note_protocol.md`, Section 3), with the same two halves as SH.
- **SEED, BOOT and SH** fit t-SNE with the display's settings (perplexity 30, PCA
  initialization) and seed.
  - SEED uses the other two of the three display seeds.
  - BOOT draws 50 resamples. In each, every object's R trials are drawn with replacement,
    independently across objects. It is the neighborhood-level analogue of bootstrap
    stability tools such as DynamicViz (Sun, Ma & Zou 2023), with trials resampled
    instead of observations.
  - SH uses the halves of the RSA pipeline, balanced over sessions (EEG) or strata
    (mouse), R/2 trials each. It is the check the note's Discussion suggested without
    testing.

**Why these.**
- *Plug-in scores against the averages* are what quality criteria for dimensionality
  reduction compute (Lee & Verleysen 2009; Venna et al. 2010). scDEED is the published
  per-point tool closest to the question; the Gate 1 audit named it the competitor whose
  labels should be measured against the truth.
- *Cross-validated distances* are the RSA recommendation.
- *Recurrence across seeds, resamples and halves* is the usual stability check.

## 3. Data, maps and references

- **Data, blocks and representation:** as in the confirmatory analyses. The splits reuse
  the confirmatory seeds: 4_100_000 + 1000·s + split (EEG) and 4_500_000 + 1000·s + split
  (mouse). The maps are therefore the confirmatory maps.
- **The primary cell only:** pipeline P1, the whitened representation, the full response
  window, t-SNE with perplexity 30, PCA initialization and the three display seeds, and
  k = 10, ε = 0.2, k′ = 20.
- **Budgets:** R ∈ {8, 16, 32} (EEG) and {8, 16} (mouse). The claims are at R = 8.
- **References:** sub-blocks of block H, as before: R_h ∈ {8, 20, 40} (EEG) and {8, 12, 24}
  (mouse).
- **The checks' own randomness** comes from seed sequences, where j indexes the display
  seed:
  - bootstrap resamples: (5_100_000, s, split, R, j) in EEG and (5_500_000, s, split, R, j)
    in mouse;
  - permutations: (5_110_000, …) and (5_510_000, …).
- **Reproduction check.** The plug-in counts of every re-created map are compared with the
  committed confirmatory rows. Any difference is reported, and the analysis proceeds on the
  re-created maps.

## 4. Quantities and estimation

For a check X and the truth label T:
- **PASS_X** = P(X passes). It is observed.
- **FC_X** = P(T = 0 | X passes): the false comfort of X. Primary.
- **CATCH_X** = P(X fails | N passes, T = 0): the share of the plug-in check's
  false-comfort neighborhoods that X flags. Secondary.
- **YIELD_X** = P(X passes, T = 1): the share of all points whose neighborhood passes X and
  is faithful. Secondary.
- **FC_NX** = P(T = 0 | N and X pass): the false comfort of X used on top of the plug-in
  check. Secondary.

**Estimation** extends Section 4 of `docs/note_protocol.md` from one label to joint ones.
- Per subject, pooled over splits, display seeds and reference sub-blocks, three
  reference rates are computed:
  - a_X(R_h) = P̂(X passes, T_ref = 0);
  - b_X(R_h) = P̂(N passes, X fails, T_ref = 0);
  - c(R_h) = P̂(N passes, T_ref = 0).
- Each rate is extrapolated with the two extrapolants, linear in R_h^(−1/2) and in 1/R_h.
- Then, per extrapolant: FC_X = â_X / PASS_X, CATCH_X = b̂_X / ĉ, YIELD_X = PASS_X − â_X,
  and FC_NX = (ĉ − b̂_X) / P(N and X pass).
- For N, FC_N reproduces the confirmatory estimator exactly.
- Estimates are medians over subjects, with 10,000 subject-bootstrap resamples (seed
  5_300_000 for EEG, 5_700_000 for mouse).
- **L_X** is the 2.5th percentile over resamples of the smaller extrapolant's median.
  **U_X** is the 97.5th percentile of the larger one's.
- **Subjects with no passes.** A subject whose check passes no neighborhood at a budget
  has no FC_X there and is left out of that median. If fewer than 5 subjects remain, FC_X
  is not evaluable at that budget.
- **Descriptive quantities.** CATCH_X is a ratio of two extrapolated rates. It can fall
  outside [0, 1] when the plug-in check's false comfort is rare; the synthetic tests show
  this. CATCH_X, YIELD_X and FC_NX are reported without a decision rule.

## 5. Estimator check (semi-synthetic, before the claims are read)

This extends `docs/note_protocol.md`, Section 5, to the checks.
- **Setting.** Split 0. The truth is H's means. The analyst's trials and the reference
  sub-blocks are drawn by resampling U's residuals within each image. The resampled
  residuals are the same in the measured features and in the representation, so scDEED's
  null sees the analyst's own averages.
- **Draws.** R = 8, 10 draws per subject. Each draw has the three display seeds and every
  check.
- **Seeds.** Draws come from (5_200_000, s, d) in EEG and (5_600_000, s, d) in mouse. The
  checks' randomness comes from (5_210_000 or 5_220_000, s, d, j) and (5_610_000 or
  5_620_000, s, d, j).
- **Output.** The exact FC_X next to both extrapolants, per subject.
- **Bias.** The median over subjects of (smaller estimate − exact) and of (larger estimate −
  exact).

**Bias guard.**
- If the smaller estimate overstates FC_X by more than 0.05 (median bias > +0.05), L_X is
  lowered by that bias.
- If the larger estimate understates FC_X by more than 0.05 (median bias < −0.05), U_X is
  raised by that amount.

## 6. Classification and decision rules

For each dataset and each check X other than N, at R = 8:
- **leaves most false comfort** if L_X > 0.5;
- **removes most false comfort** if U_X < 0.5;
- **undetermined** otherwise;
- **not evaluable** as defined in Section 4.

**What the manuscript may say:**
- **A check that leaves most false comfort in a dataset:** most neighborhoods that pass it
  there are not faithful.
- **A check that removes most false comfort in a dataset:** it is reported as a partial
  remedy there, with its pass share and yield. Nothing is claimed for it in the other
  dataset unless it is classified the same way there.
- **"None of the common checks removed most false comfort"** may be said only if no check
  removes most false comfort in either dataset.
- **"All six checks left most false comfort in both datasets"** may be said only if all 12
  classifications are "leaves most false comfort".
- **Everything else is descriptive:** CATCH_X, YIELD_X, FC_NX, and every budget above 8.

## 7. Expectation, recorded before the run

This is a prediction, not a decision rule.
- **The argument.** Under Gaussian trial noise with a shared covariance, the analyst's
  averages carry all the information their trials hold about the noise-free means. The
  within-image residuals add information only about the size of the noise.
  - Split halves are the averages plus and minus a noise-only difference: c ± d, where d is
    independent of the averages.
  - Bootstrap resamples are the averages plus resampled noise.
  - Cross-validated distances are the plug-in squared distances minus a noise-only term.
  - Every check here is therefore a function of the averages, the noise level and
    randomization. It can be stricter than the plug-in check, but it cannot see where the
    averages are wrong.
- **The expectation.** The certificates of the Gate 3 pilot used the averages and the
  noise level rigorously and resolved no neighborhood even at 72 trials. So all six checks
  are expected to leave most false comfort in both datasets.
- **In mouse V1**, where essentially no neighborhood is faithful, any passing neighborhood
  is expected to be wrong. The informative quantity there is how few neighborhoods a check
  passes.
- **Stability is not faithfulness.** The stability checks (SEED, BOOT, SH) are expected to
  pass some neighborhoods that N fails, because a map's distortions can recur across seeds
  and resamples.

## 8. Integrity

As in `docs/note_protocol.md`, Section 7:
- the code at the pre-registration commit is the code that runs;
- after a failed run, only errors that are not analysis choices may be fixed, and each fix
  is logged;
- any other change is a dated amendment, reported next to the original;
- exploratory work is labeled as such.

## 9. Compute

About 65,000 t-SNE fits in EEG and 49,000 in mouse, one thread each. On 36 cores this is
roughly 25 and 10 minutes.
