# Confirmatory protocol: false comfort in 2-D maps of trial-averaged EEG

Version 1, 2026-09-29. Pre-registered: this protocol and the code that implements it were
committed before the analysis was run. Plan Amendment 3 points here.

## 1. Purpose, and what the data can and cannot confirm

The reframe's prior-art check (`docs/reframe_prior_art.md`, commit `97680d1`) found that
only the empirical finding is new. That finding: most neighborhoods that look faithful in
t-SNE maps of trial-averaged EEG are not faithful in truth. This protocol tests it with
estimators and thresholds fixed in advance, and answers the objections the check raised.

**The data are not new.** The analysis uses the same THINGS-EEG2 test set (10 subjects,
200 images × 80 repetitions) that the need pilot analyzed. The false-comfort figures to be
confirmed (84% at R = 8, 47% at R = 32) came from exploratory analyses of those data. The
new random splits and seeds guard against flexible analysis choices, but not against
features of this particular dataset. Independent confirmation needs a second dataset
(Section 8).

## 2. Quantities

For a display, a neighborhood size k, a tolerance ε, and k′ = 2k:
- **N_i** is the *naive label*: faithful if the reach computed from the user's own distances
  is at most k′.
- **T_i** is the *truth label*: the same, computed from the noise-free means.

The four quantities:
- **NF** = P(N = 1): the naive faithful share.
- **TF** = P(T = 1): the true faithful share.
- **q** = P(N ≠ T): the mislabel rate.
- **FC** = P(T = 0 | N = 1): *false comfort*, the share of the naive map's faithful calls
  that are wrong.

FC is primary. The other three are secondary.

## 3. Design

**Blocks and representation** follow plan Amendment 2:
- F has 2 repetitions per session and fits the linear map; U has 8 per session and serves
  as the user; H has 10 per session and serves as the reference.
- New splits: 10 per subject, with seeds 4_100_000 + 1000·s + split.
- User budgets R ∈ {8, 16, 32}, taken as session-stratified prefixes of U.
- Reference sub-blocks of size R_h ∈ {8, 20, 40}, disjoint and session-stratified (5, 2 and
  1 of them).

**Pipelines.**
- **P1 (plug-in, primary).** The display is fitted to the user's means with Euclidean
  distances. Naive labels use plug-in distances between the user's means.
- **P2 (cross-validated; the RSA objection).**
  - The user's R repetitions are split into two session-stratified halves A and B.
  - Cross-validated squared distances are ⟨c_i^A − c_j^A, c_i^B − c_j^B⟩, unbiased for the
    true squared distances.
  - The display is metric MDS of √max(0, d²_cv), the standard RSA display.
  - Naive labels use the cross-validated distances.
- **Reference labels, in both pipelines,** use plug-in distances between each reference
  sub-block's means. The extrapolation in Section 4 removes the reference's noise and its
  bias.

**Cells.** The primary cell is P1, the whitened representation, samples 0–0.79 s, t-SNE with
perplexity 30, k = 10 and ε = 0.2. Every other cell changes one factor:

| # | cell | change from the primary cell |
|---|---|---|
| 1 | primary | none |
| 2 | rsa_cv_mds | pipeline P2: MDS of cross-validated distances |
| 3, 4 | k5, k20 | k = 5 or k = 20 (k′ = 2k) |
| 5, 6 | eps0.1, eps0.3 | ε = 0.1 or ε = 0.3 |
| 7, 8 | tsne_perp15, tsne_perp50 | t-SNE perplexity 15 or 50 |
| 9 | pca2 | PCA-2 display |
| 10 | mds | metric MDS of Euclidean distances |
| 11 | early_window | samples 0.05–0.50 s |
| 12 | unwhitened | PCA-20 of F's means, without whitening |

**Display seeds.**
- t-SNE (PCA initialization): seeds 4_000_101, 4_000_102 and 4_000_103.
- MDS (one initialization each): seeds 4_000_201, 4_000_202 and 4_000_203.
- PCA-2 is deterministic.

**Integrity.** Every data file is checked against its SHA-256 in
`data/things_eeg2/manifest.json` before it is loaded.

## 4. Estimation

**Disagreement statistics.** For each subject, cell, user budget R and reference size R_h,
these are pooled over splits, display seeds and reference sub-blocks:
- d_q(R_h) = the mean of 1{N ≠ T_ref};
- d_FC(R_h) = Σ #(N = 1, T_ref = 0) / Σ #(N = 1), a ratio of sums within the subject;
- d_TF(R_h) = the mean of 1{T_ref = 1}.

**Extrapolation.** Each statistic is extrapolated to an infinite reference with two
extrapolants: linear in R_h^(−1/2), and linear in 1/R_h. Each gives a per-subject estimate.
The estimate reported is the median over subjects.

**Intervals.** 95% intervals come from 10,000 bootstrap resamples of subjects (seed
4_300_000).

**Conservative lower bound.** For a quantity, L is the 2.5th percentile over resamples of
the smaller of the two extrapolants' medians.

## 5. Estimator check (semi-synthetic, before the claims are read)

- **Setting.** Per subject, on the first split: the truth is H's means; user blocks and
  reference sub-blocks are drawn independently by resampling U's residuals within each
  image.
- **Scope.** The primary cell only, with 20 draws per subject.
- **Output.** Exact FC and q, next to both estimators.
- **Bias.** For the conservative combination: the median over subjects of the smaller
  estimate minus the exact value.

## 6. Claims and decision rules

**Claims, all in the primary cell unless stated:**
- **C1 (headline).** FC(8) > 0.5. Supported if L > 0.5.
- **C1-RSA.** FC(8) > 0.5 in cell 2. Supported if L > 0.5.
- **C2.** FC(32) > 0.25. Supported if L > 0.25.
- **C3.** q(R) > 0.10 for every R ∈ {8, 16, 32}. Supported if L > 0.10 at each R.

**Generality.** The headline is stated as general if C1's criterion holds in all 12 cells
except at most two. Otherwise it is stated only for the cells where it holds, and those
are listed.

**Bias guard.** Suppose the check of Section 5 shows the conservative estimator
overestimating FC(8) by more than 0.05. Then C1, C1-RSA and the generality count are
re-evaluated with L reduced by that bias.

**What the note may say:**
- If C1 and C1-RSA are both supported: the headline holds for standard maps and for
  RSA-style maps.
- If C1 is supported but C1-RSA is not: the headline is restricted to maps built from
  plug-in distances, and the RSA result is reported as its limit.
- If C1 is not supported: there is no headline. The estimates are reported, and the
  recommendation is to stop.

## 7. Integrity

- **Code.** The code at the pre-registration commit is the code that runs. After a failed
  run, only errors that are not analysis choices may be fixed, and each fix is logged.
- **Changes.** Any change to an analysis choice is a dated amendment, reported next to the
  original.
- **Exploratory work** is labeled as such.

## 8. Out of scope here

A second, independent replicate-rich dataset in another modality. It needs a download and
the owner's permission.
