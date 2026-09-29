# Confirmatory analysis: results

2026-09-29. Protocol: `docs/note_protocol.md`, pre-registered at commit `598562e`.

**How it was run.**
- At that commit, with the code unchanged.
- Every data file matched its SHA-256 in the manifest before loading.
- 466 s on 10 workers.
- No deviation from the protocol.

The printout is `results/things_confirmatory/confirmatory_analysis.txt`, and the numbers
below come from `confirmatory_summary.json`. Every estimate is a median over the 10
subjects. Pairs are given as (linear in R_h^(−1/2)) / (linear in 1/R_h). Brackets are 95%
subject-bootstrap intervals. L is the registered conservative lower bound.

## 1. Verdict (protocol Section 6)

| claim | estimate | L | threshold | result |
|---|---|---|---|---|
| **C1** false comfort, 8 trials, primary cell | 0.811 [0.642, 0.880] / 0.865 [0.741, 0.915] | 0.642 | 0.5 | **supported** |
| **C1-RSA** the same, MDS of cross-validated distances | 0.623 [0.550, 0.826] / 0.720 [0.669, 0.882] | 0.550 | 0.5 | **supported** |
| **C2** false comfort, 32 trials, primary cell | 0.476 [0.231, 0.576] / 0.635 [0.439, 0.710] | 0.231 | 0.25 | not supported |
| **C3** mislabel rate at 8 / 16 / 32 trials, primary cell | 0.169 / 0.208 / 0.210 (smaller extrapolant) | 0.159 / 0.180 / 0.158 | 0.10 each | **supported** |

- **Generality.** C1's criterion holds in 10 of the 12 cells, and the rule allows two
  exceptions, so the headline is stated as general.
- **Estimator check (Section 5).** The conservative estimator's median bias is −0.008 for
  false comfort and −0.003 for the mislabel rate, so the bias guard did not trigger.
- **What the note may say:** the headline holds for standard maps and for RSA-style maps.

## 2. The primary cell

t-SNE (perplexity 30), plug-in distances, whitened representation, k = 10, ε = 0.2,
k′ = 20:

| trials | naive faithful | truly faithful | mislabel rate | false comfort |
|---|---|---|---|---|
| 8 | 0.188 | 0.066 / 0.046 | 0.169 / 0.172 | 0.811 / 0.865 |
| 16 | 0.261 | 0.138 / 0.097 | 0.208 / 0.219 | 0.646 / 0.751 |
| 32 | 0.325 | 0.233 / 0.163 | 0.210 / 0.243 | 0.476 / 0.635 |

- **At 8 trials** the map calls about 19% of points faithful, while about 5–7% truly are.
  Over four in five of its faithful calls are wrong, which is 15–16% of all points.
- **At 32 trials** false comfort is still estimated at 0.48–0.64. The registered bound
  (0.231) falls just short of C2's 0.25, so the note cannot claim C2.
- **The pilot's exploratory figures replicate on new splits and seeds:** 84% at 8 trials
  and 47% at 32 then; 0.81–0.87 and 0.48–0.64 now.

## 3. The RSA-style pipeline

MDS of cross-validated, noise-normalized distances:

| trials | naive faithful | truly faithful | mislabel rate | false comfort |
|---|---|---|---|---|
| 8 | 0.020 | 0.035 / 0.025 | 0.039 / 0.033 | 0.623 / 0.720 |
| 16 | 0.058 | 0.066 / 0.048 | 0.067 / 0.064 | 0.600 / 0.706 |
| 32 | 0.095 | 0.087 / 0.063 | 0.076 / 0.081 | 0.445 / 0.598 |

Cross-validated distances do not remove false comfort. They make faithful calls rare: at 8
trials the map calls only 2% of points faithful, and most of those calls are still wrong.
The false-comfort share stays high, but the points it affects are few, about 1–1.5% of the
map. That is also why the mislabel rate is low here (3–8%). The note must report both
facts, not only the share.

## 4. Robustness at 8 trials (one factor at a time)

| cell | naive faithful | false comfort | L | C1 holds |
|---|---|---|---|---|
| primary | 0.188 | 0.811 / 0.865 | 0.642 | yes |
| MDS of cross-validated distances | 0.020 | 0.623 / 0.720 | 0.550 | yes |
| k = 5 | 0.313 | 0.884 / 0.917 | 0.753 | yes |
| k = 20 | 0.210 | 0.493 / 0.633 | 0.312 | **no** |
| ε = 0.1 | 0.067 | 0.918 / 0.943 | 0.798 | yes |
| ε = 0.3 | 0.363 | 0.687 / 0.770 | 0.526 | yes |
| t-SNE perplexity 15 | 0.150 | 0.825 / 0.874 | 0.665 | yes |
| t-SNE perplexity 50 | 0.169 | 0.774 / 0.837 | 0.624 | yes |
| PCA-2 | 0.041 | 0.512 / 0.637 | 0.437 | **no** |
| MDS, Euclidean | 0.035 | 0.730 / 0.804 | 0.582 | yes |
| early window 0.05–0.50 s | 0.192 | 0.785 / 0.845 | 0.611 | yes |
| unwhitened representation | 0.160 | 1.000 / 1.000 | 0.996 | yes |

- **The two failing cells** (k = 20 and PCA-2) have point estimates of 0.49–0.64. They fail
  only on the conservative bound.
- **Without noise whitening** the true faithful share is zero to three decimals at every
  budget, so every faithful call is wrong. The noise-normalized representation of the
  primary cell is the favorable case, not the extreme one.

## 5. What the note can and cannot say

**Can say.** For t-SNE maps of THINGS-EEG2 image responses averaged over 8 trials, most
neighborhoods that a plug-in check calls faithful are not faithful to the noise-free
responses:
- in t-SNE maps: an estimated 81–87%, with a conservative lower bound of 64%;
- in MDS maps of cross-validated distances: 62–72%, with a bound of 55%, although those
  maps make faithful calls for only 2% of points;
- the finding holds in 10 of 12 pre-specified variations.

**Cannot say.**
- That false comfort exceeds a quarter at 32 trials: the estimate is 0.48–0.64, but the
  registered bound is 0.23.
- Anything beyond this dataset. The data are the same as the pilot's, with new splits and
  seeds. That guards against flexible analysis, not against features of THINGS-EEG2.

## 6. Next

- **Independent confirmation** needs a second replicate-rich dataset in another modality,
  which needs a download and the owner's permission.
- **Then the note itself.** It would carry:
  - the estimator, presented as an application of classification accuracy against true
    scores with a SIMEX-type correction;
  - the citations required by `docs/reframe_prior_art.md`;
  - this analysis as its main result.
