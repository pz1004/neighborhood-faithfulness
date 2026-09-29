# Second dataset: results (Allen two-photon, mouse V1)

2026-09-29. Protocol: `docs/note_protocol_allen.md` (commit `985098d`), with code
pre-registered at `ea36f1d`.

**How it was run.**
- The analysis code was unchanged from that commit; later commits touched only the data
  manifest and `.gitignore`.
- Every recording matched its SHA-256 before loading.
- 158 s on 10 workers.
- No deviation from the protocol.

The printout is `results/allen_confirmatory/allen_analysis.txt`. The exploratory check of
Section 4 is `results/allen_confirmatory/explore_snr.txt`, from `explore_allen_snr.py`.

Every estimate is a median over 10 mice. Pairs are given as (linear in R_h^(−1/2)) /
(linear in 1/R_h). L is the registered conservative lower bound.

## 1. Verdict (protocol Section 5)

| claim | estimate | L | threshold | result |
|---|---|---|---|---|
| **C1** false comfort, 8 trials, primary cell | 1.000 / 1.000 | 0.997 | 0.5 | **supported** |
| **C1-RSA** the same, MDS of cross-validated distances | 1.000 / 1.000 | 1.000 | 0.5 | **supported** |
| **C3** mislabel rate at 8 / 16 trials | 0.161 / 0.172 | 0.154 / 0.158 | 0.10 each | **supported** |

- **Generality.** C1's criterion holds in all 12 cells.
- **Estimator check.** The conservative estimator's median bias is +0.003 for false comfort
  and +0.000 for the mislabel rate, so the bias guard did not trigger.
- **The pre-data concern did not bind.** The protocol noted that the bound would be strongly
  conservative under this design's shorter reference. The bounds still sit at 0.93–1.00.
- **Joint conclusion under the protocol:** the note's headline is independently confirmed
  in a second species and modality, for standard and RSA-style maps.

## 2. The primary cell and the RSA-style pipeline

| cell | trials | naive faithful | truly faithful | mislabel rate | false comfort |
|---|---|---|---|---|---|
| primary (t-SNE) | 8 | 0.161 | 0.000 / 0.000 | 0.161 / 0.161 | 1.000 / 1.000 |
| primary (t-SNE) | 16 | 0.172 | 0.002 / 0.002 | 0.172 / 0.172 | 0.994 / 0.996 |
| MDS of cross-validated distances | 8 | 0.009 | 0.000 / 0.000 | 0.010 / 0.010 | 1.000 / 1.000 |
| MDS of cross-validated distances | 16 | 0.022 | 0.001 / 0.001 | 0.024 / 0.023 | 1.000 / 1.000 |

- **What the t-SNE maps show.** With 8 or 16 trials, the maps show essentially no truly
  faithful neighborhood, while the plug-in check calls 16–17% of points faithful.
- **What cross-validated distances change.** They again make faithful calls rare (1–2%),
  but the calls that remain are no more accurate.

## 3. Robustness at 8 trials

The naive faithful share, the false comfort (both extrapolants) and L, for each cell:

| primary | cv-MDS | k = 5 | k = 20 | ε = 0.1 | ε = 0.3 | perp. 15 | perp. 50 | PCA-2 | MDS | 0–0.25 s | unwhitened |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.161 | 0.009 | 0.199 | 0.278 | 0.063 | 0.310 | 0.117 | 0.127 | 0.017 | 0.089 | 0.156 | 0.460 |
| 1.000 | 1.000 | 0.990–0.994 | 0.997–0.998 | 1.000 | 0.991–0.994 | 1.000 | 0.997–0.998 | 1.000 | 1.000 | 1.000 | 1.000 |
| 0.997 | 1.000 | 0.969 | 0.993 | 1.000 | 0.978 | 0.988 | 0.993 | 1.000 | 1.000 | 0.998 | 0.999 |

The true faithful share is at most 0.014 in every cell at both budgets. Without whitening,
the naive check calls 46% of points faithful and none are.

## 4. Why false comfort is total here (exploratory, post hoc)

The data are split 0 of each mouse, U + H, 40 trials, in the F-fitted 20-dimensional
space. THINGS-EEG2 is shown alongside, with 72 trials.

| | mouse V1 (Allen) | human EEG (THINGS-EEG2) |
|---|---|---|
| image pairs with a positive cross-validated distance | 0.937 (0.780–0.971) | 1.000 |
| split-half reliability of the distance matrix | 0.726 (0.469–0.906) | 0.879 |
| distance to the 8th true neighbor | 0.397 | 1.545 |
| noise between two 16-trial means | 1.653 | 1.683 |
| **ratio** | **0.24** | **0.92** |

- **The result is not vacuous.** The mouse data have reliable, image-specific structure:
  most pairs of images are distinguishable, and the distance matrix replicates across
  independent trials.
- **The neighbor distances are the problem.** A point's 8th neighbor is about a quarter of
  the distance noise between 16-trial means away. For 42 of the 118 images, even 40 trials
  cannot separate the 8th neighbor at all.
- **So displays at these budgets carry no true neighborhoods.** A plug-in check still
  reports 16% of them faithful.
- **The same mechanism operates in both datasets, and its severity differs.** In THINGS-EEG2
  the ratio is about 0.9 and 5–23% of neighborhoods are truly faithful. In mouse V1 it is
  0.24 and essentially none are.
- **Status of the ratio.** It is an exploratory account of the difference, not a tested
  claim.

## 5. What the note can say now

**Supported by pre-registered tests on two independent datasets:**
- In 2-D maps of trial-averaged responses to images, most neighborhoods that a plug-in
  check calls faithful are not faithful to the noise-free responses.
- This holds for human EEG (THINGS-EEG2: 81–87% at 8 trials) and for mouse V1 two-photon
  imaging (Allen Brain Observatory: essentially 100% at 8 and 16 trials).
- It also holds for RSA-style maps of cross-validated distances, where faithful calls are
  rare.
- It holds in 10 of 12 and 12 of 12 pre-specified variations respectively.

**Exploratory, to be labeled as such:** the severity tracks the ratio of neighbor distance
to the noise between trial-averaged means.

**Limits:**
- Two datasets, both of visual responses to images.
- One neighborhood definition: k = 10, ε = 0.2, k′ = 2k, with variations one factor at a
  time.
- The claim at 32 trials (C2) failed on THINGS-EEG2 and could not be tested here.
