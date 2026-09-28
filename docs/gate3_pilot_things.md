# Gate 3 need pilot: THINGS-EEG2

2026-09-29. Protocol: plan Amendment 2 (commit `5024fe6`, made before any data value was
loaded). Raw results are in `results/things_pilot/`, and the decision printout is
`results/things_pilot/pilot_analysis.txt` from `analyze_things_need_pilot.py`. Section 5 is
exploratory and post hoc; its printout is `results/things_pilot/explore_snr.txt`, from
`explore_things_snr.py`.

## 1. Outcome

The outcome is the protocol's third branch: need passes, use fails.

| criterion (medians over 10 subjects) | result | verdict |
|---|---|---|
| **need**: naive mislabel rate at k′ = 2k, ε = 0.2, t-SNE | q(8) = 0.172 [0.159, 0.197]; q(16) = 0.194; q(32) = 0.195; q(80), extrapolated, = 0.215 [0.166, 0.267] | **passes**, at R = 8 and at native R |
| **use**: E-pair-t resolved share at k′ = 20, R ≤ 16 | 0.000 at R = 8, 16, 32 and 72 | **fails** |
| **validity**: E-pair-t failure rate under resampled real noise | 0/300 draws for every subject | holds |

The pre-registered recommendation is therefore that the need–power tension holds on real
data, and that the method must be reframed before any further data.

On this dataset the naive neighborhood map is wrong for about one point in five, at every
number of repetitions measured. The certificate cannot identify a single faithful or
distorted point, even with 72 repetitions per image. Section 5 shows why: the obstacle is
signal-to-noise, not the multiplicity correction.

## 2. Data

- **Source.** Gifford, Dwivedi, Roig, Cichy (2022), THINGS-EEG2, preprocessed EEG, from
  OSF node `anp5v`.
- **What was fetched.** Only `preprocessed_eeg_test.npy` for each of the 10 subjects, by
  HTTP range requests into the per-subject zips (1.07 GB each); the training files were not
  downloaded. That is 208.6 MB transferred and 217.6 MB on disk per subject, 2.18 GB in
  total.
- **Verification.** Each file's CRC-32 was checked against the zip. The SHA-256s are in
  `data/things_eeg2/manifest.json`; the data files themselves are not tracked.
- **Loading.** The files are pickled dicts. They were loaded through a restricted
  unpickler that admits only NumPy array reconstruction; `tests/test_things.py` checks that
  it blocks a hostile pickle.
- **Format, as the authors' code documents.**
  - 200 images × 80 repetitions × 17 channels × 100 samples, from −0.2 to 0.79 s.
  - Channels whitened per session (MVNN).
  - 4 sessions of 20 repetitions, shuffled across sessions with seed 20200220. Replaying
    that shuffle gives the session of every repetition.
  - Consistent with the reconstruction, the session main effect is large (median
    F = 27.5). Session-stratified blocks cancel it.
- **One property of the preprocessing to keep in mind.** The authors' whitening matrix was
  estimated from all repetitions, so it depends weakly on the data it whitens. This pilot's
  own representation is fitted on block F only.

## 3. Results

### 3.1 Need

**Primary estimator.** Reference-extrapolated naive mislabel rate at k′ = 20. Disagreement
with plug-in labels from H sub-blocks of size R_h:

| user R | d(R_h = 8) | d(20) | d(40) | intercept q(R) | per-subject range |
|---|---|---|---|---|---|
| 8 | 0.174 | 0.173 | 0.173 | **0.172** | 0.131–0.206 |
| 16 | 0.254 | 0.237 | 0.226 | **0.194** | 0.166–0.233 |
| 32 | 0.315 | 0.267 | 0.243 | **0.195** | 0.152–0.255 |

- **Stability.** At R = 8 the disagreement does not depend on the reference's size, so the
  estimate needs no extrapolation.
- **The semi-synthetic check agrees.** With exact truth (H means) and resampled real noise,
  the naive mislabel rate is 0.161 (R = 8), 0.182 (16) and 0.188 (32).
- **The certified bounds are uninformative,** as the protocol allowed. E-pair-t on H
  (R = 40) resolves no point, so the bounds are [0, 1].
- **PCA-2 displays** (secondary) are distorted almost everywhere. The naive map calls only
  4–9% of points faithful, and its disagreement with the reference is 4–7%.

### 3.2 Use and validity

**No resolution at any budget.** E-pair-t's resolved share on the real user blocks is 0.000
at k′ = 20 for R = 8, 16, 32 and 72, and at most 0.0001 at k′ = 40. Meanwhile the naive map
calls 17.5% (R = 8) to 38.1% (R = 72) of points faithful.

**Validity.** In the semi-synthetic setting E-pair-t never made a false statement: 3000
draws, all pairwise intervals, and reach under two displays at three tolerances. The noise
here is resampled from the data itself.

### 3.3 Diagnostics (G3-diag)

Split 0, residuals of U + H in the 20-dimensional space; medians over subjects with ranges:

| diagnostic | median | range | reading |
|---|---|---|---|
| mean residual variance | 1.13 | 0.88–1.25 | F's whitening carries over to within about 25% |
| eigenvalue spread, max/min | 1.72 | 1.50–2.62 | mild anisotropy, well inside E-pair-t's robust range (κ ≤ 4) |
| eigenvalue spread, max/mean | 1.32 | 1.23–1.97 | |
| excess kurtosis, median coordinate | 0.35 | 0.11–1.11 | mildly heavy-tailed |
| excess kurtosis, worst coordinate | 1.02 | 0.20–13.3 | one subject has an outlying coordinate |
| dispersion of s_i² relative to homoscedastic | 1.63 | 1.21–3.45 | moderately heteroscedastic, within M_het |
| session × image interaction, F | 1.14 | 1.05–1.34 | small image-specific session effects |

Nothing here explains the lack of resolution. The noise is close to the model E-pair-t is
valid under.

## 4. Deviations and flaws, disclosed

- **A smoke test was run first.** Subject 1, one split, identical code. Its outputs were
  overwritten by the full run. No code or parameter changed after it.
- **The pre-registered "need fails" branch rested on a false premise.** It assumed that q
  falls with R. It does not: q rises from 0.172 at R = 8 to 0.195 at R = 32. The display
  itself improves with R, which brings more points near the tolerance boundary; the true
  faithful share grows from 0.05 to 0.20 (Section 5).
  - The flaw did not bind: the pass branch applied, with q(32) and q(80) both ≥ 0.05.
  - The extrapolation to R = 80 fits that rising trend. The robust statement is that q is at
    least 0.13 in every subject at every measured R.
- **Section 5 is exploratory,** chosen after the decision was computed.

## 5. Why nothing resolves (exploratory)

**Signal-to-noise** (`explore_snr.txt`: split 0, U + H, cross-validated distances from two
independent 36-repetition halves):
- **The images are distinguishable.** Every pair's cross-validated squared distance is
  positive. The median true distance is 2.47 whitened units, and the typical 8th-neighbor
  distance is 1.49.
- **The gap is small.** The certificate must tell the 8th true neighbor from the 20th (k = 10,
  ε = 0.2, k′ = 20). Those distances differ by only **16%**, about 0.23 units.
- **The intervals are wide.** At R = 72, E-pair-t's median half-width is **1.10**, about five
  times that gap.
- **Multiplicity is not the cause.** Unadjusted single-point 95% intervals would be about
  ±0.35, still three times too wide. Resolving 16% gaps would take thousands of
  repetitions per image. Relaxing the error control, for example to FDR, cannot close a gap
  of this size.

**False comfort** (semi-synthetic, t-SNE, k′ = 20, medians over subjects):

| R | truly faithful | naive map calls faithful | of those calls, wrong |
|---|---|---|---|
| 8 | 0.052 | 0.167 | **84%** |
| 16 | 0.126 | 0.217 | 63% |
| 32 | 0.200 | 0.273 | 47% |

The naive map's errors are mostly false comfort. It reads noise-driven neighborhoods as
faithful.

## 6. What this does to the project

**The two ends of signal-to-noise.** Together with Gate 2, the pilot covers both ends:

| | digits, injected noise (Gate 2) | THINGS-EEG2, real noise (this pilot) |
|---|---|---|
| naive map mislabels | 1.3–4.1% | 17–20% |
| naive "faithful" calls that are wrong | 1–3% | 47–84% |
| certificate resolves | 33–64% | 0% |

- Per-point certification is informative only where the naive map is already nearly right.
  Where the naive map misleads, it certifies nothing.
- This is review finding 3, the need–power tension, measured at both ends. On EEG the
  obstacle is not the simultaneous error control, so no relaxation of it rescues per-point
  statements.

**Recommendation.**
- *Pre-registered:* reframe before any further data.
- *In my reading:* stop the per-point certificate. What the pilot built instead is a
  dataset-level measurement: how wrong a display's neighborhoods are, and how much of its
  apparent fidelity is false comfort, estimated from replicate splits with subject-level
  intervals.
- *If continued,* that is the reframing: a neighborhood noise ceiling for 2-D displays,
  analogous to the noise ceilings of representational similarity analysis. It needs its
  own prior-art check before any work, covering RSA noise ceilings, split-half reliability
  of embeddings, and Lin & Fukuyama 2024.
