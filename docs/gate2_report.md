# Gate 2: development validity, robustness, power and need

2026-09-28. Plan version 1 (commit `ce69f4e`) with Amendment 1. Code at commit `f53aba7`.
The raw results are in `results/gate2/`, and every number below is printed by
`analyze_gate2.py` into `results/gate2/gate2_analysis.txt`. Two extra tabulations, in
Sections 5 and 6, have their own scripts and printouts in the same directory.

**Word use.** "Exact" in the plan means *valid in finite samples, not asymptotically*. It
does not mean tight. The Gate 1 audit rightly asked for "valid" instead, and this report
uses that word.

## 1. Decision

**Gate 2 passes.**
- **G2-V (validity).** It passes under Amendment 1 and fails under the original wording,
  for E-ball only (Section 3).
- **Default engine.** E-pair-t, the studentized pairwise engine. It is the one engine with a
  finite-sample guarantee that survives the robustness set (Section 3.2).
- **G2-P (power).** E-pair-t clears the floors: 0.326 at R = 8 (floor 0.20) and 0.504 at
  R = 16 (floor 0.30).

**Two results matter more than the gate.**
1. E-pair-t beats the ceiling of every ball-based certificate on the same panels, including
   on the predecessor's own target (Section 5). That narrows the predecessor's conclusions.
2. On this development panel the naive map is already right for 96–99% of points at
   k′ = 2k. The certificate never corrects it: it abstains on every point the naive map
   gets wrong (Section 6). Need is therefore the main risk, and Gate 3's need criterion
   decides the project.

## 2. Before Gate 2 proper

- **Structural tests.** `tests/test_reach.py` (plan Section 5.7) passes with 0 failures
  (`results/gate2/test_reach.log`).
- **Pilot reproduction.** The code reproduces every cell of the 2026-09-28 pilot exactly:
  16 synthetic cells and the digits t-SNE cells.
- **One new number seen early, as the plan allowed.** At the pilot's single configuration
  (digits, R = 8, t-SNE seed 42), E-pair-t certified 0.261 faithful and 0.066 distorted at
  k′ = 40. The Gate 2 mean for that configuration, over 15 displays, is 0.259 and 0.067.
- **Test seed.** The tests use 3_090_000, which is outside every Gate 2 base.

## 3. Validity (7900 draws, 24 configurations)

A draw is a failure for an engine if any statement the engine makes is false:
- every pairwise interval;
- every ball, for E-ball;
- every reach interval, at three tolerances and under two displays;
- for E-cmp, every elementary comparison for the displayed members, and all triples on
  digits200.

| configuration (R = 8 unless stated) | E-ball | E-pair-t | E-pair-pool | E-cmp |
|---|---|---|---|---|
| digits, gauss (1000 draws) | 57 | 6 | 11 | 0 |
| digits, gauss, R = 16 | 13 | 5 | 5 | 0 |
| digits, gauss_het / R = 16 | 10 / 12 | 2 / 4 | 299 / 300 | 0 / 0 |
| digits, t5 / R = 16 | 13 / 12 | 0 / 3 | 63 / 30 | 0 / 0 |
| digits, t3 / R = 16 | 15 / 10 | 6 / 9 | 293 / 285 | 0 / 0 |
| digits, rep_t5 / R = 16 | 6 / 7 | 0 / 4 | 230 / 154 | 0 / 0 |
| digits, aniso4 / R = 16 | 66 / 65 | 6 / 8 | 13 / 11 | 0 / 0 |
| digits, aniso16 / R = 16 | 241 / 246 | 38 / 63 | 74 / 83 | 0 / 0 |
| surface: gauss, gauss_het, t5, aniso4 | 14, 15, 16, 68 | 2, 1, 2, 2 | 5, 291, 32, 2 | 0 |
| digits, γ = 0.5: gauss, t5 | 13, 15 | 0, 1 | 0, 23 | 0 |
| digits200: gauss, t5, aniso4, rep_t5 | 13, 10, 40, 5 | 4, 6, 9, 1 | 7, 23, 5, 92 | 0 |

Every cell is failing draws out of 300 unless stated.

**Reach claims were never wrong.**
- No reach interval of E-ball, E-pair-t or E-cmp excluded the true reach in any of the
  7900 draws, under any noise model, either display, and all three tolerances.
- E-pair-pool failed on reach only under t3 noise, in 3 draws at each R.
- The power runs (Section 4) add 0 false labels over 1120 records (280 displays × 4
  engines, each at 3 tolerances × 5 values of k′).

### 3.1 G2-V: validity under each engine's own model

- **E-pair-t and E-pair-pool pass everywhere,** with failure rates of 0.000–0.023.
- **E-ball fails the original clauses in two configurations:**
  - digits gauss: 57/1000, rate 0.057;
  - surface gauss_het: 15/300, upper bound 0.076.
- **Why this is not a bug.** Every one of E-ball's failures is a ball event. No pairwise
  interval or reach interval of E-ball was false in any draw, even at κ = 16, where 80% of
  draws had a failing ball. The ball layer's design rate is exactly
  1 − (1 − δ/n)^n = 0.0488, because the n ball events are independent. The observed rates
  are consistent with it: p = 0.17 and 0.54 against δ = 0.05.
- **The flaw is in the clauses.** They fail a correct engine of this kind about 40% of the
  time per configuration (Amendment 1).
- **Verdict.** Under Amendment 1 (exact one-sided binomial test at level 0.01), every engine
  passes in every configuration. The original wording's verdict is recorded here and in
  `gate2_analysis.txt`.

### 3.2 G2-R: robustness, default engine and scope

The robustness set is digits with R = 8. The criterion is a failure rate ≤ 0.10 on each
member.

| engine | gauss_het | t5 | rep_t5 | aniso4 | verdict |
|---|---|---|---|---|---|
| E-ball | 10/300 | 13/300 | 6/300 | **66/300** | not robust |
| E-pair-t | 2/300 | 0/300 | 0/300 | 6/300 | **robust** |
| E-pair-pool | **299/300** | **63/300** | **230/300** | 13/300 | not robust |
| E-cmp | 0/300 | 0/300 | 0/300 | 0/300 | robust, but has no finite-sample guarantee |

- **Default engine: E-pair-t.** E-cmp is excluded from the default by the plan (Sections 2.2
  and 4), since its guarantee is approximate under every model. The exclusion is also moot:
  E-pair-t has the higher primary metric in all 12 digits configurations. On the synthetic
  panels it is higher at k′ = k and at least as high at k′ = 2k.
- **What the table shows about the engines:**
  - Studentizing makes the pairwise engine insensitive to heavy tails. Its failure rate
    under t5 and rep_t5 noise was no higher than under Gaussian noise, as expected of
    self-normalized statistics. The crude single-jump estimate behind review finding 1 was
    far too pessimistic for pairs.
  - The ball layer's own events break under anisotropy (22% at κ = 4, 80% at κ = 16).
  - E-cmp never failed anywhere, because its variance bounds are very conservative.
- **Scope for the claim:**
  - E-pair-t's guarantee is valid under M_het.
  - It holds empirically under heavy-tailed, replicate-outlier and moderately anisotropic
    noise.
  - It breaks at strong anisotropy: 38/300 (13%) at κ = 16 with R = 8, and 63/300 (21%)
    with R = 16.
  - Gate 3's residual diagnostic therefore has to measure the eigenvalue spread of the
    pooled residual covariance.
  - The immediate remedy for spread beyond about 4 is to inflate τ̂ by √λ_max(Σ̂₀),
    estimated from the pooled residuals. That is valid under a common-shape model, and it
    is the first follow-on if a real dataset needs it.

## 4. Power (plan Section 5.6)

**Primary metric:** the resolved share (faithful + distorted) at k′ = 40, ε = 0.2, on digits
with t-SNE (3 seeds × 5 draws).

| noise | R | E-ball | E-pair-t | E-pair-pool | E-cmp | display faithful (truth) |
|---|---|---|---|---|---|---|
| gauss, γ = 0.15 | 8 | 0.092 | **0.326** | 0.364 | 0.242 | 0.710 |
| gauss, γ = 0.15 | 16 | 0.226 | **0.504** | 0.521 | 0.456 | 0.708 |
| gauss, γ = 0.15 | 32 | 0.391 | **0.642** | 0.648 | 0.623 | 0.707 |
| gauss, γ = 0.30 | 8 | 0.003 | 0.049 | 0.075 | 0.019 | 0.705 |
| gauss, γ = 0.30 | 16 | 0.023 | 0.186 | 0.203 | 0.131 | 0.710 |
| gauss, γ = 0.30 | 32 | 0.106 | 0.355 | 0.366 | 0.316 | 0.711 |
| gauss_het, γ = 0.15 | 8 / 16 / 32 | 0.078 / 0.210 / 0.369 | 0.284 / 0.458 / 0.602 | 0.294 / 0.455 / 0.599 | 0.197 / 0.403 / 0.575 | 0.705 |

**G2-P passes:** 0.326 ≥ 0.20 at R = 8, and 0.504 ≥ 0.30 at R = 16.

**Where the resolution goes, at γ = 0.15 and ε = 0.2.** E-pair-t's faithful / distorted
shares, with the display's true faithful share:

| k′ | R = 8 | R = 16 | R = 32 | truth |
|---|---|---|---|---|
| 20 (the predecessor's target) | 0.015 / 0.263 | 0.045 / 0.413 | 0.077 / 0.531 | 0.24 |
| 40 | 0.259 / 0.067 | 0.398 / 0.107 | 0.494 / 0.148 | 0.71 |
| 80 | 0.658 / 0.021 | 0.763 / 0.032 | 0.816 / 0.049 | 0.90 |

- **Distorted statements are available early.** Where faithful statements are rare, at
  k′ = k, E-pair-t already certifies 26–53% of points as distorted.
- **Interval width.** The median reach interval spans a factor of 3.4 (R = 8), 2.4 (R = 16)
  and 1.8 (R = 32). For E-ball the factors are 6.4, 4.1 and 2.9.
- **Cost.** One E-pair-t call takes 0.2 s at n = 1000. The pooled variant gains at most 0.04
  in resolved share, and only under the homoscedastic noise it needs.

## 5. What this does to the predecessor's conclusions

**Judgment call.** This section records a result that narrows the ground earlier
conclusions stand on, so it is a section rather than a caveat.

**The predecessor's conclusions.** It concluded that roughly two thirds of coverage is lost
to the width of the confidence sets. Its best valid sharpening of those sets was worth
+0.066 in ceiling coverage (its gate report, Section 13), which withdrew the confidence-set
paper as an option. Every confidence set it studied was a ball around each object.

**The same-panel measurement.** `results/gate2/synthetic_ceiling.txt` uses the power runs'
synthetic panels (4 families × 2 noise levels, PCA-2, R = 8), at the predecessor's own
target k′ = k = 20, ε = 0.2:

| | ball interval | ball witness ceiling (any ball-based certificate) | E-pair-t |
|---|---|---|---|
| pooled over the 8 panels | 0.139 | 0.303 | **0.412** |
| range over panels | 0.003–0.492 | 0.026–0.737 | 0.043–0.797 |

- **Every panel.** E-pair-t exceeds the ball ceiling on all 8 panels: +0.109 pooled, and
  +0.273 over the ball interval.
- **Digits.** At k′ = 40 and R = 8, E-pair-t certifies 0.259 faithful against a ball ceiling
  of 0.200 on the same displays.
- **Why.** A pairwise statement only has to absorb the noise projected on the one direction
  joining the pair. Over n(n−1)/2 events that costs about 5.7 standard deviations of the
  pair's projected noise, a half-width of 8.0 σ/√R. A ball must hold in every direction at
  once, so its radius is √(Dq) · σ/√R = 7.9 σ/√R for each object, and the ball interval's
  half-width ρ_i + ρ_j is 15.7 σ/√R.

**What it changes.**
- The predecessor's finding that the confidence-set width binds holds for *ball-shaped*
  confidence sets only. The radius-sensitivity analysis ruled out sharpening balls, not
  the confidence-set paper in general.
- Its digits result at the exact target stands. E-pair-t reaches 0.015, 0.045 and 0.077
  faithful at k′ = k for R = 8, 16 and 32, still under the plan's 0.10 floor. Rank
  tolerance, not the engine, is what makes digits informative.
- The predecessor's records (its gate report and the unsubmitted report) are not edited
  here. That repository is closed.

## 6. Need (G2-N, recorded; decisive in Gate 3)

**Naive plug-in mislabel rate, digits with t-SNE, `gauss` noise, ε = 0.2.** In parentheses:
the share of naive-faithful points that are actually distorted.

| γ, R | k′ = 20 | k′ = 40 | k′ = 80 |
|---|---|---|---|
| 0.15, 8 | 0.035 (0.072) | 0.020 (0.014) | 0.006 (0.004) |
| 0.15, 16 | 0.027 (0.049) | 0.019 (0.013) | 0.005 (0.004) |
| 0.30, 8 | 0.057 (0.128) | 0.040 (0.029) | 0.015 (0.009) |
| 0.30, 16 | 0.044 (0.089) | 0.031 (0.022) | 0.010 (0.006) |

- **Other cases.** The bias-corrected naive map gives the same numbers to within 0.001.
  Heteroscedastic noise adds at most 0.005. On the synthetic families, whose PCA-2 displays
  are nearly faithful at intrinsic dimension 1–2, the naive map mislabels at most 0.5% of
  points.
- **What the certificate adds.** `results/gate2/need_crosstab.txt` covers 15 t-SNE displays
  per configuration at k′ = 40:
  - **It never corrected a naive error.** At γ = 0.15, R = 8 the naive map is wrong on 293
    of 15000 point-labels, and E-pair-t abstains on all 293. The same holds in all four
    configurations, with 278 to 602 naive errors in each.
  - **It confirmed a third to a half of the correct labels:** 33% of them at γ = 0.15,
    R = 8, and 51% at R = 16. At γ = 0.30 it confirms 5% (R = 8) and 19% (R = 16).
  - **It never made a wrong statement.**

**What this means.**
- At these noise levels the naive errors lie inside the certificate's noise margin by
  construction. The margin is about 5.7 noise standard deviations, while the naive map's
  errors come from distance perturbations of about one.
- So certification cannot find the naive map's mistakes. It can only mark the part of the
  map where no mistake is possible, and leave the rest (including every mistake) unresolved.
- On this panel the unresolved part is about 97% correct anyway.

**Consequence for the plan.**
- G3-need requires a plug-in mislabel rate ≥ 0.05 at k′ = 2k on at least one real dataset.
  Digits under injected Gaussian noise would fail it (0.020–0.040). If real replicate data
  looks like this panel, the project stops at Gate 3.
- Where it could pass: real trial-to-trial noise is often large relative to stimulus
  spacing. Split-half noise ceilings in representational-similarity studies are typically
  well below 1, and the Gate 1 audit found downsampled single-cell kNN graphs sharing only
  about 2.9 of 50 neighbors with deep-sequencing ones. In those regimes naive maps err more,
  and the certificate needs more replicates to resolve anything. Gate 3 measures where real
  data sits.

## 7. Planning rule (recorded, not gated)

**The quantity.** The inflation on the distance-ratio scale,
Δ_i = (r_hi/k)^{1/d} − (reach_i/k)^{1/d}, divided by the engine's half-width in units of the
point's spacing, h/d0_i. Its median and quartiles:

- **Digits, t-SNE.** c = 1.27–1.59 for E-pair-t across R ∈ {8, 16, 32} and γ ∈ {0.15, 0.30},
  with interquartile ranges inside [0.85, 2.21]. E-ball gives the same c, 1.26–1.47, in
  units of its own half-width.
- **Synthetic, PCA-2, E-pair-t.** c = 1.39–1.71 at γ = 0.30, but 0.46–1.07 at γ = 0.15. At
  low noise many points' upper bounds equal their true reach exactly, so the continuum model
  does not apply.

**Verdict.** Against the plan's usefulness criterion (c within a factor of 2 across
geometries and R), the rule holds on digits and at γ = 0.30. It fails at low noise on the
synthetic families, where c spans a factor of 3.7.

**The planning rule as supported.** Set r_hi/k ≈ [(reach/k)^{1/d} + 1.4 · h/d0]^d. Here h is
closed-form in R: h = t · √(2σ²/R), with t between 5.4 and 5.7 at n = 1000 (it depends on
R through ν). For points far from the noise
floor the rule is a conservative, order-of-magnitude guide.

## 8. Deviations and notes

- **Amendment 1** (the G2-V criterion) was written after the validity runs and before the
  power runs. It does not affect the choice of engine.
- **The pooled engine passed G2-V but is not a candidate.** It breaks under heteroscedastic
  and heavy-tailed noise, so the plan's robustness rule excludes it.
- **A refinement from the Gate 1 audit is deferred.** Inverting the noncentral χ² would
  tighten E-pair-t's lower bound by about 1.3τ out of about 6.8τ at digits R = 8. It needs a
  studentized worst-case argument not yet written, and it would not change any gate
  decision.
- **Stepwise baseline, deferred.** The Mogstad et al. stepwise procedure is the engine
  baseline Gate 1 asked for. It is left for after Gate 2.
- **Wall time.** Validity took 544 s on 30 workers. Power took 81 s on 32 workers.

## 9. Next

1. **Gate 3 needs real replicate data, and every candidate needs a download.** Candidates
   from the Gate 1 audit (sizes are the audit's estimates, unverified):
   - Allen Brain Observatory natural scenes: 118 images × 50 trials; about 0.1–1 GB per
     two-photon session.
   - THINGS-EEG2 test set: 200 images × 80 repetitions per subject; the preprocessed split
     is small.
   - TVSD: 100 images × 30 repetitions, V1/V4/IT multi-unit activity.
2. **Because need decides the project, the cheapest decisive step is a need pilot on one of
   them.**
   - The pilot certifies a reference labeling from the held-out replicates (R = 40 or more)
     with E-pair-t itself.
   - It then measures the naive R = 8 map's mislabel rate on the points that reference
     resolves.
   - Before any of it: the residual diagnostics of G3-diag.
