# Certified neighborhood reach: research plan

Version 1, 2026-09-28. Written before any Gate 2 measurement. The gate criteria in Section 6
are fixed by this version. A change made after data is seen is an amendment: it is dated,
gives its reason, and leaves this text in place.

## 0. Provenance

This project succeeds `setwise-certificates-neighborhood-fidelity`, which was dropped on
2026-09-23 (its gate report, Section 15). Its evidence base is frozen at commit `9bc9f03`
(tarball SHA-256 `8e7336bc597c340f493720576c0df019f6aea3e7892d9eaf5625f9d3eab5ddb9`). On
2026-09-28 its negative-results report was withdrawn from submission and five successor
items were derived from it. The combined item, reviewed and revised in Section 2, is the
subject of this plan.

A scratch pilot on 2026-09-28 supplied the numbers the proposal rested on. Its scripts and
log are archived unchanged in `pilot/2026-09-28/`. They import the predecessor repository
by absolute path and are not rerun here; Gate 2 re-derives their numbers with this
repository's code (Section 5.7).

## 1. Question and claim

**Question.** A two-dimensional display (t-SNE, UMAP, PCA) of noisy measurements is read
neighbor by neighbor. Which displayed neighborhoods reflect the neighborhoods of the
noise-free objects? Which certainly do not? How far off is each one?

**Claim to be earned by the gates.** Take data with R replicate measurements per object and
any display fitted to the replicate means. The method returns, for every point i at once,
an interval `[r_lo(i), r_hi(i)]`. With probability at least 1 − δ over the noise, every
interval contains that point's *neighborhood reach*: the ⌈(1−ε)k⌉-th smallest true rank
among its k displayed neighbors. Coverage is exact in finite samples under a stated noise
model. It holds for any display, including one fitted to the same data, and for all
tolerances (ε, k′) simultaneously. At any tolerance the intervals give a three-way map:
*certified faithful* (`r_hi ≤ k′`), *certified distorted* (`r_lo > k′`), or *unresolved*.

**Not claimed:**
- an inequality of independent interest;
- an engine beyond elementary constructions (Section 4);
- dataset-level quality curves (summing per-point simultaneous bounds gives needlessly
  conservative dataset bands).

The contribution is the statement, its exactness, and the measured characterization of when
it is informative (Sections 5.6 and 5.8).

## 2. Review of the 2026-09-28 proposal

### 2.1 What the proposal rested on

The setting is digits with a t-SNE display fitted to replicate means, noise γ = 0.15,
tolerance ε = 0.20, and three t-SNE seeds on one noise draw. The table gives the share of
1000 points.

| k′ | R | display faithful | ball interval, faithful | ball ceiling | comparison, faithful | comparison, distorted |
|---|---|---|---|---|---|---|
| 20 | 8 | 0.237 | 0.002 | 0.009 | 0.004 | 0.135 |
| 40 | 8 | 0.702 | 0.075 | 0.194 | 0.188 | 0.053 |
| 40 | 16 | 0.705 | 0.183 | 0.335 | 0.351 | 0.099 |
| 40 | 32 | 0.708 | 0.313 | 0.452 | 0.481 | 0.143 |

"Comparison" is the pilot's Bonferroni test on bias-corrected squared-distance differences.
The same test with plug-in variances and a Bonferroni count over the used comparisons only,
which is not valid for a data-dependent display, reached 0.375 faithful at k′ = 40, R = 8.
On the synthetic panels (PCA-2, R = 8) the ball ceiling rose from 0.32 at k′ = 20 to 0.91 at
k′ = 40. No false statement occurred in any cell.

### 2.2 Findings and revisions

| # | finding | severity | revision |
|---|---|---|---|
| 1 | The proposed engine is valid only approximately. It needs n·C(n−1,2) ≈ 5×10⁸ simultaneous two-sided tests at n = 1000, a per-test level near 5×10⁻¹¹. At that level the pure-noise bilinear term has a tail up to about 10⁵ times the Gaussian one (Chernoff bound 9×10⁻⁶ against 5×10⁻¹¹). A crude single-jump estimate puts a t5 tail about 10⁵ times too heavy for comparisons, about 500 times for pairs, and about 2 times for balls. | major | An engine ladder by family size (Section 4): balls (n events), pairwise distance intervals (n² events), comparisons (n³ events). The default is the most powerful engine that is exact under the stated model and passes the robustness checks. The approximate comparison engine stays only as a reference. |
| 2 | "Balls are the fallback valid by construction" holds only under isotropic Gaussian noise. Every engine inherits a noise model, and real replicate noise is heteroscedastic, anisotropic and heavy-tailed. | major | Gate 2 measures family-wise error under seven noise models. The real-data gates must check the model on residuals and state the scope. |
| 3 | Need is not established. Nothing yet measures how wrong a naive map (reach computed from the noisy means) is against the truth. If naive maps are nearly right, the inference layer solves a non-problem. | major | Naive and bias-corrected plug-in maps are baselines everywhere. Gate 3 has a need criterion. |
| 4 | The target was defined through overlap counts at each k′. It is simpler as an order statistic of true ranks, which makes the bounds closed form and shows they are uniform in (ε, k′). Choosing a tolerance after seeing the map therefore cannot break validity, only reporting. | moderate | Definition in Section 3. Box engines get exact closed-form reach bounds (members at their upper bounds, outsiders at their lower bounds), tighter than counting individually certified members. The full interval is reported, not a chosen k′. |
| 5 | The required rank tolerance grows exponentially with intrinsic dimension, while the noise adds roughly a constant on the distance-ratio scale (k′/k)^{1/d}. | moderate | The planning rule is stated and calibrated on that scale (Section 5.8). |
| 6 | The earlier prior-art scan missed two literatures. Visual-analytics distortion maps (Aupetit 2007; CheckViz 2011; Martins et al. 2014; the Nonato and Aupetit survey 2019; Jeon et al.) have the same output form, without noise. Uncertain-database possible/certain kNN queries are the interval engine under another name. | moderate | Gate 1 audits the combination and both literatures before any claim of novelty. |
| 7 | Data thinning does not belong in the first paper. The truth it recovers (rates in a linear representation) is not what log-normalized displays are built on, and fold independence needs a known count model. | moderate | Thinning moves to a follow-on. The real-noise gate uses a dataset with many genuine replicates, so held-out replicates can form a near-truth reference (Gate 4). |
| 8 | The proposed usefulness floor (resolved share ≥ 0.30 at k′ = 2k, R ≤ 16) is already passed by the development dataset (0.45), so on digits it tests nothing. | minor | It is kept as a sanity floor on development data (Gate 2) and as the confirmatory floor on new data (Gate 3). |
| 9 | Reach is conditional on the display the user has, including its random seed. | minor | Stated in the claim. Power is reported over three display seeds. |

### 2.3 The revised item

**Certified neighborhood reach.** For 2-D displays of replicate-measured data, the method
gives simultaneous finite-sample intervals on how far each displayed neighborhood reaches
into the true neighbor ranking. From them it draws a faithful / distorted / unresolved map
valid at every tolerance. The engine is the most powerful exact member of a
ball / pairwise / comparison ladder that survives the robustness checks. Two measurements
back it: how often naive maps mislabel, and a planning rule for the replicates a tolerance
needs. Adaptive re-measurement of unresolved points and thinning for count data are
follow-ons.

## 3. Setting and target

- **Objects.** i = 1..n have true means μ_i ∈ ℝ^D. Replicates are x_ir = μ_i + e_ir for
  r = 1..R, independent across r and i.
- **Model M_het.** e_ir ~ N(0, σ_i² I_D), with unknown object-specific σ_i.
- **Model M_hom.** M_het with σ_i ≡ σ.
- **Statistics.**
  - c_i is the mean of x_i1..x_iR;
  - s_i² = Σ_r ‖x_ir − c_i‖² / (D(R−1)), with ν = D(R−1) degrees of freedom;
  - the pooled variance s̄² is the mean of the s_i², with nν degrees of freedom.
- **Display.** Any map from the data to ℝ². E_i is the set of i's k displayed nearest
  neighbors, and it may depend on the same noise. Every engine's statements hold
  simultaneously over all pairs, or all triples, so no selection correction is needed.
- **True rank.**
  - rank_i(j) is j's 1-based position in the order of ‖μ_i − μ_l‖ over l ≠ i, with ties
    broken by object ID as in the predecessor.
  - With m = ⌈(1−ε)k⌉, **reach** is reach_i = the m-th smallest of {rank_i(j) : j ∈ E_i}
    = m + #{l ∉ E_i ∪ {i} ranked before the m-th displayed neighbor}.
- **Faithful at (ε, k′).** reach_i ≤ k′, which is equivalent to |E_i ∩ Top_{k′}(μ_i)| ≥ m.
  At k′ = k this is the predecessor's target, realized-ok.
- **Distance-ratio scale.** (reach_i/k)^{1/d}, where d is a local intrinsic dimension. It is
  used for planning, never for certification.
- **Naive maps.**
  - *Plug-in*: the same reach computed with c in place of μ.
  - *Bias-corrected*: ranks by ‖c_i − c_l‖² − (D/R) s_l².
  - Both label every point, with no abstention. Their *mislabel rate* at (ε, k′) is the
    share of points whose naive label disagrees with the truth.

## 4. Engines

Every engine yields, on an event of probability at least 1 − δ, one of two things:
- intervals [L_ij, U_ij] containing ‖μ_i − μ_j‖ for all pairs (box engines);
- certified comparisons for all triples (comparison engines).

Reach bounds follow, with U_(m)(i) and L_(m)(i) the m-th smallest U_ij and L_ij over j ∈ E_i.

- **Box engines:**
  - `r_hi(i) = m + #{l ∉ E_i ∪ {i} : L_il ≤ U_(m)(i)}`;
  - `r_lo(i) = m + #{l ∉ E_i ∪ {i} : U_il < L_(m)(i)}`.
  - These are the exact extremes of reach over the box: members at their upper bounds,
    outsiders at their lower bounds, and the reverse. They are valid with ID tie-breaking.
- **Comparison engines:** r_hi(i) and r_lo(i) are the m-th smallest over j ∈ E_i of
  1 + #{l : l not certified farther than j} and 1 + #{l : l certified closer than j}.

Every (ε, k′) label derives from the same event.

**E-ball** is the predecessor's F-pivot balls:
- ρ_i = s_i √(Dq/R), with q = F⁻¹_{D,ν}(1 − δ/n) verified in arbitrary precision;
- L_ij = max(0, ‖c_i − c_j‖ − ρ_i − ρ_j) and U_ij = ‖c_i − c_j‖ + ρ_i + ρ_j.

It is exact under M_het.

**E-pair-t** (new, pairwise, studentized). For every pair:

    U_ij = ‖c_i − c_j‖ + t · τ̂_ij
    L_ij = max(0, sqrt(max(0, ‖c_i − c_j‖² − q_W · τ̄_ij²)) − t · τ̂_ij)

The ingredients:
- τ̂_ij² = (s_i² + s_j²)/R, and τ̄_ij² = (b_i + b_j)/R, where b_i = s_i² ν / χ²_ν⁻¹(δ_σ/n)
  bounds σ_i² from above;
- t = T_ν⁻¹(1 − δ_Z/(n(n−1))) and q_W = χ²_{D−1}⁻¹(1 − 2δ_W/(n(n−1)));
- the error budget splits as δ = δ_Z + δ_W + δ_σ = 0.040 + 0.005 + 0.005;
- every quantile is verified in arbitrary precision, as the predecessor's F quantile was.

Validity argument, under M_het:
1. Let τ² = (σ_i² + σ_j²)/R. Then c_i − c_j = (μ_i − μ_j) + τ(Z u + √W v), where u is the
   unit vector along μ_i − μ_j and v is a unit vector orthogonal to it.
2. Z ~ N(0,1) and W ~ χ²_{D−1} are independent of each other and of (s_i², s_j²).
3. **Upper bound.** ‖c_i − c_j‖ ≥ d_ij + τZ, so U_ij holds whenever −τZ ≤ t τ̂_ij.
4. **Lower bound.** ‖c_i − c_j‖² ≤ (d_ij + τ|Z|)² + τ²W, so L_ij holds on the event
   {τ|Z| ≤ t τ̂_ij, W ≤ q_W, τ ≤ τ̄_ij}.
5. **Why the Student quantile is conservative.** τ̂²/τ² = wA + (1−w)B, where A and B are
   i.i.d. χ²_ν/ν and w = σ_i²/(σ_i² + σ_j²). The map v ↦ Φ̄(t√v) is convex on (0, ∞): its
   second derivative is t⁴ φ(u)(u² + 1)/(4u³) with u = t√v. By Jensen,
   P(Z > t √(wA + (1−w)B)) ≤ P(T_ν > t) for every w. The Student quantile is therefore
   conservative whatever σ_i and σ_j are. (A 4×10⁶-draw Monte Carlo at ν = 20 confirms
   the worst case is w ∈ {0, 1}.)
6. A union over pairs, both sides and the n variance bounds gives the family-wise statement.

Studentizing costs a factor 1.056 at R = 8, against 1.340 for bounding each σ_i separately.
The upper half-width is 8.0 σ/√R, against 15.7 σ/√R for the ball interval.

**E-pair-pool** is E-pair-t with every s_i² replaced by s̄² and ν replaced by nν. It is exact
under M_hom only.

**E-cmp** (reference only, approximate) is the pilot's conservative comparison engine:
- bias-corrected squared-distance differences;
- a Bonferroni count over all n·C(n−1,2) comparisons, with δ split three ways;
- variances from σ upper bounds and ball-inflated distances, and a Gaussian quantile.

It is exact under no model (finding 1). It measures what exactness costs.

**Reference quantity, not a certificate.** The predecessor's ball witness ceiling at k′ = 2k
(t-SNE displays only) keeps continuity with the predecessor's numbers.

## 5. Gate 2 design (development data only)

### 5.1 Geometries

- **`digits`.** Exactly the predecessor's real-data panel: scikit-learn's bundled digits,
  split seed 20260923, zero-variance drop, standardization and PCA20 fitted on the
  797-object development split, n = 1000 evaluated objects, D = 20, and d0 = the median
  true 20th-neighbor distance.
- **`surface`.** The predecessor's `curved_surface` generator (n = 1000, D = 20). Its
  layout seed is 3_310_000 and its d0 comes from a development draw with seed 3_300_000.
- **`digits200`.** The first 200 objects of the digits panel. It is used only to check
  every elementary comparison of E-cmp.

### 5.2 Noise models

The noise scale is σ = γ d0/√D, as in the predecessor, with γ = 0.15 unless stated.

| name | noise |
|---|---|
| `gauss` | M_hom |
| `gauss_het` | M_het with σ_i = σ w_i, log w_i ~ U(−ln 2, ln 2), redrawn each draw |
| `t5`, `t3` | independent Student-t coordinates scaled to unit variance |
| `rep_t5` | whole-replicate scale mixture: e_ir = σ √g_ir z_ir, with g_ir ~ 3/χ²_5 and z_ir ~ N(0, I_D) (one bad acquisition inflates every coordinate of one replicate) |
| `aniso4`, `aniso16` | N(0, σ² Σ₀), with Σ₀ = Q diag(λ) Qᵀ, λ geometric with max/min 4 or 16 and mean 1, and Q a fixed rotation (seed 3_400_000 + κ) |

A noise-dominated stress runs `gauss` and `t5` at γ = 0.5 on digits with R = 8.

### 5.3 What the validity runs check

- **Box engines:** every pairwise interval (n(n−1)/2 statements) and every point's reach
  interval.
- **E-ball:** additionally, every ball.
- **E-cmp:** at n = 1000, every derived rank and reach statement; at n = 200 (digits200),
  every elementary comparison.
- **Displays:** reach statements use two displays per draw, PCA-2 of c and the noisy-kNN
  display (the k nearest in c). The second follows the noise exactly. Validity does not
  depend on the display.

A draw is a *family-wise failure* for an engine if any of its statements is false.

### 5.4 Validity configurations and draws

| geometry | noise | R | draws |
|---|---|---|---|
| digits, γ = 0.15 | gauss | 8 | **1000** (primary) |
| digits, γ = 0.15 | gauss | 16 | 300 |
| digits, γ = 0.15 | gauss_het, t5, t3, rep_t5, aniso4, aniso16 | 8, 16 | 300 each |
| surface, γ = 0.15 | gauss, gauss_het, t5, aniso4 | 8 | 300 each |
| digits, γ = 0.5 | gauss, t5 | 8 | 300 each |
| digits200, γ = 0.15 | gauss, t5, aniso4, rep_t5 (E-cmp elementary checks) | 8 | 300 each |

### 5.5 Power configurations

- **Digits.** Noise `gauss` and `gauss_het`, γ ∈ {0.15, 0.30}, R ∈ {8, 16, 32}. Each has
  5 draws, with a t-SNE display (perplexity 50, PCA initialization, seeds 17, 29, 31) and a
  PCA-2 display, both fitted to c.
- **Synthetic.** The predecessor's four primary families with `gauss`, γ ∈ {0.15, 0.30},
  R = 8, a PCA-2 display, and 5 draws.
- **Metrics:**
  - faithful, distorted and unresolved shares at ε ∈ {0.1, 0.2, 0.3} and
    k′ ∈ {20, 30, 40, 60, 80};
  - the median of log2(r_hi/r_lo);
  - naive and bias-corrected mislabel rates;
  - wall time per engine;
  - the ball witness ceiling at k′ = 40 for t-SNE displays at γ = 0.15.

### 5.6 Primary metrics

- **Primary power metric.** The *resolved share* (faithful + distorted) at k′ = 2k = 40,
  ε = 0.2, on digits with t-SNE and `gauss`, averaged over draws and display seeds.
- **Need metric.** The plug-in mislabel rate at the same tolerance.

### 5.7 Structural tests, before any Gate 2 run

1. Reach computed from its definition equals the closed form, with ties broken by ID.
2. Box reach bounds contain the truth whenever the box does. They equal it when the box is
   degenerate. The adversarial configuration attains `r_hi` when there are no ties.
3. Zero-noise control: identical replicates give exact reach bounds in every engine.
4. Monotonicity: intervals widen as δ falls.
5. Reproduction of the pilot at (digits, γ = 0.15, acquisition seed 710_000, t-SNE seed 42,
   R = 8, k′ ∈ {20, 40}, ε = 0.2):
   - E-cmp's faithful and distorted shares equal the pilot's comparison bounds exactly;
   - E-ball's faithful share equals the pilot's interval bound exactly;
   - E-ball's distorted share is at least the pilot's, because the box bound is tighter
     than the pilot's per-member count.

   The one number this run reveals for a new engine is `pair_t` at the pilot's single
   configuration. It is logged in the Gate 2 report as seen before Gate 2 proper.

### 5.8 Planning-rule calibration (recorded, not gated)

- On the distance-ratio scale, the inflation Δ_i = (r_hi(i)/k)^{1/d} − (reach_i/k)^{1/d},
  with d the panel's median Levina–Bickel dimension, is modeled as c · h/d0_i.
- h is the default engine's typical upper half-width, and d0_i is the true k-th-neighbor
  distance of i.
- The distribution of c_i over points is recorded per geometry and R. The rule is useful if
  c stays within a factor of 2 across geometries and R.

## 6. Gate criteria (pre-registered)

**Gate 1 (prior-art audit of the combination).**
- It fails if one prior work states per-point simultaneous confidence statements about a
  display's neighborhood fidelity against the noise-free truth (C1 together with C2 or C3
  in `docs/gate1_prior_art.md`).
- If only components are anticipated, it passes. The novelty statement is then narrowed to
  what was not found.

**G2-V (exactness; a failure is a bug, not a result).**
- Each engine is checked under every model it claims exactness for:
  - E-ball and E-pair-t: `gauss` and `gauss_het`;
  - E-pair-pool: `gauss`.
- The criterion is: estimated family-wise failure rate ≤ δ = 0.05, and one-sided 95%
  Clopper–Pearson upper bound ≤ 0.075.
- A failure stops Gate 2 until the cause is found and fixed. The fix is logged.

**G2-R (robustness; sets the default engine and the scope).**
- The robustness set is `gauss_het`, `t5`, `rep_t5` and `aniso4` on digits with R = 8.
- An engine is robust if its failure rate is ≤ 0.10 on each member of the set.
- The default engine is the robust engine with the highest primary power metric at R = 16.
- If no engine is robust, the default is the most powerful engine exact under M_het. The
  claim is then restricted to Gaussian noise with a mandatory residual diagnostic, and a
  robust engine becomes a precondition for Gate 3.

**G2-P (power; kill criterion).**
- The default engine's primary power metric must be at least 0.20 at R = 8 and at least
  0.30 at R = 16.
- If it is not, the project stops: new data cannot rescue a method that misses its floor
  on the development panel.

**G2-N (need; recorded).** The plug-in mislabel rate on digits with t-SNE, at γ = 0.15 and
γ = 0.30, R ∈ {8, 16}.

**Gate 2 passes** if G2-V and G2-P pass.

**Gate 3 (real data, genuine replicates; registered now, specifics amended before any data
is loaded).**
- Scope: at least two datasets with genuine replicates (R_total ≥ 8), with t-SNE and PCA-2
  displays of c. The downloads need the owner's permission.
- **G3-need:** plug-in mislabel rate ≥ 0.05 at k′ = 2k, ε = 0.2 on at least one dataset at
  its native R.
- **G3-use:** resolved share ≥ 0.30 at k′ = 2k, ε = 0.2 with R ≤ 16 on at least two
  datasets.
- **G3-diag:**
  - residual kurtosis, anisotropy (the eigenvalue spread of the pooled residual covariance)
    and the dispersion of s_i² are recorded;
  - claims are restricted to what the diagnostics support.

**Gate 4 (validity against real noise).**
- The dataset needs R_total ≥ 40. Certificates from R ∈ {8, 16} subsamples are checked
  against a reference built from the held-out replicates, counting only the statements the
  reference itself resolves.
- The criterion is a failure rate ≤ 0.10 over subsamples.

## 7. Seeds and data hygiene

- **Seed bases:**
  - 3_100_000 for validity;
  - 3_200_000 for power;
  - 3_300_000 and 3_310_000 for synthetic development and layout draws;
  - 3_400_000 for the anisotropy rotations;
  - t-SNE seeds 17, 29 and 31.

  A configuration's draws use base + 10_000 · (configuration index) + draw.
- **Disjoint from the predecessor.** The predecessor's bases are 100_000, 123_000, 710_000
  and 900_000, plus the digits split seed 20260923 and small fixed seeds (t-SNE seeds 42,
  123 and 456). This repository reuses:
  - the split seed, to reproduce the digits panel;
  - acquisition seed 710_000, once, for the reproduction test;
  - the predecessor's ambient-map seed, as part of its frozen geometry.
- **Reserved.** Seeds from 4_000_000 up are for confirmatory work and are not drawn in
  Gate 2.
- **Data.** Only data bundled with installed libraries is used before Gate 3. Nothing is
  downloaded without the owner's permission.

## 8. Stop conditions

- Gate 1 finds the combination anticipated.
- The default engine misses G2-P.
- On real data, G3-use fails, or G3-need fails on every dataset.

## 9. Follow-ons, not in the first paper

- Adaptive, anytime-valid replicate allocation to unresolved points.
- Count data through data thinning, with its own model checks.
- An exact comparison engine (n³ family). It is worth building only if the comparison
  family survives G2-R and the n³ family's shared-query advantage exceeds its multiplicity
  cost.

## Amendments

### Amendment 1 (2026-09-28): the G2-V criterion

**When.** Made after the Gate 2 validity runs had been seen, and before any power run. It
does not affect the choice of default engine (Section 6, G2-R and G2-P).

**What went wrong.** G2-V as written has two clauses: an estimated failure rate ≤ δ, and a
one-sided 95% Clopper–Pearson upper bound ≤ 0.075. Both are miscalibrated for an engine
whose family-wise failure rate equals its design value.

- **The ball layer is such an engine by construction.** Its n ball events are independent,
  and each is exact at δ/n, so their union has probability exactly
  1 − (1 − δ/n)^n = 0.0488 at n = 1000 (and at n = 200).
- **How often a correct engine fails.** Its estimated rate exceeds 0.05 with probability
  0.39 in each configuration, at both 300 and 1000 draws. At 300 draws its upper bound
  exceeds 0.075 about half the time.
- **Across configurations.** E-ball is checked in eight configurations, so the clauses
  were close to certain to fail a correctly working ball engine.
- **Where the slack is.** The Bonferroni engines built on dependent pairs (E-pair-t,
  E-pair-pool) have slack. The ball layer has none.

**The replacement.** An engine fails G2-V under its own model if an exact one-sided binomial
test rejects "failure rate ≤ δ" at level 0.01 in any configuration. The Gate 2 report
gives the verdict under the original clauses next to the amended one. The amended criterion
also applies to Gate 4.

### Amendment 2 (2026-09-29): Gate 3 need pilot on THINGS-EEG2

**When and why.** Written before any THINGS-EEG2 data value was loaded. The owner authorized
this one dataset for a need pilot, the cheapest decisive step named in the Gate 2 report.
Section 6 requires Gate 3 specifics to be fixed before data is loaded, and this amendment
fixes them. One dataset cannot pass Gate 3, which needs two, but it can stop the project.

**The data, as documented by the authors' code (not yet inspected):**
- Source: the preprocessed test set of Gifford et al. (2022), OSF node `anp5v`. Only the
  test file inside each subject's zip was fetched; `fetch_things_eeg2_test.py` checks each
  file's CRC-32 against the zip.
- Shape: 10 subjects, each with 200 images × 80 repetitions × 17 occipital/parietal
  channels × 100 samples (−0.2 to 0.79 s at 100 Hz).
- Preprocessing: baseline-corrected; channels whitened within each session (MVNN, one
  covariance averaged over time points).
- Repetition order: 4 sessions of 20 repetitions, merged, then shuffled across sessions with
  seed 20200220. That shuffle reconstructs each repetition's session.

**Design (per subject, per split; 10 random splits with seeds 3_500_000 + 1000·s + split):**
1. **Blocks.** Repetitions are divided within each session into F (2 per session, 8
   total), U (8 per session, 32) and H (10 per session, 40). Every block therefore draws
   equally on the four sessions.
2. **Representation, fitted on F only.** Take the post-stimulus samples (t ≥ 0: 80 samples ×
   17 channels = 1360 features).
   - PCA of F's image means to 100 dimensions.
   - F's residual covariance in those dimensions (1400 degrees of freedom, Ledoit–Wolf
     shrinkage), then whitening.
   - PCA of the whitened means to D = 20.

   The truth is each image's expected response, averaged over the four sessions and
   projected into that space.
3. **User.** Session-stratified prefixes of U with R ∈ {8, 16, 32}.
   - Displays: t-SNE (perplexity 30, PCA initialization, seeds 17/29/31) and PCA-2, both of
     the user's means.
   - Tolerances: k = 10 (n = 200) and ε = 0.2. The primary k′ is 2k = 20; secondary
     k′ ∈ {10, 40}.
4. **Primary need estimator: mislabel rate extrapolated in the reference.**
   - Each user labeling is compared with plug-in labels from disjoint, session-stratified
     sub-blocks of H of size R_h ∈ {8, 20, 40} (5, 2 and 1 sub-blocks).
   - The disagreement is averaged over splits and seeds, per subject.
   - q̂(R) is the intercept of a least-squares line in R_h^(−1/2).
   - As R_h → ∞ the reference labels converge to the truth. To first order their own error
     rate scales with their noise standard deviation.
5. **Secondary need estimators:**
   - (a) *Certified bounds.* E-pair-t is run on H (R = 40). Naive errors counted on the
     points H resolves give a lower bound. Adding H's unresolved share gives an upper bound.
   - (b) *Semi-synthetic.* The truth is split 0's H means. User data is that truth plus U
     residuals resampled within image, for R ∈ {8, 16, 32}, 30 draws per subject. The
     naive mislabel rate is then exact.
6. **Use.** E-pair-t's resolved share at k′ = 20 on the real user blocks, R ∈ {8, 16, 32},
   and on U ∪ H (R = 72).
7. **Validity under real noise.** In the semi-synthetic setting, E-pair-t's family-wise
   failure rate over 300 draws per subject at R = 8: all pairs, plus reach under the PCA-2
   and noisy-kNN displays.
8. **Diagnostics (G3-diag).** From split 0, on the U ∪ H residuals in the 20-dimensional
   space:
   - excess kurtosis;
   - eigenvalue spread of the pooled residual covariance (max/min and max/mean);
   - dispersion of s_i² relative to its value under homoscedasticity;
   - a session × image interaction statistic, using the reconstructed sessions;
   - mean residual variance, as a check that F's whitening carries over.

**Decision (medians over the 10 subjects; fixed in advance):**
- **Need at native R (80):**
  - *Fails* if q̂(32) < 0.05, since q falls with R.
  - *Passes* if q̂(80) ≥ 0.05, extrapolated with q = a + b/√R over R ∈ {8, 16, 32}.
  - *Inconclusive* otherwise.
  - Need at R = 8 is reported against the same 0.05.
- **Use:** E-pair-t's resolved share at k′ = 20 must be ≥ 0.30 with R ≤ 16.
- **Validity:** the use numbers stand only if E-pair-t's semi-synthetic failure rate is
  ≤ 0.10.

**Recommendations, fixed now:**
- Need fails: stop. A second dataset would have to be argued to differ systematically.
- Need passes or is inconclusive, and use passes: go on to a second dataset.
- Need passes and use fails: the need–power tension holds on real data. No further data
  until the method is reframed.

**Deviation.** The plan states G3-need at native R. The pilot measures R ≤ 32 directly and
reaches native R only by extrapolation.

### Amendment 3 (2026-09-29): confirmatory analysis for the empirical note

**When and why.** Written after the reframe's prior-art check (`docs/reframe_prior_art.md`)
ended in its empirical-note branch, when the owner asked for the confirmatory analysis to
be written and pre-registered.

**What is registered.**
- The protocol: `docs/note_protocol.md`.
- The code that implements it: `run_things_confirmatory.py`, `analyze_things_confirmatory.py`,
  and the additions to `reach/things.py`.
- Its synthetic tests: `tests/test_things_confirmatory.py`, with its log in
  `results/things_confirmatory/`.

All of it is committed together, before the analysis touches the THINGS-EEG2 files.

**Observation from the synthetic check, made before any data.** With a known truth, the
reference disagreement moved toward the exact false comfort as the reference grew. The
extrapolant linear in 1/R_h landed closer to the exact value in most cells. The one linear
in R_h^(−1/2) came out low. The protocol's conservative bound already takes the smaller of
the two, so this observation changes nothing in it.

### Amendment 4 (2026-09-29): second dataset for the note

The owner authorized downloading the Allen Brain Observatory two-photon natural-scenes
data, and asked for the analysis to be pre-registered. The protocol is
`docs/note_protocol_allen.md`. It was committed after only the experiment metadata had
been fetched, before any recording file was downloaded or any response value read. The
code that implements it is committed separately, before the analysis runs.

### Amendment 5 (2026-10-01): common checks against the truth

**When and why.** The owner chose the Journal of Neuroscience Methods for a full manuscript
and asked for a comparison with existing methods, pre-registered like the earlier
analyses. The protocol is `docs/checks_protocol.md`. It measures six common checks
against the truth, on the confirmatory maps of both datasets:
- per-point trustworthiness;
- scDEED's "trustworthy" label;
- cross-validated distances;
- recurrence across t-SNE seeds, bootstrap resamples and split halves.

**What is registered, committed together before the analysis touches a data file:**
- the protocol;
- its code: `reach/checks.py`, `run_checks.py` and `analyze_checks.py`;
- its synthetic tests: `tests/test_checks.py`, with its log in `results/checks/`.

**Observation from the synthetic tests, made before any data.** scDEED's null was first
written to permute the 20 coordinates of the representation. That null kept the signal's
unequal spread across coordinates, and scDEED then called only 11% of points trustworthy on
strongly structured data. Permuting the measured features before the fixed representation,
as scDEED permutes genes before its PCA step, gave 74%. The protocol registers the
measured-feature null.

### Amendment 6 (2026-10-04): validating the estimator, in response to a colleague's review

**When and why.** On 2026-10-04 a colleague sent a revision plan for the JNM manuscript
(`docs/review_comments.md`). It asks that the estimator be validated end to end before the
paper decides what to claim. The audit of that plan (`docs/review_audit.md`) confirmed the
concern from the committed rows (`analyze_review_audit.py`):
- the estimator had been checked at one budget and for one map only;
- its larger extrapolant overstates at 16 and 32 trials;
- the mouse references are saturated;
- the subject-bootstrap bound is under-covered for a population median.

The protocol is `docs/revision_protocol.md`. It specifies:
- a simulation study in three layers;
- observable statuses for every result;
- reporting rules fixed before any new result;
- descriptive revision analyses on the real data.

**What is registered, committed together before any of its analyses run:**
- the protocol;
- the audit script;
- the correctness tests: `tests/test_estimand.py`, with its log in `results/revision/`.

A dated addendum, fixing what the development pilot decides, is committed before any
validation seed is used. The original claims keep their recorded verdicts. The manuscript
at tag `jnm-baseline` is the archived baseline.
