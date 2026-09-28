# Prior-art check: a neighborhood noise ceiling for 2-D displays

2026-09-29. This is the first step of the reframing that the THINGS-EEG2 pilot's
pre-registered recommendation requires (plan Amendment 2; `docs/gate3_pilot_things.md`,
Section 6). Section 1 was written and committed before any search result was seen. The
searches' findings and the verdict follow in Section 2.

## 1. Protocol (fixed before the searches reported)

### 1.1 The reframed item

The per-point certificate is replaced by a dataset-level measurement. Given replicate
measurements, the method estimates how wrong a 2-D display's neighborhoods are and how much
of their apparent fidelity is spurious. It also estimates how both change with the number
of replicates.

| id | claim |
|---|---|
| NC1 | **Target.** Two dataset-level quantities, per display and tolerance, stated against the noise-free truth. (a) The *mislabel rate* of the naive fidelity map, which labels a point faithful if at least (1−ε)k of its k displayed neighbors are among its k′ nearest in the noisy data. (b) *False comfort*: the share of the naive map's faithful calls that are not faithful in truth. |
| NC2 | **Estimator.** From replicate splits. User labels from R replicates are compared with plug-in labels from disjoint held-out reference sub-blocks of increasing size R_h. The disagreement is extrapolated to an infinite reference, linearly in R_h^(−1/2), with uncertainty over subjects and splits. |
| NC3 | **Neighborhood noise ceiling.** The true fidelity a display of R-replicate means can reach, and its growth in R, estimated the same way and used to plan replicate budgets. |
| NC4 | **Validation.** A semi-synthetic check in which the truth is the mean of held-out replicates and the noise is resampled within-object residuals. |
| NC5 | **Empirical finding.** On THINGS-EEG2, most apparently faithful t-SNE neighborhoods of condition means are spurious at typical budgets: 84% of naive faithful calls are wrong at R = 8, and 47% at R = 32. |

### 1.2 Searches

Three were run in parallel. Each had about 30 queries and was asked to open every source it
cites.
1. Neuroscience and representational-space reliability: RSA noise ceilings,
   noise-corrected model comparison, and neighbor-level alignment measures.
2. Statistics of measurement error: SIMEX and extrapolation to zero noise, imperfect
   reference standards, reproducibility of rankings and top-k sets, and nearest neighbors
   under error.
3. Embedding evaluation and single-cell reproducibility: DR quality metrics against the
   truth, and replicate-, downsampling- or thinning-based reliability of kNN graphs and
   embeddings.

### 1.3 Decision rule

- **Drop the reframe** if a single prior work estimates, from replicates or resampling, how
  often neighborhoods computed from noisy data (a kNN graph, an embedding or a display)
  disagree with the noise-free neighborhoods, and corrects for the reference's own noise.
  That is NC1 together with NC2 in substance.
- **Continue, with narrowed claims,** if only components are found. Examples: whole-RDM
  noise ceilings, SIMEX in other settings, or kNN reproducibility between replicates with
  no truth target.
- **A short empirical note** is the recommendation if NC1–NC4 are each anticipated, so that
  their combination is routine and only the finding NC5 remains new. It replaces a methods
  paper.

## 2. Findings and verdict

The three searches each ran 28–30 queries and opened most sources through publisher,
Europe PMC, arXiv or CrossRef records. Some were read at abstract level only, and each
search marked what it did not open.

### 2.1 Verdicts, consolidated

| claim | verdict | closest work | what it lacks |
|---|---|---|---|
| NC1 target | partly anticipated | **Psychometrics** already defines this kind of target: the classification accuracy and consistency of decisions made from noisy scores, against true scores, split into false positives and negatives, and estimated without the truth (Livingston & Lewis 1995; Rudner 2001, 2005). **Dimensionality reduction** has the framing of judging a display against the signal, not the noisy data (Lin & Fukuyama 2024; DEMaP, Moon et al. 2019). | Set-valued, display-dependent neighborhood labels. And false comfort: a self-check computed on the same noise as the display. Its mechanism is known, as conditional dependence between test and reference errors (Vacek 1985; Flegal et al. 1991). |
| NC2 estimator | partly anticipated; in substance a SIMEX-type extrapolation | Simulation extrapolation to zero error (Cook & Stefanski 1994); its replicate-based form (Devanarayan & Stefanski 2002); misclassification SIMEX (Küchenhoff et al. 2006); Spearman–Brown step-up of half-test pass/fail consistency (Woodruff & Sawyer 1989). Ceilings extrapolated to unlimited data or subjects (David & Gallant 2005; Schrimpf et al. 2021). Held-out repetition halves as a reference (Gifford et al. 2022). kNN overlap against a deep-sequencing reference, uncorrected (Ahlmann-Eltze & Huber 2023). | Its use for neighborhood labels, and the extrapolant q + b·R_h^(−1/2) for a disagreement rate. |
| NC3 noise ceiling in R | partly anticipated | The dependence of reliability and decision accuracy on the number of measurements is standard: Spearman–Brown; effective test length; the per-voxel noise ceiling of NSD, NC(n) = 100·ncsnr²/(ncsnr² + 1/n) (Allen et al. 2022); sample sizes for correct ranking (Bechhofer 1954). | Not found for neighborhoods or displays. Every noise ceiling found concerns whole-RDM correlation, single voxels or sensors, or decoding. |
| NC4 validation | anticipated, as a design | Plasmode simulation (Schreck et al.); deep data with binomial downsampling (Ahlmann-Eltze & Huber 2023); residual bootstrap for RDM reliability (Raju 2026) | Its use to validate display labels, which is routine |
| NC5 finding | **not found** | Neighborhood loss of 2-D maps is well documented against the noisy data (Chari & Pachter 2023; Cooley et al. 2019; Lause, Berens & Kobak 2024). | Nobody measured it against the truth, for the naive map's own faithful calls, as a function of replicates, or on EEG. |

### 2.2 Decision under the rule of Section 1.3

**The drop trigger is not met.** No work estimates how often noisy neighborhoods disagree
with the noise-free ones while correcting for the reference's own noise. The closest,
Ahlmann-Eltze & Huber 2023, treats its deep reference as the truth.

**Branch: a short empirical note, not a methods paper.** Each of NC1–NC4 has a direct
precedent in method type:
- the target is classification accuracy against true scores;
- the estimator is a SIMEX-type extrapolation, done by growing the reference instead of
  adding noise;
- the ceiling's growth in R is Spearman–Brown logic;
- the validation is a plasmode simulation.

Their combination for neighborhood labels is new as an application but routine as a
method. What remains new is the finding NC5, together with one conceptual point: a
fidelity check computed from the same noise as the display is biased toward false comfort.

**Judgment call.** In a visualization or dimensionality-reduction venue the combination
could be presented as a methods-and-finding paper. But a reviewer who knows psychometrics
or measurement error would read NC1–NC3 as "classification accuracy with a SIMEX-type
correction, for kNN labels". The honest claim is the finding, with the estimator described
as an application of those methods.

### 2.3 Threats and required citations

- **Lin & Fukuyama 2024** was the most dangerous work for two of the three searches. It
  says to judge t-SNE/UMAP against the signal, not the data. The differences to state:
  - a truth from held-out replicates instead of PCA truncation;
  - per-point labels;
  - mislabel and false-comfort rates;
  - dependence on R and planning.
- **Livingston & Lewis 1995** was the most dangerous for the statistics search: it defines
  classification accuracy against true scores.
- **The SIMEX family.** The estimator should be named as a SIMEX-type extrapolation and
  cited as such.
- **Ahlmann-Eltze & Huber 2023** is the design to contrast: a deep reference without
  correction, and no display or R series.
- **The RSA objection** (Walther et al. 2016; Guggenmos et al. 2018). Practitioners use
  cross-validated, noise-normalized distances, so the headline has to hold for naive maps
  built from those too.
- **GSN** (Kay et al. 2025) estimates signal and noise covariances. Simulating R-trial
  displays from it is a parametric alternative to benchmark against.

### 2.4 The estimator's failure modes, and what they do to the pilot

- **The extrapolant.** The statistics search argued that the reference term can be O(1/R_h)
  instead of O(R_h^(−1/2)). This happens when reference noise flips labels symmetrically
  around a continuous margin. A fit linear in R_h^(−1/2) then puts the intercept too low.
  The pilot refitted with the 1/R_h form (medians over subjects, t-SNE, k′ = 20):

  | user R | linear in R_h^(−1/2) (the pilot's primary) | linear in 1/R_h |
  |---|---|---|
  | 8 | 0.172 | 0.172 |
  | 16 | 0.194 | 0.216 |
  | 32 | 0.195 | 0.231 |

  The pilot's estimator was the lower of the two, so the need verdict is unaffected. Future
  estimates should report both, or estimate the exponent (Bertail, Politis & Romano 1999).
- **Variance and range.** Extrapolation variance is high and depends on the available
  range of R_h (Lage-Castellanos et al. 2019). The unit of uncertainty is the subject.
- **Shared error.** Error shared by user and reference becomes part of the "truth".
  Session-stratified blocks cancel the large global session offsets. The image-specific
  session interaction is small (F = 1.14) but not zero.
- **Ties.** Resampling rank quantities undercovers near ties (Hall & Miller 2009).

### 2.5 What an empirical note needs before it is written

A confirmatory, pre-registered analysis, using the reserved seed range ≥ 4_000_000:
1. **False comfort on the real data**, estimated with the reference method. The 84% and 47%
   figures so far are exploratory and semi-synthetic.
2. **The RSA objection**: naive maps from cross-validated, noise-normalized distances as
   well as from plug-in distances.
3. **Both extrapolants**, with subject-level bootstrap intervals.
4. **Robustness**:
   - k ∈ {5, 10, 20} and ε ∈ {0.1, 0.2, 0.3};
   - t-SNE at perplexities 15, 30 and 50, plus PCA-2 and MDS;
   - two time windows;
   - an unwhitened representation.
5. **Preferably a second replicate-rich dataset** in another modality, which needs a
   download and the owner's permission.
