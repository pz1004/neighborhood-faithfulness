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
