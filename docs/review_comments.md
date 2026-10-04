# Concrete scientific revision plan

**Revise the paper around a validated neighborhood-reliability assessment—not around preserving the claim that most neighborhoods are spurious.** The first decision is whether the extrapolation and uncertainty procedure support latent-error estimates in the settings studied. The manuscript should be rewritten only after that question is resolved.

:chatgpt-content-reference{index="29"}[Download the complete revision plan](sandbox:/mnt/data/JNM_scientific_revision_plan.md)

**Proposed working title**

> **Assessing neighborhood reliability in maps of trial-averaged neural responses using held-out repetitions**

The experimental settings below are proposed additions, not existing results. They preserve the paper’s neighborhood definition, F/U/H structure, and two datasets while concentrating new work on the unresolved scientific issues.

## 1. Revision sequence

| Stage | Task | Required output | Decision |
|---|---|---|---|
| **WP0** | Freeze the current analysis and specify the estimand | Revision protocol, estimand specification, correctness tests | Establish exactly what is being estimated |
| **WP1** | Validate the estimator and complete uncertainty procedure | Calibration results and identified failure regimes | Determine which latent-error claims are supportable |
| **WP2** | Audit preprocessing and reference dependence | Trial-provenance record and controlled sensitivity analyses | Determine which reference constructions are defensible |
| **Gate A** | Select the scientific direction | A latent-error method, model-dependent estimation study, or finite-reference benchmark | Decide before expanding the benchmark |
| **WP3** | Rebuild comparisons around error, acceptance, and yield | Revised benchmark and reference-depth results | Establish practical utility |
| **WP4** | Rewrite the scientific argument and results presentation | Revised manuscript and supplement | Align every conclusion with its evidence |
| **Gate B** | Conduct a final claim-to-evidence audit | Completed scientific checklist | Determine submission readiness |

**Do not add another dataset or a large collection of embedding algorithms before Gate A.** Those additions would not resolve the central estimator and reference-construction questions.

## 2. WP0: Freeze the baseline and specify the estimand

### 2.1. Preserve the current analysis as an auditable baseline

Keep the existing analyses and outputs unchanged. Write a revision-specific protocol that separates original prespecified analyses, revision analyses designed after seeing the original results, and subsequent exploratory analyses.

Preserve the unsupported C2 result and the disclosure that the EEG analysis followed a pilot using the same recordings. New splits or seeds do not turn those recordings into an independent confirmation. The manuscript already documents both points. manuscript manuscript

### 2.2. Distinguish map-level, protocol-level, and subject-level quantities

For a fixed map \(M\), let \(A_i\) indicate that check \(X\) accepts neighborhood \(i\), and let \(T_i\) denote its latent faithfulness. Define

\[
\operatorname{FC}_X(M)
=
\frac{\sum_i A_i(1-T_i)}
{\sum_i A_i},
\]

provided at least one neighborhood is accepted. With zero acceptance, false comfort is **undefined**, not zero.

The manuscript currently pools counts across splits and seeds within subjects and then reports medians across subjects. Those summaries are not automatically the false comfort of one displayed map. manuscript

For the pooled procedure, explicitly define a protocol-level target such as

\[
\theta_{sX}(R)
=
\frac{
\mathbb E_Z\!\left[\sum_i A_{siZ}(1-T_{siZ})\right]
}{
\mathbb E_Z\!\left[\sum_i A_{siZ}\right]
},
\]

where \(Z\) includes the declared trial-sampling, split, representation-fitting, and map-seed mechanism. This ratio weights maps according to how many neighborhoods they accept; it generally differs from averaging their individual false-comfort ratios.

Then distinguish the descriptive median for the ten observed subjects from inference about a population median. State which one each interval addresses.

Also define the latent means **within the representation fitted on F for the relevant split**, rather than as a representation-independent neural geometry. The current representation is refitted on F and then applied to U and H. manuscript

### 2.3. Add inexpensive correctness tests before new experiments

Test zero-noise data, exact references, zero acceptance, all acceptance, deterministic ties, and joint-count identities.

In a zero-noise plug-in analysis, \(N=T\), so false comfort must be zero whenever neighborhoods are accepted. This does **not** imply that every neighborhood in a two-dimensional map is faithful; embedding distortion can remain.

Verify within each subject, before taking medians, that

\[
\text{acceptance}
=
\text{faithful yield}
+
\text{unfaithful yield}.
\]

Check that additional reference sub-blocks do not inadvertently change the weighting of maps, and do not treat splits, seeds, or overlapping reference sizes as independent biological observations.

**Deliverables:** `revision_protocol.md`, `estimand_specification.md`, and `unit_test_report.md`.

## 3. WP1: Validate the estimator and its uncertainty

This is the highest-priority work package. The current validation uses 20 draws per subject for the primary estimator and 10 for the comparison, with emphasis on median bias. The revision needs to evaluate the complete procedure used to support the reported bounds. manuscript

Organize the simulation protocol around aims, data-generating mechanisms, estimands, methods, and performance measures—the framework described by Morris, White, and Crowther [1]. The application-specific designs below are proposed for this paper. [Wiley Online Library](https://onlinelibrary.wiley.com/doi/10.1002/sim.8086?utm_source=chatgpt.com)

### 3.1. Use three validation layers

**Layer A: Conditional reference validation.** Hold the known mean responses, fitted representation, analyst data, displayed maps, and acceptance decisions fixed. Regenerate only the held-out references. Compare the extrapolated estimate against exact labels.

This isolates the reference-correction problem without repeatedly fitting every map.

**Layer B: End-to-end point estimation.** Generate new F, U, and H data, fit the representation and maps, and compute exact latent labels for those same maps. Transform the known means through each simulated F-fitted representation when calculating the truth.

This tests whether success in the fixed-map experiment survives the actual analysis pipeline.

**Layer C: Study-level uncertainty.** Generate independent virtual studies with ten subjects from an explicitly specified heterogeneous population. Reproduce the published interval procedure, including the split/seed design, both extrapolants, subject bootstrap, zero-denominator handling, and any retained guard.

First evaluate the subject-median interval using true subject-level protocol quantities; then evaluate it using estimated quantities. Where those true quantities require oracle simulation, compute them independently and sufficiently precisely. A realized ratio from a finite collection of maps is not automatically a subject’s protocol expectation.

**Do not assess a population-median interval against only the realized sample median.** That would validate a different inferential target.

### 3.2. Use a small, informative simulation grid

Begin at \(R=8\), retaining the present reference grids:

\[
R_h^{\mathrm{EEG}}=\{8,20,40\},
\qquad
R_h^{\mathrm{mouse}}=\{8,12,24\}.
\]

These are the current experimental reference depths. manuscript

| Component | Proposed settings | Purpose |
|---|---|---|
| **Truth family 1** | Original plasmode truth based on the H mean | Reproduce the existing validation |
| **Truth family 2** | Low-dimensional, well-separated means embedded in the same ambient representation | Establish a favorable recovery regime |
| **Truth family 3** | Means with small distance gaps around the \(k'\)-neighbor cutoff | Test difficult neighborhood ordering |
| **Noise scale** | Multipliers \(\{0.1,0.5,1,2\}\), plus a separate zero-noise control | Span substantially different signal-to-noise regimes |
| **Dependence** | Independent residual resampling first; one targeted dependent-noise extension informed by WP2 | Test the relevant assumption without an excessive factorial design |

This gives **12 nonzero-noise screening cells per dataset-sized design**. During a separate development pilot, verify that they actually span low, intermediate, and high exact false comfort. Adjust the signal/noise scale before freezing the validation grid when necessary. Include cases near the 25% and 50% decision thresholds.

Add one particularly important component stress test: **known faithful displayed neighborhoods evaluated against shallow, noisy references**. Determine whether a high disagreement plateau produces an incorrectly confident high-false-comfort estimate. Label this as a conditional-reference diagnostic, not as an estimate of how often that scenario occurs in the real datasets.

### 3.3. Set replication counts and computational priorities

Use **100 independent repetitions per screening cell**. For the study-level scenarios supporting retained inferential claims, plan **1,000 independent virtual studies**, starting with the primary plug-in analysis.

At 1,000 studies, the Monte Carlo standard error is approximately

\[
\sqrt{0.95(0.05)/1000}=0.0069
\]

for 95% coverage, and approximately 0.0049 for a 2.5% tail-error probability. Report Monte Carlo intervals alongside estimated coverage and error rates.

Time the first 20 runs before allocating the full computational budget. Cache maps, acceptance masks, and reference labels for reuse across extrapolants and summaries. Avoid running 50-map bootstrap recurrence in every screening condition; add its expensive calibration only where a latent-error bound for that check will be retained.

### 3.4. Replace “small median bias” with a complete performance assessment

Report bias, mean absolute error, root-mean-square error, an upper quantile of absolute error, interval width, coverage, erroneous threshold declarations, out-of-range estimates, and the rate of non-evaluable outputs.

For the current nominal tails, explicitly evaluate

\[
\Pr(L>\Theta),\qquad \Pr(U<\Theta),
\]

against 0.025, and joint interval coverage against 0.95. Also report how frequently the procedure makes an informative declaration; an always-uninformative interval is not a useful method.

As a proposed **screening rule**, absolute bias above five percentage points or RMSE above ten percentage points should trigger investigation. These tolerances are not evidence of calibration, and much smaller errors can invalidate a claim near its threshold.

Evaluate the existing five-percentage-point guard as part of the archived baseline. Any replacement correction must be developed on one simulation set and assessed on independent validation simulations. **Do not tune a correction using the validation truths against which it is subsequently judged.**

### 3.5. Add reference-depth diagnostics and explicit unresolved outputs

Fit the two shallower sizes and predict the largest observed size:

\[
\{8,20\}\rightarrow40 \quad\text{for EEG},
\qquad
\{8,12\}\rightarrow24 \quad\text{for mouse}.
\]

In simulations, additionally generate references at twice and four times the current maximum depth. Fit only the original depths, then compare the extrapolation against those deeper references and the exact truth.

For real data, call this a **finite-depth predictive check**, not validation of the infinite-reference limit. The nearly flat mouse curves in Figure 1d are precisely why this distinction matters: a ceiling can persist over observed depths while providing limited information about the ultimate limit. manuscript

Return explicit estimation statuses, for example:

- **Estimable:** the procedure operates within its empirically supported domain.
- **Insufficient acceptance:** too little information supports the conditional error estimate.
- **Reference extrapolation unresolved:** the available reference does not support the latent estimate.

Near-ceiling curves, divergent extrapolants, impossible joint probabilities, and poor deeper-reference prediction are candidate warning criteria. Fix any reporting rule using development data and validate performance both overall and conditional on reporting. An applicability domain cannot depend exclusively on true rank margins that users cannot observe.

**Deliverables:** subject-level validation results, study-level calibration results, reference-depth diagnostics, and a report showing both successful and unsuccessful regimes.

## 4. WP2: Audit preprocessing and reference dependence

### 4.1. Reconstruct the complete trial provenance

Create a trial manifest containing subject, stimulus, original session/run, onset time, feature-window support, preprocessing identifier, split, and F/U/H assignment. Record which observations estimated each shared transform or normalization.

The supplement currently reconstructs EEG session membership, not complete chronological trial identity. Establish whether raw event information can be recovered before claiming that temporal purging has been implemented. supplementary

### 4.2. EEG: separate preprocessing and session-dependence sensitivities

First reproduce the current analysis. Then reconstruct inputs before shared noise normalization and apply a common transformation fitted only on F or an explicitly independent calibration block. Apply it unchanged to U and H.

Do not normalize U and H independently into different coordinate systems.

For a concrete session-disjoint sensitivity at \(R=8\) and 16, use:

| Block | Proposed allocation |
|---|---|
| **F** | Eight repetitions per image from one session |
| **U** | Up to sixteen repetitions per image from another session |
| **H** | Forty repetitions per image from the remaining two sessions |

Rotate the ordered F/U session choices, giving twelve possible assignments. Keep all rotations clustered within subject. This allocation uses the reported four sessions and twenty repetitions per image per session; actual trial integrity still needs verification. manuscript supplementary

Compare against a budget-matched random-split analysis with the same preprocessing. Report the session mixture defining each response target. Across-session changes can reflect state or signal changes as well as dependence, so interpret this as a sensitivity analysis—not as a pure causal measurement of leakage.

**This design does not validate the \(R=32\) result.** Either implement a separately validated lower-reference-depth design for that budget or explicitly restrict the session-disjoint evidence to \(R=8\) and 16.

### 4.3. Mouse V1: split the continuous timeline, not each image separately

The manuscript uses a 0.5-second response window for presentations separated by 0.25 seconds, so window overlap must be addressed explicitly. manuscript

Assign contiguous acquisition blocks globally across stimuli. Remove trials whose complete feature support crosses a split boundary. Add a guard interval informed by residual temporal dependence estimated from F/calibration data, not selected to obtain a preferred result.

Use timestamps and stimulus counts—without examining disagreement outcomes—to determine feasible common trial budgets after purging. Attempt \(R=8\) first. When the original F/U/H counts cannot be maintained, lower the common budgets and rerun matched controls rather than silently comparing unequal information budgets.

Report direct window overlap, residual dependence, retained repetitions, finite-reference disagreement, and any supported latent estimates. Reduced overlap improves the design but does not prove complete independence.

**Deliverables:** trial manifest, preprocessing-dependency table, and paired subject-level sensitivity results.

## 5. Gate A: Decide which scientific claim to retain

| Validation outcome | Revision direction |
|---|---|
| Point estimation and uncertainty perform adequately in a declared domain, and the reference audit supports the empirical interpretation | **Proceed with a scoped latent-error methods paper.** |
| Point estimation remains informative, but the uncertainty procedure is not calibrated | **Retain explicitly model-dependent estimates and sensitivity results; withdraw conservative-bound language and the affected threshold declarations.** |
| Infinite-reference extrapolation is unreliable, but finite-reference disagreement is reproducible | **Reframe as a finite-reference benchmarking paper with a different estimand.** |
| Neither the estimator nor the controlled empirical comparison supports a distinctive, actionable conclusion | **Stop expanding the current manuscript.** |

A dataset-specific decision is acceptable. For example, latent estimates might be supported for some EEG settings while mouse results remain finite-reference measurements. Those targets must be labeled separately.

Under the fallback direction, define the target as disagreement with a reference of stated depth \(r\). **Do not rename finite-reference disagreement “latent false comfort” or treat the deepest available mean as ground truth.** The manuscript itself acknowledges that finite reference labels remain noisy. manuscript

## 6. WP3: Rebuild the benchmark around acceptance and yield

### 6.1. Preserve native operating points and add matched-acceptance comparisons

Retain the six checks and plug-in comparator at their existing operating points. Add descriptive comparisons at proposed acceptance levels of **5%, 10%, 20%, 30%, and 50%**, where their scores permit them.

Use the existing score, reach, or recurrence fraction. Fix score directions and tie handling without consulting H labels. Do not force a smooth curve on a method with only a few attainable operating points.

For each subject, jointly calculate

\[
A_X=\Pr(X\text{ accepts}),
\]

\[
Y_X=\Pr(X\text{ accepts},T=1),
\qquad
U_X=\Pr(X\text{ accepts},T=0),
\]

with

\[
A_X=Y_X+U_X,\qquad
\operatorname{FC}_X=U_X/A_X
\]

when \(A_X>0\).

The supplement already contains acceptance and yield results, so this is partly a reorganization of existing evidence rather than an entirely new experiment. Compute joint quantities before taking medians; do not multiply separately reported medians. supplementary

For zero acceptance, report zero yield and undefined false comfort. For rare acceptance, expose the denominator, contributing subjects, and estimation status.

### 6.2. Add a prevalence-aware comparator

For each map, randomly select exactly the same number of neighborhoods as the evaluated check. This preserves map-specific acceptance counts and the weighting they induce.

Compare faithful yield and enrichment, not only absolute false comfort. When faithful neighborhoods are extremely rare, many checks can have false comfort near one despite differing in discrimination. The mouse results make this issue particularly relevant. manuscript manuscript

Do not use a single pooled prevalence as the comparator when both faithfulness and acceptance vary between maps.

### 6.3. Audit the implemented comparators

For **scDEED**, document the fixed-representation adaptation, permutation count, null pooling, score calculation, and implementation tests. Avoid selecting a null procedure solely because it labels more structured points trustworthy. The current adaptation is described in Supplementary Section S1.4. supplementary

For **cv-MDS**, report negative squared-distance frequency, zero-distance ties after truncation, the ranking convention, and MDS fit. Restrict conclusions to the implemented transformed-distance pipeline, and keep them separate from the cross-validated-distance check on the t-SNE map. manuscript manuscript

Finally, report the total repetition requirement beside the analyst’s \(R\). The held-out assessment uses additional data; it should not be presented as a cost-matched algorithmic improvement over checks restricted to U.

## 7. WP4: Rewrite the theory and interpretation

### 7.1. Replace the impossibility argument with a limited convergence statement

Retain the Gaussian algebraic identities where useful, but remove the conclusion that they rule out informative same-trial inference. The manuscript currently makes that stronger claim while also acknowledging the possibility of per-point certificates. manuscript

A suitable replacement theoretical result is:

> For a fixed displayed map, consistent reference-distance estimates recover the true top-\(k'\) neighbor sets when the relevant cutoff gaps are strictly positive. Consequently, the corresponding reach labels converge to their latent labels.

A short proof can use the gap

\[
\Delta_i=d_{i,(k'+1)}-d_{i,(k')}>0.
\]

With finitely many objects, sufficiently small uniform distance error preserves cutoff membership and hence the reach label.

**State explicitly what this does not establish:** a convergence rate, either proposed finite-depth linear extrapolation, or confidence-interval coverage. Exact ties require separate treatment; deterministic index-based tie breaking alone does not establish convergence of noisy rankings.

### 7.2. Make these section-level changes

| Location | Required revision |
|---|---|
| **Title and abstract** | Remove unsupported universal wording; distinguish analyst and reference budgets; identify model-dependent estimates as such. |
| **Introduction** | Emphasize the assessment problem, calibration, and practical limits rather than novelty of the generic conditional error measure. |
| **Sections 2.1–2.3** | Add precise estimands, assumptions, layered validation, and unresolved-output rules. |
| **Sections 2.4 and 2.6** | Document end-to-end preprocessing and the dependence sensitivities. |
| **Section 2.7** | Separate original prespecification from revision analyses; replace broad “generality” claims with within-study robustness. |
| **Sections 2.8 and 3.5** | Separate false comfort, acceptance, abstention, faithful yield, and enrichment. |
| **Section 3.1** | Present calibration and failure regimes, not only median bias. |
| **Section 3.6** | Treat the neighbor-distance/noise ratio as descriptive rather than an established mechanism. |
| **Section 4.1** | Remove the claim that dissimilarity reliability or an RSA noise ceiling certifies the displayed map’s global geometry. |
| **Section 4.2.4 and S1.3** | Replace the universal impossibility interpretation with the limited statement above. |
| **Conclusion** | State supported settings and practical trade-offs without requiring the original percentages to survive. |

Keep the existing structured abstract headings—**Background, New method, Results, Comparison with existing methods, Conclusions**—and write its Results portion last.

### 7.3. Resolve the remaining factual and logical inconsistencies

Correct the Figure 2a caption: the EEG cv-MDS row does not show more accepted than estimated faithful neighborhoods. Correct the statement that four checks accept at least as many neighborhoods as the plug-in check. State that false-comfort percentages come from joint labels rather than implying that they follow directly from marginal percentages. manuscript manuscript manuscript manuscript

Also remove the inference from positive cross-validated point estimates to statistically established pairwise distinguishability. Describe the plasmode truth as exact within the stipulated simulation, while discussing its limited biological representativeness. These are different issues in the current Results and Limitations. manuscript manuscript

For the pilot certificates, report how often certificates were issued in simulations, not only that no false statements occurred. Remove the claim that multiplicity correction was not responsible unless an analysis isolates that effect. supplementary

## 8. Revised figures, schedule, and completion criteria

### Results presentation

Replace rather than simply add to the current figure set. Aim for five main figures:

| Figure | Content |
|---|---|
| **1** | Estimand and F/U/H workflow, retaining the reach schematic and showing total repetition requirements |
| **2** | Estimator calibration: exact versus estimated false comfort, coverage, and an identified failure regime |
| **3** | Observed reference-depth curves, with extrapolation shown only where supported |
| **4** | Acceptance-versus-error curves, native operating points, and random-selection baselines |
| **5** | One representative neural-response map illustrating observed neighborhood evidence without claiming per-point certification |

Add a continuous overlap summary,

\[
\frac{|E_i\cap N_{k'}(\text{reference},i)|}{k},
\]

to distinguish slight failures of the binary criterion from near-total disagreement. On real data, label this as a finite-reference measure unless a latent interpretation is separately validated.

Move the large configuration robustness display and full subject-level sensitivity results to the supplement.

**The revision should establish when this assessment supports a neighborhood interpretation and when it should decline to do so. A smaller, calibrated effect would support a stronger paper than the current larger but insufficiently calibrated headline.**

### Methodological reference

[1] T. P. Morris, I. R. White, and M. J. Crowther, “Using simulation studies to evaluate statistical methods,” *Statistics in Medicine*, vol. 38, no. 11, pp. 2074–2102, 2019, doi: 10.1002/sim.8086. [pubmed.ncbi.nlm.nih.gov](https://pubmed.ncbi.nlm.nih.gov/30652356/?utm_source=chatgpt.com)
