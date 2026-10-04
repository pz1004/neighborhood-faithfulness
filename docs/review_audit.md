# Audit of the colleague's revision plan (2026-10-04)

**What was audited.**
- `docs/review_comments.md`, as received on 2026-10-04. It is a revision plan in five work packages (WP0–WP4) with two decision gates, rather than a list of point comments.
  - It was generated with ChatGPT. Its link to a "complete revision plan" points into a ChatGPT sandbox and resolves nowhere, so only the text present was audited.
  - Its bare "manuscript" and "supplementary" tags are citation markers that lost their targets. Every claim was therefore located in our sources here, with line numbers.
- The manuscript audited is the one at tag `jnm-baseline` (commit `43d636c`).

**How.**
- Every factual claim of the review was checked against the manuscript, the supplement, the protocols and the code.
- Every number below was recomputed from the committed result rows by `analyze_review_audit.py`, using the analysis scripts' own estimators.
  - Its output is `results/review_audit/audit.txt` and `audit_summary.json`.
  - It reproduces the published R = 8 plasmode bias exactly (−0.008 EEG, +0.003 mouse).
  - Nothing was simulated, and no data file was opened.

**Verdict on the whole.**
- **The central concern is right.** The estimator was validated too narrowly to carry the claims made with it, and the committed results show this more sharply than the review could (Section 2).
- **The errors it names are real.** The six factual or logical errors in the text that it identifies are confirmed.
- **We take a smaller route in three places,** where it asks for more than the problem needs: re-processing the raw EEG, re-splitting the mouse timeline, and a representative map as a main figure (Section 3).
- **What is not yet known.** The direction it proposes, deciding the paper's claims only after an end-to-end validation, is adopted. The validation and the revision analyses are specified in `docs/revision_protocol.md` (plan Amendment 6). Which claims survive is not known until they run.

## 1. What the review asks for

**Its thesis.** Revise the paper around a validated assessment of neighborhood faithfulness, not around defending the headline that most faithful-looking neighborhoods are spurious.

**Its sequence:**
- **WP0:** freeze the baseline and define the estimands.
- **WP1:** validate the estimator and its intervals in three layers.
- **WP2:** audit the trial provenance and dependence between blocks.
- **Gate A:** choose the scientific direction.
- **WP3:** rebuild the benchmark around acceptance and yield.
- **WP4:** rewrite the paper.
- **Gate B:** a final audit that matches every claim to its evidence.

## 2. Evidence recomputed from the committed results

Below, "sqrt" is the extrapolant linear in R_h^−1/2 and "inv" the one linear in 1/R_h. The published "smaller" estimate is almost always sqrt.

**E1. The plasmode was analyzed at one budget only.**
- The simulation rows for R = 16 and 32 exist but were never analyzed.
- The table gives the median over subjects of estimate minus exact false comfort, with the RMSE in brackets.

  | | R = 8 | R = 16 | R = 32 |
  |---|---|---|---|
  | EEG, sqrt | −0.008 (0.030) | +0.006 (0.043) | +0.026 (0.071) |
  | EEG, inv | +0.042 (0.046) | +0.098 (0.102) | +0.179 (0.179) |
  | mouse V1, sqrt | +0.003 (0.012) | +0.038 (0.045) | — |
  | mouse V1, inv | +0.004 (0.011) | +0.048 (0.050) | — |

- **Per-subject errors.** For sqrt at R = 32 they run from −0.077 to +0.135.
- **The mislabel rate q** is recovered well everywhere: the RMSE is ≤ 0.027.
- **Consequences:**
  - Every reported range of false comfort, "the range of the two extrapolants", takes its upper end from inv. At 16 and 32 trials that extrapolant overstates by 10 and 18 points (for example, 48–64% at 32 trials).
  - The checks' upper bound U also uses the larger extrapolant.

**E2. Mouse references are saturated.**
- **Every mouse cell.** The disagreement with the reference stays at about 1 for every reference size (8, 12, 24):
  - at R = 8 in the primary cell, the median disagreement is 1.000, 1.000 and 1.000 across the three sizes;
  - in every one of the 24 mouse cells and budgets, 7 to 10 of the 10 mice have a curve that ends at ≥ 0.95 after falling by ≤ 0.02.
- **The EEG unwhitened cell.** It is saturated in the same way (10 of 10 subjects at R = 8).
- **The other EEG cells.** Their curves fall by 0.04–0.38.
- **Why the mouse estimate is unresolved.** A random selection of as many neighborhoods per map gives ≥ 0.997 too (E4). So "essentially all, bound 99.7%" is the intercept of a flat line, not a resolved value.
- **Where the plasmode leaves room** (mouse, R = 16, exact false comfort 0.875–0.988), the estimator already overstates by 4–5 points, and by up to 9 in a single mouse.

**E3. Finite-depth predictive check.**
- **What it does:** a line through the two shallower reference sizes predicts the deepest.
- **Median error (predicted minus observed) of the sqrt extrapolant:**

  | | R = 8 | R = 16 | R = 32 |
  |---|---|---|---|
  | real EEG, primary cell | +0.011 | +0.023 | +0.040 |
  | EEG plasmode | +0.011 | +0.026 | +0.052 |

- **The plasmode reproduces the curvature** of the real EEG curves, budget by budget.
- **Prediction error tracks bias.** Within the plasmode, a larger prediction error goes with a larger upward extrapolation bias (E1), so the check is a usable warning signal.
- **It cannot detect saturation.** In mouse V1 the error is about 0 because the curves are flat.

**E4. Prevalence: a random selection matched per map.**
- **The comparator.** A random selection of as many neighborhoods as the check passes in each map would have false comfort 1 − (share faithful), weighted by each map's passes.
- **Enrichment** is the check's precision divided by the random selection's precision. EEG medians at 8 trials:

  | | False comfort | Random selection | Enrichment |
  |---|---|---|---|
  | plug-in, t-SNE | 0.81–0.86 | 0.93–0.95 | 2.8 |
  | plug-in, cv-MDS | 0.62–0.72 | 0.96–0.97 | 10.2 |

- **Enrichment of the checks in EEG at 8 trials:**

  | Check | Enrichment |
  |---|---|
  | trustworthiness | 3.1 |
  | scDEED | 2.2 |
  | cross-validated distances | 4.2 |
  | seed recurrence | 2.0 |
  | split-half recurrence | 5.7–6.1 |

- **Same-trial labels carry information.** They are better than chance, sometimes by a factor of ten. What makes most of their passes wrong is that faithful neighborhoods are rare.
- **Mouse V1.** Every check has a false comfort near 1 and enrichment cannot be determined, because the random selection is wrong 99.7–100% of the time.

**E5. The lower bound L with 10 subjects.**
- **What L is.** In EEG, L equals the 3rd-smallest subject estimate:

  | Claim | L | 2nd-smallest |
  |---|---|---|
  | C1 | 0.642 | 0.605 |
  | C1-RSA | 0.550 | 0.530 |
  | C2 | 0.231 | 0.215 |

- **Its coverage.** As a lower bound for a population median from 10 subjects, the 3rd-smallest estimate errs with probability 0.055, even if every subject's value were exact; the nominal level is 0.025. The 2nd-smallest errs with probability 0.011.
- **The EEG checks at 8 trials,** judged by the 2nd-smallest estimate:

  | Check | 2nd-smallest | Classification |
  |---|---|---|
  | trustworthiness | 0.591 | still "leaves most" |
  | scDEED | 0.688 | still "leaves most" |
  | cross-validated distances | 0.521 | still "leaves most" |
  | seed recurrence | 0.693 | still "leaves most" |
  | split-half recurrence | 0.311 to 0.795 (lower to upper) | stays undetermined |

**E6. Pooling over maps, and maps that pass nothing.**
- **Pooled ratio vs per-map mean.**
  - The pooled ratio of sums and the mean of per-map false comfort have medians over subjects within 0.003 of each other for every check.
  - Within a single subject they differ by at most 0.003–0.007 for most checks.
  - The exceptions: 0.020 for cross-validated distances, and 0.090 for split-half recurrence, where many maps pass nothing.
- **Maps that pass nothing** (R = 8, out of 300 maps per dataset):
  - bootstrap recurrence: 293 EEG maps and all 300 mouse maps;
  - split-half recurrence: 96 EEG maps (a median of 17% per subject, up to 93%) and 296 mouse maps;
  - cross-validated distances: 1 EEG map, and 44 mouse maps (a median of 13% per mouse, up to 37%).

**E7. The plasmode checked only the primary t-SNE map.** All four simulation files carry only that map.
- The cv-MDS estimator behind C1-RSA was never checked, and neither were the 11 variations.
- The bias guard applied to them used the t-SNE map's bias.

**E8. The pilot's certificates never resolved a neighborhood in simulation.**
- **Issuance:** in the 900 simulated draws that recorded issuance, no neighborhood was resolved at k′ = 10, 20 or 40.
- **Validity:** the 3,000 validity draws recorded only violations, and found none.
- **Consequence:** "made no false statement in 3000 draws" (S1.5) is therefore uninformative for reach claims. Only the coverage of the pairwise intervals is a real result.
- **The multiplicity claim is untested.** "Not the multiplicity correction" rests on an unscripted estimate in `docs/gate3_pilot_things.md` §5. S1.5 cites Table 4 for it, which does not test it.

## 3. Verdict on each comment

**Verdicts:**
- **Accept**: the concern and the proposed fix are both right.
- **Accept, changed**: the concern is right, but the fix differs.
- **Partly**: part of the concern holds.
- **Done**: already addressed.
- **Decline**: not adopted, for the reason given.

Line numbers refer to `manuscript/manuscript.tex` at `jnm-baseline` unless another file is named. "Step" refers to the revision plan in `docs/revision_protocol.md`.

| # | Comment (§ of the review) | Verdict | Why | Step |
|---|---|---|---|---|
| R1 | Revise around a validated assessment, not the "spurious" headline (top) | Accept | E1, E2 and E7: the validation covers one budget and one map, and the mouse estimate rests on a saturated reference | Gate A, rewrite |
| R2 | Working title "Assessing neighborhood *reliability*…" (top) | Accept, changed | Reliability means reproducibility. That is what recurrence checks measure, and the paper shows it is not faithfulness (§4.2.2). The title should say "faithfulness"; it is chosen at Gate B | rewrite |
| R3 | Work packages with two gates; no new datasets or embedding methods before Gate A (§1) | Accept | Adopted as the stages of the revision protocol | all |
| R4 | Freeze the baseline; separate original, revision and exploratory analyses (§2.1) | Accept | Baseline tagged `jnm-baseline`. The protocol also discloses that E1–E8 were seen before it was written | 0, protocol §1 |
| R5 | Keep C2's failure and the disclosure that the pilot reused the EEG recordings (§2.1) | Done | `:467`, `:393`, `:684` | — |
| R6 | False comfort of one map; undefined, not zero, when nothing is accepted (§2.2) | Accept | Already implemented, as NaN plus the not-evaluable rule (`analyze_checks.py:58, 138`), but the paper never states it | rewrite |
| R7 | The pooled ratio is a protocol-level target θ_s(R), not one map's false comfort (§2.2) | Accept | Defined in the protocol. Numerically it matters little (E6) | protocol §2 |
| R8 | Say whether intervals target the 10-subject median or a population median (§2.2) | Accept | E5 sharpens it: L is under-covered for a population median | Layer C |
| R9 | The truth lives in each split's representation fitted on F (§2.2) | Accept | `:183` says "fixed representation". The truth is representation-specific: EEG truly faithful share 4.6–6.6% whitened, 0.0% unwhitened | rewrite |
| R10 | Correctness tests: zero noise, exact references, zero and full acceptance, ties, count identities, sub-block weighting (§2.3) | Accept | Partly covered by the synthetic tests (count consistency; accuracy within 0.10–0.15 in an easy regime). Sub-block weights are uniform by construction (`analyze_things_confirmatory.py:28–34`). Only subjects are resampled | `tests/test_estimand.py` |
| R11 | Validation too thin (20 or 10 draws, median bias only); use the ADEMP framework (§3) | Accept | E1, E2 and E7 | Stage 1 |
| R12 | Three layers: references only, end to end, study-level coverage (§3.1) | Accept | Layer A first (cheap, with deep references). Layer C only for claims still standing after A and B | Stage 1 |
| R13 | Three truth families × noise ×{0.1, 0.5, 1, 2} plus zero noise; dependent noise; a development pilot (§3.2) | Accept, changed | Families built from the data: T1, the published truth; T2, low-rank; T3, matched to the measured 8th-to-20th-neighbor gaps (0.16 EEG, 0.30 mouse). cv-MDS included (E7). The dependent noise is a shared session or stratum component | Stage 1 |
| R14 | Stress test: truly faithful neighborhoods judged against shallow, noisy references (§3.2) | Accept, top priority | E2: the mouse curves are flat, and the bias is already +0.04 at R = 16 | Stage 1 |
| R15 | 100 repetitions per cell; 1,000 virtual studies; Monte Carlo errors; timing; caching; skip bootstrap recurrence (§3.3) | Accept | A t-SNE map takes under 1 s here; Layer C costs about 5 h of wall time per scenario and map type | Stage 1 |
| R16 | Full performance measures; Pr(L > Θ) and Pr(U < Θ); the guard; no tuning on validation truths (§3.4) | Accept | E1 argues for revising the estimator. Any revision is chosen on development simulations only and judged on fresh ones | Stage 1 |
| R17 | Finite-depth predictive check; 2× and 4× deeper references in simulation (§3.5) | Accept | Computed (E3). It cannot detect saturation, so it needs a companion flag | Stage 1 |
| R18 | Explicit statuses (estimable, insufficient acceptance, unresolved) from observable criteria (§3.5) | Accept | Criteria: too few acceptances; saturation; finite-depth error; the extrapolants diverging. Thresholds fixed on development data | protocol §4 |
| R19 | Trial manifest; the EEG files give only session membership (§4.1) | Accept | True (S1.1). The authors' public code may show whether the within-session order is recoverable; no download is needed | Stage 2 |
| R20 | Rebuild the EEG before the dataset's noise normalization (MVNN); refit it on F (§4.2) | Decline, disclose | **Cost:** it needs the raw recordings (a new, large download) and the authors' pipeline. **Expected effect:** negligible. The per-session MVNN is a residual covariance over thousands of trials, applied identically to every block (blocks are session-balanced), and it cannot move image-specific signal between blocks. **What the review assumes:** U and H are not in different coordinate systems. **The real gap:** the paper does not say that MVNN used all trials (`:265`); our pilot report does → disclose | rewrite |
| R21 | EEG blocks from separate sessions at R = 8, 16 (12 rotations), against a budget-matched random split (§4.2) | Accept | Cheap with the recovered session labels. Descriptive; does not cover R = 32 | Stage 2 |
| R22 | Mouse: re-split the continuous timeline, with purging and guard intervals (§4.3) | Partly | **Overlap:** the pre-specified 0–0.25 s window has no overlap between presentations and gives the same result (Table S6: 100.0%, L 99.8%). About 0.07 adjacent U–H presentation pairs are expected per image pair. **Feasibility:** blocks from separate epochs are infeasible at R = 8 (in the first experiment, some images have 7 presentations in an epoch). **Direction:** shared noise biases toward agreement (N4). → report the overlap accounting, without re-splitting | Stage 2 |
| R23 | Gate A by dataset; never call finite-reference disagreement "latent false comfort" (§5) | Accept | The rules are written into the protocol before any new result | protocol §5 |
| R24 | Native operating points plus matched acceptance at 5–50% (§6.1) | Accept, scoped | Matched curves only for checks with a continuous score: plug-in (by naive reach), trustworthiness, the scDEED score, the cross-validated reach and the bootstrap support. Seed and split-half recurrence have 2–3 attainable points. Needs a rerun that keeps per-point scores | Stage 2 |
| R25 | Joint acceptance and yield within subject before medians; expose denominators (§6.1) | Done, extend | Already computed within subject (`analyze_checks.py:49–61`). The tables will add accepted counts and contributing subjects | Stage 2 |
| R26 | Random comparator matched per map; yield and enrichment (§6.2) | Accept, high value | E4. Computable from the stored per-map counts | Stage 2 |
| R27 | Audit scDEED; don't pick a null because it passes more points (§6.3) | Accept | **Undocumented:** one permutation per map, null pooled over that map's n scores (`reach/checks.py:66–76`). **Wrong rationale:** `:397` justifies the null by the share it passes; the right reason is fidelity to scDEED's gene permutation. **Sensitivity:** add the coordinate-permutation null | Stage 2 |
| R28 | Audit cv-MDS: negative distances, zero ties, ranking convention, MDS fit (§6.3) | Accept | Naive labels rank the signed d̂² (`run_things_confirmatory.py:110`); the map uses √max(0, d̂²) (line 67). Report the frequencies and the stress | Stage 2 |
| R29 | Report the total repetitions beside R; not a cost-matched improvement (§6.3) | Accept | EEG uses 8 + R + 40 of 80 repetitions; mouse 10 + R + 24 of 50. "Which the estimator makes practical" (`:146`) must state the cost | rewrite |
| R30 | Replace the impossibility argument with a limited convergence statement (§7.1) | Accept, changed | **(a) Overstated argument:** the same-trial argument (`:145`, `:652–658`, `:698`, S1.3) is overstated, and our own results contradict it: split-half on top of plug-in at 32 trials, 20–43% against 48–64% (`:570`), and the enrichment in E4. → state only what sufficiency implies. **(b) Misplaced lemma:** the proposed lemma is about the reference converging, not about same-trial checks → it belongs in S1.2, replacing a premise that the pilot contradicts (N7) | rewrite |
| R31 | Section-by-section changes (§7.2) | Accept | **Exception:** §3.6 is already flagged "not … a tested claim" (`:604`), but its title "Why the check fails" and practice item 5 of §4.1 use it causally → retitle and reword. **The rest** is accepted row by row | rewrite |
| R32 | Keep the five structured-abstract headings; write Results last (§7.2) | Accept | Required by the journal | rewrite |
| R33 | The Figure 2a caption is wrong for EEG cv-MDS (§7.3) | Accept | Called faithful 2.0% vs truly faithful 2.5–3.5% (Table 2); the caption (`:457`) says "far more" | rewrite |
| R34 | "Four checks passed as many as plug-in or more" is wrong (§7.3) | Accept | `:555`. Cross-validated distances pass 7.9% vs 18.8% (EEG) and 2.3% vs 16.1% (mouse); scDEED passes 8.8% vs 16.1% (mouse) | rewrite |
| R35 | False comfort comes from joint labels, not marginal shares (§7.3) | Accept | "Hence" (`:466`): the marginal shares only imply false comfort ≥ 1 − 6.6/18.8 ≈ 65% | rewrite |
| R36 | A positive cross-validated point estimate is not established distinguishability (§7.3) | Accept | "Reliably distinguishable" (`:615`) | rewrite |
| R37 | The plasmode truth is exact by construction; its representativeness is a separate issue (§7.3) | Accept | `:681` says the opposite | rewrite |
| R38 | Pilot certificates: report how often they were issued; drop the multiplicity claim unless tested (§7.3) | Accept | E8 | rewrite |
| R39 | Five main figures, including a representative map (§8) | Partly | **Adopted:** calibration, reference-depth and acceptance-versus-error figures. **Declined as a main figure:** the representative map is illustrative and reads like per-point certification. It goes to the supplement at most | rewrite |
| R40 | Continuous overlap \|E_i ∩ N_k′(reference, i)\|/k, as a finite-reference measure (§8) | Accept | Cheap. Most informative where the binary label saturates (mouse V1) | Stage 2 |
| R41 | Move the robustness figure and subject-level sensitivity to the supplement (§8) | Accept | — | rewrite |
| R42 | "A smaller, calibrated effect beats a larger, uncalibrated headline" (§8) | Accept | — | — |
| R43 | Cite Morris, White & Crowther (2019), doi 10.1002/sim.8086 | Accept | To be verified against the publisher's record when added | rewrite |

## 4. Comments added by the audit

- **N1 (E1).** Report the plasmode at every budget and for both extrapolants. At 16 and 32 trials, stop presenting the 1/R_h extrapolant as the upper end of an honest range.
- **N2 (E2).** The mouse false comfort, and the EEG unwhitened cell, are saturated.
  - Unless the validation shows that they can be resolved, report them as disagreement with a finite reference, together with the exploratory evidence on signal-to-noise.
  - "C1 held in 12 of 12 mouse variations" then is not a robustness result.
- **N3 (E5).** "Conservative bound" (`:236–241`) presumes a coverage that was never shown. Call it the pre-specified lower bound until Layer C shows its coverage.
- **N4. Shared noise pushes the estimate down, not up.**
  - Noise shared by blocks U and H (a session, a drift, overlapping response windows) becomes part of the reference's limit.
  - That makes analyst and reference agree more, which is expected to bias false comfort downward. Such dependence therefore cannot manufacture a high false comfort.
  - The dependent-noise cell of the validation tests this, and Limitations will say so. It lowers the priority of the review's WP2.
- **N5 (E7).** Validate the cv-MDS estimator that C1-RSA rests on, and state that the bias guard applied to it used the t-SNE bias.
- **N6.** The dataset's authors estimated MVNN from all trials of a session. Our pilot report says so (§2); the manuscript does not.
- **N7. S1.2's premise is contradicted.**
  - S1.2 rests on the trace terms being "nearly equal across images".
  - The pilot measured the spread of per-image noise variance at 1.63 times its homoscedastic expectation (1.21–3.45 across subjects).
  - The convergence statement of R30(b) does not need this premise.
- **N8. Highlights and cover letter.** They repeat claims under review:
  - the highlight "~100% (mouse V1)";
  - the cover letter's "Methodological innovation: an estimator", which sits against our own "not new in kind".
- **N9. Results that bound the scope get their own section.** E1–E3 (where the estimator can and cannot resolve false comfort) narrow the reach of the earlier results. They get a Results section of their own rather than a limitations bullet. This is a judgment call, and a smaller edit is possible.

## 5. What this means for the claims, before any new analysis

This is provisional; the validation decides.

| Claim | Status now | What decides it |
|---|---|---|
| **C1, EEG** (t-SNE, 8 trials) | Point estimate 81% (sqrt). Recovered well at this budget (plasmode bias −0.008, RMSE 0.030). The 2nd-smallest subject estimate is 0.605 | Layer C |
| **C1-RSA, EEG** | Margin thin (2nd-smallest 0.530); estimator never validated | Layer B, then C |
| **C2** | Not supported, and stays so: L is from sqrt, biased by +0.026 at R = 32 | — |
| **C3, EEG** | The mislabel rate is recovered within an RMSE of 0.03 at every budget; robust | — |
| **Mouse V1 (C1, C1-RSA, C3)** | Rest on saturated references; latent claims unsupported as they stand | Stress test; finite-reference statement |
| **Checks, EEG at 8 trials** | "Leaves most" for trustworthiness, scDEED, cross-validated distances and seed recurrence survives the stricter order statistic. Enrichment over random is 2.0–4.2 | Layer C, if the classification is kept as an inferential claim |
| **Checks, mouse** | Uninformative: no check can be distinguished from a random selection | — |
| **Ranges at 16 and 32 trials** | Their upper ends are overstated (E1) | Revise the reporting |

The validation, the revision analyses and the rules for this decision are in `docs/revision_protocol.md`.
