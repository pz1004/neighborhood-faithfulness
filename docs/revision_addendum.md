# Addendum 1 to the revision protocol: what the development pilot decided

**Status.** Dated 2026-10-04 and committed before any validation seed is used, as `docs/revision_protocol.md` (Section 8) requires. The validation run refuses to start until this file and `results/revision/dev/decisions.json` are committed (`run_validation.py`).

**Sources.**
- The pilot's rows are in `results/revision/dev/`, and its printout is `analysis.txt` (from `analyze_validation.py --phase dev`).
- The pilot used the development seeds (from 6,000,000), 2 repetitions per base subject (20 per cell), and both datasets.
- **Development data are never reused for validation.** Every number below is development evidence. The validation measures performance afresh.

## 1. Decisions

Each decision is listed with the protocol rule that made it.

**r (T2's rank).** 5, unchanged.

**α (T3).** α was chosen on a grid of 0.05 steps to match the median relative gap between the 8th and 20th true neighbor (`results/revision/dev/alpha.json`).

| Dataset | Target (cross-validated) | T2 (α = 0) | T1 (α = 1) | Chosen α | Gap at chosen α |
|---|---|---|---|---|---|
| EEG | 0.159 | 0.295 | 0.115 | 0.45 | 0.163 |
| Mouse V1 | 0.298 | 0.268 | 0.101 | 0 | 0.268 |

- In mouse V1 no α reaches the target, so the nearest value is used, as the protocol says, and recorded. **T3 is T2 there.** The mouse T3 cells are therefore replicates of T2 with their own seeds.
- In both datasets, T1 (the published plasmode's truth) has smaller relative gaps than the cross-validated estimate from the real data.

**The noise grid (λ).** λ ∈ {0.1, 0.25, 0.35, 0.5, 1, 2}.
- **The protocol's test.** With {0.1, 0.5, 1, 2}, the exact false comfort did not span the protocol's ranges:
  - in EEG, no cell lay near 0.25 (0.15–0.35);
  - in mouse V1, no cell lay between 0.25 and 0.5, or near 0.5.
- **The change.** A development extension with λ = 0.25 and 0.35 (`run_extension.log`) closed both gaps.
- **Added, not replaced.** The protocol says to *replace* λ values. They were added instead, so that the regimes 0.1 and 2 already cover, the lowest and highest false comfort, stay in the grid. This is a deviation in form, recorded here.
- **Size of the validation.** There are now 26 cells per dataset.

**The estimator.** The published pair is kept.
- **The selection rule (protocol Section 3.5)** asks for |median bias| ≤ 0.03 in every unsaturated development cell (here, cells whose default status is estimable or extrapolation unresolved).
- **No candidate qualified:**

  | Candidate | Pooled RMSE | Largest \|median bias\| |
  |---|---|---|
  | (a) R_h^−1/2 alone | 0.175 | 0.747 |
  | (b) more reference sizes | 0.179 | 0.769 |
  | (c) cross-validated references | 0.202 | 0.813 |

  The largest biases are in the stress cells.
- **A narrower reading gives the same answer.** Restricted to the screening cells whose default status is estimable, all three still fail. Examples:
  - EEG cv-MDS T1: −0.037;
  - EEG T3 at λ = 0.35: −0.094;
  - mouse V1 at λ = 0.1: −0.032 and −0.033.
- **Consequence (Section 3.5).** Latent estimates are model-dependent by default. How Section 5 still applies cell by cell is fixed in Section 2 below.

**The status thresholds.** The defaults stand: A_min 20, c_sat 0.95, f_sat 0.02, c_fd 0.03, c_div 0.10.
- **No combination met the tuning rule.** The rule tunes the thresholds so that every estimable development cell meets |median bias| ≤ 0.03 and RMSE ≤ 0.08. The grid search (243 combinations) found none, for two reasons:
  - **The dependent-noise cells cannot be flagged by any observable.** Shared noise looks like signal. They are estimable under every threshold that keeps the real data's regime estimable.
  - **The low-noise mouse cells miss by 0.002–0.003.**
- **Not specified by the protocol.** The protocol does not say what happens when no threshold set qualifies. The pre-specified defaults are kept, and this addendum records that choice.

## 2. How Gate A is applied

The reporting rules of protocol Section 5 leave the following points open. They are fixed here, before the validation, and stored in `decisions.json` (`gate_a`).

1. **Whose accuracy.** With the published pair kept, accuracy is assessed on the *smaller* of the two extrapolants, on which the bound L rests. The published range [smaller, larger] is reported with its *bracket coverage*: the share of repetitions whose exact value lies in it. This is a new performance measure of the published reporting.
2. **Matching cells.** These are validation cells of the same dataset and the same map type, with the same status, whose median d(R_h,max) and fall are both within 0.05 of the real cell's.
   - The dependent-noise (D) and stress (S) cells are diagnostics and never match.
   - **D's levels are far from the data.** D tests the direction of the bias under shared noise at ρ = 0.1 and 0.3. The pilot's session × image interaction in EEG (F = 1.14 with about 18 repetitions per session and image, `docs/gate3_pilot_things.md`) puts the shared share near 0.01.
   - **S is flagged anyway.** S shows the failure under noisy references, and the statuses flag it.
3. **Worst case.** Every matching cell must meet the accuracy criterion; the worst matching cell decides.
4. **Two readings.** Two classifications are reported, and the author chooses at Gate A:
   - **strict:** Section 3.5 applies. No candidate qualified, so every latent estimate is model-dependent and carries no threshold claim.
   - **local:** Section 5 applies cell by cell. A real cell may carry a threshold claim if its matching cells meet the accuracy criterion and Layer C calibrates its bound.
5. **Layer C** runs under the local reading only, for real cells whose matching cells are accurate.
   - **Projected cost.** About 4.4 h of wall time per scenario for EEG t-SNE without the guard, and about 6.3 h with it.
   - **The guard.** Following Section 3.4, it is computed in the first 200 studies only, and the rate at which it triggers is reported.
   - **Size.** 1,000 studies per scenario and map, since 4.4 h is under the 8 h limit.

## 3. Timing check (Section 9)

Median seconds per repetition, with 35 workers on 36 threads (`time.jsonl`):

| | EEG | Mouse V1 |
|---|---|---|
| One budget, t-SNE, with Layer A | 70 | 27 |
| Every budget (3 in EEG, 2 in mouse) | 175 | 53 |
| cv-MDS | 21 | 6 |

The first 20 runs had a median of 53 s. The validation is projected at about 90 min for EEG and 35 min for mouse V1.

## 4. What the development pilot showed

This is development evidence only, not a result.

**EEG, the regime of the real data** (λ = 1, 8 trials):
- The smaller extrapolant's median bias was +0.003 (T1), −0.017 (T2) and −0.021 (T3), against an exact false comfort of about 0.8.
- 1/R_h overstated by 0.04–0.05.

**Moderate false comfort** (λ 0.25–0.5, exact 0.13–0.50): R_h^−1/2 underestimated by 0.04–0.14, while 1/R_h was within −0.005 to +0.035.

**Neither extrapolant is right across regimes.** The exact value usually lay between them:

| | Smaller | Larger |
|---|---|---|
| 16 trials | about 0 | +0.09 to +0.12 |
| 32 trials | +0.05 to +0.08 | +0.18 to +0.23 |

**Stress cells.**
- The estimator overstated by +0.18 and +0.21 in EEG, and by +0.55 and +0.75 in mouse V1.
- Every stress cell was flagged "extrapolation unresolved".
- None of them was saturated: the deepest reference's disagreement stayed at or below 0.934, and it fell by at least 0.058.

**Dependent noise.** At ρ = 0.3 the estimate fell below the exact value (EEG −0.115; mouse −0.204), as audit item N4 expects.

**Saturated cells.**
- In mouse V1, every cell with λ ≥ 1 had an exact false comfort of at least 0.977 at 8 trials.
- Across both datasets, every saturated development cell had an exact false comfort of at least 0.915 (mouse T1 at 16 trials).
- In the pilot, no simulated mechanism produced a saturated curve with a low exact value.

## 5. Changes to the code since its commit (`817f138`)

- `run_validation.py`:
  - development-only options `--lams` and `--append`, for the extension;
  - rows written as gzipped JSON lines (the 24 MB files compress to 1 MB).
- `analyze_validation.py`:
  - when no threshold set qualifies, the defaults are kept;
  - the narrower reading is reported as a sensitivity;
  - the Gate A clarifications are written into the decisions file;
  - the bracket coverage is added.
- `tests/test_validation.py` still passes (`results/revision/test_validation.log`).
