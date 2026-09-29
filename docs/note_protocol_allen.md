# Confirmatory protocol, second dataset: false comfort in 2-D maps of mouse V1 responses

Version 1, 2026-09-29. Plan Amendment 4 points here.

**What was known when this was written.** Only Allen's experiment metadata
(`ophys_experiments.json`) had been downloaded. No recording file had been downloaded,
and no response value had been read. The code that implements this protocol is committed
separately, before the analysis runs.

## 1. Purpose

The THINGS-EEG2 confirmatory analysis (`docs/note_confirmatory_results.md`) supported the
note's headline: most neighborhoods that a plug-in check calls faithful in 2-D maps of
trial-averaged EEG are not faithful in truth. That analysis reused the pilot's data.

This protocol tests the same claims on independent data. The species, the recording
modality, the stimuli and the noise sources are all different:
- mouse primary visual cortex (V1);
- two-photon calcium imaging of single cells;
- 118 natural scenes, each shown 50 times.

## 2. Data and selection rule

**Source.** The Allen Brain Observatory Visual Coding two-photon dataset (de Vries et al.
2020), from Allen's public bucket (`allen-brain-observatory`, `visual-coding-2p/`). It is
used under the Allen Institute terms of use: non-commercial research, with attribution.

**Candidates.** Every experiment satisfying all of:
- `stimulus_name` is `three_session_B`, the session that contains the natural scenes;
- the targeted structure is VISp (V1);
- the donor's Cre driver is `Cux2-CreERT2`;
- the imaging depth is 175 or 275 µm.

A mouse with both depths contributes only its 175 µm experiment. The metadata lists 17
such mice.

**Order and inclusion.** Mice are taken in ascending donor ID, and each experiment file is
downloaded in turn. An experiment is included if both hold:
- its file has at least 100 cells (ROIs) with a ΔF/F trace;
- its natural-scenes stimulus table has exactly 50 presentations of each of the 118 images.

This stops once 10 experiments are included. If fewer than 10 of the 17 qualify, all that
qualify are used, and the shortfall is reported.

**Integrity.** Every file's SHA-256 is recorded in `data/allen_2p/manifest.json` when it
is downloaded and checked again before it is analyzed.

## 3. Responses

- **One trial.** A trial is one presentation of an image. Blank presentations are excluded.
- **Primary response.** Each cell's mean ΔF/F over the 0.5 s after the presentation's onset.
- **Alternative window.** The 0.25 s during the presentation. It replaces the THINGS
  analysis's early window in cell 11.
- **Frames.** Onsets and frame times come from the file's stimulus table and imaging
  timestamps.
- **Carryover.** Presentations follow each other without gaps, and calcium signals are
  slow, so a response partly carries into the next image's window. Image order is random,
  so this adds noise, not image-specific bias.

## 4. Design

Everything follows `docs/note_protocol.md` except the changes below.

- **Time strata.** Each image's 50 presentations are ranked by time and split into an
  early and a late stratum of 25. These replace THINGS-EEG2's four sessions.
- **Blocks within each stratum** (seed 4_500_000 + 1000·subject + split):
  - F: 5 repetitions, 10 in total;
  - U: 8, 16 in total;
  - H: 12, 24 in total.

  Each block is ordered round-robin over the two strata, so every prefix of even length is
  balanced.
- **Budgets.** User R ∈ {8, 16}. There are 50 repetitions, so R = 32 is not available.
- **Reference sub-blocks.** R_h ∈ {8, 12, 24}: 3, 2 and 1 disjoint, stratum-balanced
  sub-blocks of H.
- **Cross-validated halves** of a user block alternate over its round-robin rows, so each
  half holds R/4 repetitions of each stratum.
- **Representation, fitted on F only.** PCA of F's image means to min(100, cells)
  dimensions, then Ledoit–Wolf whitening, then PCA to D = 20. The unwhitened cell uses PCA
  to 20 without whitening.
- **Subjects.** Each included experiment is one subject, from a distinct mouse.
- **Splits and draws.** 10 splits per subject. The semi-synthetic check uses 20 draws per
  subject, with seeds 4_600_000 + 10_000·subject + draw. The bootstrap uses seed
  4_700_000.
- **Everything else unchanged:** the 12 cells, the display seeds, k = 10, ε = 0.2,
  k′ = 2k, the two pipelines, both extrapolants, the conservative bound L, and the
  semi-synthetic estimator check with its bias guard. In cell 11, "early window" means the
  0.25 s window of Section 3.

## 5. Claims and decision rules

**Claims, all in the primary cell unless stated:**
- **C1.** False comfort at 8 trials > 0.5. Supported if L > 0.5.
- **C1-RSA.** The same in the cross-validated MDS cell. Supported if L > 0.5.
- **C3.** Mislabel rate > 0.10 at R = 8 and at R = 16. Supported if L > 0.10 at each.

C2, the claim at 32 trials, cannot be tested with 50 repetitions and is dropped, not
replaced.

**Generality.** C1's criterion must hold in at least 10 of the 12 cells.

**Bias guard.** As in `docs/note_protocol.md`: if the semi-synthetic check shows the
conservative estimator overestimating false comfort at 8 trials by more than 0.05, L is
lowered by that bias.

**What the note may say, given both datasets:**
- **C1 holds here:** the headline is independently confirmed in a second species and
  modality. It covers RSA-style maps too if C1-RSA also holds here.
- **C1 fails here:** the headline is limited to THINGS-EEG2 or to EEG, and the note must say
  so. If it is written, it reports this result as its limit.

## 6. Integrity

As in `docs/note_protocol.md` Section 7:
- the code at the pre-registration commit is the code that runs;
- only errors that are not analysis choices may be fixed, and each fix is logged;
- every other change is a dated amendment.
