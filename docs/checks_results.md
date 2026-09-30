# Common checks against the truth: results

2026-10-01. Protocol: `docs/checks_protocol.md`, pre-registered with its code and synthetic tests
at commit `766e4b7`.

**How it was run.**
- The code ran as pre-registered.
- Every data file matched its SHA-256 digest before loading.
- Run times on 36 workers (18 cores): 3,162 s for EEG and 1,258 s for mouse V1.
- The wrapper command reported exit code 1 after both runs had finished. The cause was its
  final `tail -3` over two files, which is invalid syntax. Every expected row was written:
  2,700 (EEG), 1,800 (mouse) and 1,200 semi-synthetic rows each.
- No deviation from the protocol.

**Reproduction.** Every re-created map gave the same plug-in counts as the committed
confirmatory rows: 2,700 of 2,700 rows in EEG and 1,800 of 1,800 in mouse V1.

The printout is `results/checks/checks_analysis.txt`, and the numbers below come from
`results/checks/checks_summary.json`. Every estimate is a median over 10 subjects. Pairs are
given as (linear in R_h^(−1/2)) / (linear in 1/R_h). L and U are the registered conservative
bounds.

## 1. Verdict at 8 trials (protocol Section 6)

| check | EEG: passes | EEG: false comfort | EEG: L, U | EEG: class | mouse: passes | mouse: false comfort | mouse: L, U | mouse: class |
|---|---|---|---|---|---|---|---|---|
| N, plug-in (reference) | 18.8% | 0.811 / 0.865 | 0.642, 0.915 | (C1) | 16.1% | 1.000 / 1.000 | 0.997, 1.000 | (C1) |
| TW, trustworthiness | 19.0% | 0.792 / 0.853 | 0.643, 0.909 | leaves most | 16.3% | 1.000 / 1.000 | 0.998, 1.000 | leaves most |
| SCD, scDEED | 32.5% | 0.863 / 0.902 | 0.759, 0.938 | leaves most | 8.8% | 1.000 / 1.000 | 0.997, 1.000 | leaves most |
| CV, cross-validated | 7.9% | 0.728 / 0.805 | 0.595, 0.856 | leaves most | 2.3% | 1.000 / 1.000 | 0.979, 1.004 | leaves most |
| SEED, seed recurrence | 34.6% | 0.870 / 0.907 | 0.771, 0.943 | leaves most | 27.2% | 1.000 / 1.000 | 0.998, 1.000 | leaves most |
| BOOT, bootstrap recurrence | 0.0% | defined in 3 subjects | — | not evaluable | 0.0% | no passes | — | not evaluable |
| SH, split-half recurrence | 1.2% | 0.467 / 0.621 | 0.383, 0.722 | undetermined | 0.0% | defined in 4 subjects | — | not evaluable |

**What the manuscript may say under the rules:**
- **None of the common checks removed most false comfort in either dataset.**
- **Four checks left most false comfort in both datasets:** trustworthiness, scDEED,
  cross-validated distances and seed recurrence.
- The sentence "all six checks left most false comfort in both datasets" is **not** permitted.
  - Bootstrap recurrence passed essentially nothing at 8 trials in either dataset, so it is
    not evaluable.
  - Split-half recurrence is undetermined in EEG and not evaluable in mouse V1.

## 2. Estimator check (protocol Section 5)

The median bias at 8 trials is given as smaller estimate / larger estimate, minus the exact
value:

| check | EEG | mouse V1 |
|---|---|---|
| N | −0.002 / +0.042 | +0.003 / +0.004 |
| TW | −0.003 / +0.045 | +0.003 / +0.005 |
| SCD | −0.002 / +0.026 | +0.004 / +0.006 |
| CV | −0.007 / +0.049 | +0.002 / +0.008 |
| SEED | +0.002 / +0.028 | +0.003 / +0.005 |
| SH | −0.014 / +0.016 | +0.000 / +0.000 (4 subjects) |
| BOOT | 1 subject only | no subject |

- **The guard never applied.** No conservative bias exceeded +0.05, and no larger estimate fell
  more than 0.05 below the exact value.
- **BOOT is not validated.** In the simulations it passed too few neighborhoods to check the
  estimator: one EEG subject and no mouse subject.
- **SH is validated in EEG.** Its exact false comfort per subject ranged from 0.475 to 1.000, and
  both estimates followed it.

## 3. What the checks did (descriptive)

**The four permissive checks pass as much as the plug-in check, or more, and are wrong as often.**
- Trustworthiness was rank-matched to the plug-in check's pass count. It flagged about 27–28% of
  the plug-in check's false-comfort neighborhoods in EEG and 20% in mouse V1. Its own false
  comfort was the same as the plug-in check's (0.79–0.85 in EEG).
- scDEED and seed recurrence passed about a third of EEG neighborhoods, against 19% for the
  plug-in check. Their false comfort was higher: 0.86–0.91.
  - About 20% of all points passed them but failed the plug-in check (passes 0.325 and 0.346;
    passes shared with N 0.123 and 0.126).
  - Stability is not faithfulness: those extra passes are neighborhoods that are distorted in
    the same way across seeds or against a permutation null.
- Cross-validated distances passed fewer neighborhoods (7.9% in EEG, 2.3% in mouse V1). They
  flagged 72–74% of the plug-in check's false comfort in EEG and 93% in mouse V1, but 73–81% of
  what they passed in EEG was still wrong.

**The two strict checks pass almost nothing at 8 trials.**
- Bootstrap recurrence asks that a neighborhood recur in at least 95% of 50 bootstrap maps. It
  passed essentially no neighborhood in either dataset.
- Split-half recurrence passed 1.2% of EEG neighborhoods and essentially none in mouse V1. It
  flagged 97–98% of the plug-in check's false comfort in EEG. Of what it passed there, 47–62% was
  still not faithful; the yield was 0.4–0.6% of all points.

**With more trials in EEG, the strict checks begin to pass mostly faithful neighborhoods.** These
budgets are descriptive, not claims.
- Split-half recurrence at 32 trials passed 14.4% of points, with false comfort 0.32–0.52.
- Split-half recurrence used on top of the plug-in check at 32 trials passed 9.8%, with false
  comfort 0.20–0.43 (U = 0.498).
- Bootstrap recurrence at 32 trials passed 1.5%. Its false comfort was estimated at −0.12 /
  0.15, with U = 0.27; the negative value is an extrapolation artifact.
- Bootstrap recurrence at 16 trials passed 0.15%, with false comfort 0.00 / 0.25 (U = 0.43),
  defined in 7 subjects.
- The estimator was not validated for bootstrap recurrence (Section 2).
- The plug-in check at 32 trials had false comfort 0.48–0.64.

**Mouse V1.** Every check that passed any neighborhood was wrong about essentially all of them
(false comfort ≥ 0.97 at both budgets). The checks differ only in how many they passed.

**Yield in EEG at 8 trials.** Yield is the share of all points passed and faithful:
- the plug-in check: 2.5–3.5%;
- scDEED and seed recurrence: 3.3–4.9%;
- cross-validated distances: 1.6–2.2%;
- split-half recurrence: 0.4–0.6%.

## 4. The expectation recorded before the run (protocol Section 7)

**The prediction:** all six checks would leave most false comfort in both datasets.
- **It held for the four checks that pass many neighborhoods.**
- **It did not hold as stated for the two strict ones.** They passed too few neighborhoods at 8
  trials to be classified, except split-half recurrence in EEG. That check's false comfort,
  0.47–0.62, is lower than the plug-in check's 0.81–0.87 but not clearly below 0.5.

**The part of the argument that survives.** No check computed from the analyst's trials told
faithful from unfaithful neighborhoods at 8 trials. A check could only become stricter, trading
yield for precision:
- bootstrap and split-half recurrence pass only neighborhoods with large margins relative to the
  noise;
- at 8 trials almost none exist;
- where the strict checks do pass neighborhoods, in EEG at 16–32 trials, more of those passes are
  faithful.

**"Stability is not faithfulness" held.** Seed recurrence and scDEED passed about a fifth of all
points that the plug-in check failed. Their false comfort was higher than the plug-in check's, so
those extra passes were mostly unfaithful.

## 5. What the manuscript can say now

**Supported by the pre-registered rules on both datasets:**
- Four common checks left most false comfort at 8 trials: per-point trustworthiness, scDEED's
  trustworthy label, cross-validated distances and seed recurrence.
- None of the six removed most false comfort.

**Reported as observed, with their classification:**
- The two strict resampling checks passed almost nothing at 8 trials. Bootstrap recurrence is
  not evaluable in either dataset.
- Split-half recurrence is undetermined in EEG, with false comfort 0.47–0.62 on 1.2% of points,
  and not evaluable in mouse V1.

**Descriptive only:** at 32 trials in EEG, split-half recurrence on top of the plug-in check
passed 9.8% of points with false comfort 0.20–0.43.

**Limits:**
- the primary analysis only;
- the checks are as specified, for example bootstrap support of 95% of 50 maps;
- the estimator was not validated for bootstrap recurrence, which passed too few neighborhoods
  in the simulations.
