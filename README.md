# Certified neighborhood reach

Simultaneous, finite-sample confidence intervals on how far each displayed neighborhood of a
2-D display (t-SNE, UMAP, PCA) reaches into the true neighbor ranking, for data measured
with replicates. From the intervals comes a faithful / distorted / unresolved map that is
valid at every tolerance.

**Status: Gate 1 (prior-art audit) and Gate 2 (development validity and power) in
progress.** The protocol is [docs/research_plan.md](docs/research_plan.md), version 1 of
2026-09-28, written before any Gate 2 measurement.

## Provenance

This project succeeds `setwise-certificates-neighborhood-fidelity` (dropped on 2026-09-23;
evidence base frozen at commit `9bc9f03`). Section 2 of the plan reviews the proposal that
started it. `pilot/2026-09-28/` archives the scratch pilot the proposal rested on, unchanged.
The pilot imports the predecessor repository by absolute path and is not part of the
reproduction.

## Layout

| path | contents |
|---|---|
| `docs/research_plan.md` | the protocol: claim, review of the proposal, engines, gate criteria, seeds |
| `pilot/2026-09-28/` | the scratch pilot and the tail check behind review finding 1 |
