# Reorganization dependency map

## Targets entry points

`_targets.R` is now limited to the preserved Study 06 and Study 07 development
graphs. Completed Studies 01--05 use the common CLI and runner. The contract
smoke pseudo-study, completed-study transition dispatch, and phased launchers
were removed. No target graph was run during cleanup.

## Study-to-study dependencies

- Studies 01--06 use `R/benchmark-data.R` for the pinned qgdata panel. The
  historical cross-study dependency on Study 01 data setup has been removed.
- Study 02 and mechanical callers in Studies 05--07 now use
  `R/benchmark-data.R` for splits, scaling, summary statistics, and training LD;
  no study sources Study 02 internals.
- Study 04 now derives its matched grid from the Study 03 spec and calls shared
  data, simulation, method, extraction, checkpoint, convergence, and execution
  helpers. It has no per-study targets graph or cross-study source dependency.
- Study 06 annotation work uses shared data helpers but remains in development
  and scientifically unmigrated.
- Study 05 uses the shared specification, execution, reporting, checkpoint,
  and operator-metric layers. Its block/eigen and integrated SBayesR audit
  policies remain in `operator-design.R`.
- Paused Study 07 sources Study 05 operator helpers and shared data helpers. Its
  scientific workflow remains paused and unmigrated.
- Before Phase 2, Studies 01 and 03--07 sourced Study 02 promotion solely to
  obtain canonical text-aware MD5 behavior. Active Studies 01 and 03--06 now
  load shared provenance and capsule mechanics directly. Study 02 and paused
  Study 07 now need no compatibility alias.

Study 02's migration supplied the evidence for prediction-only data,
simulation, method, extraction, metric, and execution interfaces. Other task
types remain deferred.

## Duplicated helper families

- Reporting labels, factors, formatting, themes, replicate summaries, and
  capsule-script display live in `R/benchmark-reporting.R`; all reports now
  source that authoritative implementation directly.
- Canonical line-ending-independent MD5 was implemented in Study 02 promotion
  and reused across other studies.
- Atomic uncompressed RDS save logic was repeated in historical Study 06,
  Study 07, and developer diagnostics.
- Input hashes and strict stale-checkpoint refusal occurred in historical Study 06 and the
  supplemental diagnostics.
- The scalar convergence kernel (rank-normalized R-hat, bulk/tail ESS, mean
  MCSE, posterior SD, relative MCSE, and thresholds) was near-identical in
  Studies 05--07.
- Each study has distinct trace extraction around a common three-dimensional
  trace-array shape and retained-chain window mechanics.
- Required-file inventories, canonical checksums, staging directories,
  selected copies, semantic validation, and atomic final rename recur across
  promotion files.

The Study 03 SBayesR comparison and integrated Study 05 SBayesR evidence
formerly hashed Study 03 implementation files. They now use
`sblrbench-semantic-v2` identities containing scientific inputs, controls,
seeds, data/operator hashes, and package provenance. Old local caches are
detected as retired rather than reused or translated.

## Capsule promotion dependency on Study 02

The historical cross-study dependency was mechanical, not scientific:
promotion files needed `.study02_canonical_md5()`. Shared canonical hashing,
capsule inventory, and capsule validation removed it. Study-specific required
outputs, manifests, README text, semantic checks, and destinations remain in
the two development studies until their migrations.

## Reports and package functions

The completed reports and Study 06 status report source
`R/benchmark-reporting.R` directly. Reports otherwise use base R table
operations, `jsonlite`, `ggplot2`, `knitr`, and frozen capsule files.
The sole Study 05 report reads both operator evidence components from one
frozen capsule.

No current report reads `results/local`, `readRDS()` output, a checkpoint, or a
native fit object. No report calls `sblr`, constructs LD, simulates data, or
invokes a sampler. The textual word “simulated” in reports describes frozen
designs and is not executable dependency evidence.

## Concrete problems left for migration phases

- Example-data acquisition remains owned by Study 01; Study 02 mechanics are
  shared, but later data interfaces await their own migrations.
- Study 04 is migrated. Its spec retains only the matched grid, diagnostic
  quantities, candidates, thresholds, chain requirements, and recommendation
  rules; reusable mechanics live under `R/`.
- Study 06 and Study 07 retain development target graphs but no dedicated
  launch wrappers. Their future shared-runner interfaces remain unresolved.
- Trace extraction schemas differ and remain duplicated around study-specific
  quantities.
- Promotion functions still duplicate README/manifest composition and exact
  source selection; only safe mechanics are shared.
- Method dispatch and checkpoint naming remain study-specific.
- Study 06 annotation logic cannot be generalized safely before its migration.
