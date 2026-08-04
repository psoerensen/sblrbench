# Reorganization inventory

Study 01 fine-mapping is migrated to the shared framework. Its authoritative
source is `spec.R`, `analysis.R`, `report.qmd`, and `locus-design.R`.

The responsibility-level inventory is
[reorganization_inventory.csv](reorganization_inventory.csv). It was built
from function definitions, `source()` relationships, targets files, launch
scripts, report code, tests, and capsule promotion code at baseline commit
`8507023f2d84f98ac80c30d22067213801de840e`.

## Architectural finding

The repository already uses ordinary functions, lists, data frames, Quarto,
and targets. No R6, plugin system, universal metric class, registry framework,
workflow language, or inheritance hierarchy is needed. The only committed
`docs/implementation_plan.md` predates this reorganization and describes the
initial package implementation; it is not a Phase 0–2 reorganization plan.
The detailed approved plan supplied for this work is therefore authoritative.

## Phase 2 evidence

Six files contain useful consolidated code now:

- `R/benchmark-reporting.R`: the existing shared report-only implementation.
- `R/benchmark-provenance.R`: Git/package provenance and canonical hashes.
- `R/benchmark-checkpoints.R`: atomic RDS and strict validated reuse.
- `R/benchmark-convergence.R`: trace windows, long extraction, and scalar
  R-hat/ESS/MCSE calculations.
- `R/benchmark-capsules.R`: checksum inventories and bounded staging/copying/
  atomic promotion mechanics.
- `R/benchmark-validation.R`: structural capsule checksum validation.

These files separate repeated mechanics from study-specific required-file
lists, scientific thresholds, checkpoint identities, and semantic validators.

## Deferred proposed files

The following approved names are deliberately not created:

- `benchmark-spec.R`: study specifications remain embedded in functioning
  configs and targets files; no stable common interface is yet proven.
- `benchmark-data.R`: Study 01 and Study 02 internals are shared in practice,
  but a first study migration must define the supported data boundary.
- `benchmark-simulation.R`: architectures and truth objects differ by study.
- `benchmark-methods.R`: dispatch and controls remain scientific choices.
- `benchmark-execution.R`: targets remains the execution engine and no
  `run_benchmark()` abstraction is approved yet.
- `benchmark-extraction.R`: scalar, annotation, operator, and MT trace schemas
  differ; only the common trace-array primitive is consolidated.
- task metric files: existing metric interfaces need a first study migration
  before task-specific placement is stable.
- `metrics-operator.R`: Study 06 operator semantics remain local.

No empty placeholder or speculative abstraction is created for these items.

## Phase 3: Study 02 migration

Study 02 is the first migrated study. Its scientific contract now lives in
`studies/02_prediction/spec.R`; its exact readable entry point is
`studies/02_prediction/analysis.R`; and its frozen-capsule report is
`studies/02_prediction/report.qmd`. Shared data, simulation, method translation,
extraction, prediction metrics, specification, and execution mechanics now live
under `R/benchmark-*.R` and `R/metrics-prediction.R`.

The former `config.R`, `pilot.R`, `targets.R`, `promotion.R`, README, dedicated
launcher, and worked-example script were removed after repository callers were
updated. This is an internal clean break: the historical `_targets` store is
left untouched, while the new runner uses science-identity checkpoints.

## Obsolete candidates

## Phase 4: Study 03 migration

Study 03 parameter estimation now uses the common runner, shared full-sample
data/simulation/method mechanics, shared scalar extraction, and task-specific
parameter metrics. Obsolete targets, promotion, estimand, metric, launch,
worked-example, config, simulation, method, and pilot layers were removed.
Affected Study 03 and Study 06 diagnostics now use shared semantic checkpoint
identities; old source-hashed local caches are retired and not reusable.

Study 06 v1/v2 files are not removed
or classified as obsolete solely from their names. They remain inventory items
until callers, reproducibility roles, and migration replacements are verified.
Git history is the archive; no legacy directory is introduced.

## Phase 5: Study 04 migration

Study 04 convergence now uses an ordinary-list `spec.R`, the shared two-stage
runner, semantic checkpoints, true-trace extraction, convergence diagnostics,
recommendation mechanics, validation summaries, and table-only plotting.
Historical config, target graphs, diagnostics, launchers, promotion scripts,
examples, and compatibility loaders were removed after caller updates. The
selection and validation capsules remain byte-identical and authoritative.
