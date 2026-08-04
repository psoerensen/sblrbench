# Reorganization inventory

Study 01 fine-mapping is migrated to the shared framework. Its authoritative
source is `spec.R`, `analysis.R`, `report.qmd`, and `locus-design.R`.

## Final cleanup status

Studies 01--05 now use the common runner and clean source contracts. The
`00_contract_smoke` pseudo-study, report/five-replicate compatibility loaders,
migration refresh drivers, and per-development-study launch wrappers were
removed after caller and test audits. Root targets dispatch is retained only
for the Study 07 development graph. Study 06 uses the shared runner with
qualification failed and final fitting remains blocked. Study 07's authoritative internal identifier is
`07_mt_validation`.

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
- `benchmark-execution.R`: `run_benchmark()` is authoritative for Studies
  01--06; Study 06 adds explicit validation, qualification, and final modes.
- `benchmark-extraction.R`: scalar, prediction, parameter, fine-mapping,
  operator, and annotation extraction are shared; MT remains deferred.
- task metric files: stable metric families include prediction, parameter
  estimation, fine-mapping, operator, and annotation metrics.
- `metrics-operator.R`: reusable Study 05 matrix and action-error metrics;
  block/eigen construction policies remain local to Study 05.

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
Affected Study 03 and LD-operator diagnostics now use shared semantic checkpoint
identities; old source-hashed local caches are retired and not reusable.

Historical LD-operator targets, launchers, promotion code, duplicated helpers,
and v2 orchestration were removed during migration. The remaining operator and
SBayesR helper implementations were consolidated into the single Study 05
`operator-design.R`, which paused Study 07 now sources directly.
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
