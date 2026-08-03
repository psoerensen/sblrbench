# Reorganization inventory

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

## Obsolete candidates

Names such as `scripts/run_benchmark.R` and v1/v2 Study 06 files are not removed
or classified as obsolete solely from their names. They remain inventory items
until callers, reproducibility roles, and migration replacements are verified.
Git history is the archive; no legacy directory is introduced.
