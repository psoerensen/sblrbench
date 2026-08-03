# Study 02 prediction migration

Study 02 is the first study migrated to the common execution framework. The
scientific reference remains `results/reference/02_prediction/current`; no
capsule was regenerated or promoted.

## Responsibility map

| Previous responsibility | Authoritative implementation |
|---|---|
| `config.R` scientific constants | `studies/02_prediction/spec.R` ordinary list |
| `pilot.R` data/simulation/fitting | `R/benchmark-data.R`, `R/benchmark-simulation.R`, `R/benchmark-methods.R` |
| `targets.R` orchestration | prediction-only `run_benchmark()` in `R/benchmark-execution.R` |
| inline extraction and metrics | `R/benchmark-extraction.R`, `R/metrics-prediction.R` |
| `promotion.R` capsule checks | shared validation plus frozen capsule; no promotion in Phase 3 |
| dedicated launcher | `scripts/run_benchmark.R` and `.ps1` |
| worked example | `inst/templates/prediction-analysis.R` |
| `prediction.qmd` | `report.qmd`, still capsule-only |

The old config, pilot, targets, promotion, README, dedicated launcher, and
worked-example files were deleted after source/link searches and caller
updates. Studies 05--07 now call shared data/provenance functions directly;
their scientific specifications and orchestration were not migrated.

## Scientific preservation

Regression checks cover the 40-coordinate benchmark grid, scenario and method
order, five replicates, architecture/simulation/fit/chain seeds, split and
marker ordering, exact controls and priors, metric definitions, and output
schemas. A deterministic two-scenario simulation fixture reproduces the
pre-migration object hashes. The frozen Study 02 capsule is checksum-validated
and recursively compared with its pre-migration file inventory; all other
reference capsules receive the same byte-identity audit.

The report reads only its frozen capsule. It does not load checkpoints, invoke
targets, prepare data, simulate, construct LD, call `run_benchmark()`, or call
an `sblr` sampler.

## Checkpoint decision

The historical Study 02 fits are targets-store internals rather than a stable
checkpoint interface. They remain untouched in `_targets`. The new framework
starts a clean schema whose identity is based on scientific data and marker
order, split, simulated truth, method controls, seeds, package SHA, and data
provenance—not source paths. No legacy loader or fit-object rewrite is kept.
Validation-only execution exits before data preparation or fit dispatch and
makes zero sampler calls.

## Validation and deferred work

Focused tests, parsing, package loading, CLI validation-only execution, capsule
hash comparison, and direct report rendering are recorded in the Phase 3 handoff.
The framework intentionally supports only prediction. Study 03 must define its
parameter-estimation specification, estimands, extraction, and metrics during
its own migration; no speculative interface was added here.
