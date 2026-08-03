# Study 03 parameter-estimation migration

Study 03 now uses the common ordinary-R runner. Its complete scientific list is
`studies/03_parameter_estimation/spec.R`, its exact workflow is `analysis.R`,
and `report.qmd` continues to read the unchanged frozen capsule.

The downloadable exact workflow now prints the data and scenario design,
coordinate/seed grid, methods and resolved controls, and the authoritative
estimand formulas before calling the shared runner. It then validates fit and
oracle coverage, summarizes convergence and extracted results, assigns and
saves named recovery/bias/probability/runtime plots, and inventories outputs.
No low-level execution or fit-object parsing was moved back into the study.

## Responsibility mapping

The old targets graph is replaced by parameter task dispatch in
`R/benchmark-execution.R`. Shared data, simulation, four-method controls,
checkpointing, runtime, marker extraction, and provenance are reused from the
Study 02 framework. Draw extraction was added to `benchmark-extraction.R`; the
six realized estimands and recovery/coverage/pairing formulas are authoritative
in `metrics-parameter-estimation.R` and the spec registry.

Obsolete targets, promotion, metric, estimand, launch, worked-example, README,
smoke-test, config, simulation, method, and pilot files were removed. Study 04
now consumes the committed spec and shared simulation/truth functions
mechanically; its scientific grid is not migrated. The Study 03 comparison and
Study 06 supplemental diagnostics use the shared `sblrbench-semantic-v2`
checkpoint identity. Historical source-hashed local caches are explicitly
rejected and were not translated or deleted.

## Scientific preservation and checkpoints

The benchmark remains two architectures, five replicates, four methods, 50
causal markers, h2 0.30, and 40 fits. Seeds, controls, priors, draw-wise
heritability, draw-wise `vbs*pi_trace*m`, realized truths, interval summaries,
and paired comparisons are regression-tested. Historical targets/checkpoints
remain untouched; the new schema is a clean break based on scientific identity,
not source paths. Frozen capsules remain authoritative. No fit or capsule
promotion occurred. See [study03_checkpoint_retirement.md](study03_checkpoint_retirement.md).

All capsule checksum and byte-identity results, validation-only CLI runs,
package tests, parsing, and report renders are recorded in the handoff.
