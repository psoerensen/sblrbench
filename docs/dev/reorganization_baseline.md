# Reorganization behavioral baseline

This document freezes the behavioral and scientific reference for the
`sblrbench` structural reorganization. It records the clean committed state
before shared-helper consolidation; it is not a new scientific analysis.

## Provenance

- `sblrbench` commit: `8507023f2d84f98ac80c30d22067213801de840e`
  (`Study 6 updated`).
- `sblrbench` version: `0.0.0.9000`.
- Validated isolated `sblr`: version `0.2.0`, SHA
  `02e8c74baa906e83c4a08d42a9cc6339b4e81072`.
- qgdata source: `6cca5819e711d326cfb2614d7e9d9f34942612cd`.
- R: `R version 4.4.1 (2024-06-14 ucrt)`.
- Clean-baseline package tests: 481 passed, 0 failed, 0 warned, 0 skipped.
- Clean-baseline Quarto render: all 12 configured pages rendered successfully.

The isolated library at `results/local/current_benchmark_refresh/rlib` was
prepended before inspecting `sblr`. The older normal-user-library installation
was not used.

## Scientific studies

- [Study 01](../../studies/01_finemapping/report.qmd): completed
  separated-locus fine-mapping benchmark.
- [Study 02](../../studies/02_prediction/report.qmd): completed held-out
  single-trait prediction benchmark.
- [Study 03](../../studies/03_parameter_estimation/report.qmd):
  completed single-trait parameter-estimation benchmark.
- [Study 04](../../studies/04_convergence/report.qmd): completed
  convergence selection and five-replicate validation.
- [Study 05](../../studies/05_annotation_models/annotation-convergence.qmd):
  completed prespecified annotation convergence stop; no scientific benchmark
  followed the stop.
- [Study 06](../../studies/06_ld_operator/report.qmd): completed
  retained-low-rank/block LD-operator validation.
- [Supplemental Study 06](../../studies/06_ld_operator/sbayesr_ld_robustness/report.qmd):
  completed SBayesR exact, hard-sparse, and block-eigen robustness evidence.
- Study 07 remains reserved for paused future MTBLR validation. Its scaffold is
  not current scientific evidence and is not resumed by this reorganization.

Study 06 therefore has two preserved parts: its existing retained-low-rank
operator evidence and the supplemental SBayesR LD-operator robustness evidence.

## Frozen reference capsules

Eight authoritative capsule directories were present and all 168 recorded
checksums passed. The machine-readable inventory is
[reorganization_baseline_capsules.csv](reorganization_baseline_capsules.csv).
No capsule was changed or regenerated for this baseline.

## Quarto website

The configured pages are `index.qmd`, `framework.qmd`, `metrics.qmd`,
`reproducibility.qmd`, `studies/index.qmd`, and the seven report pages listed
above (including both Study 06 pages). All scientific reports read frozen
`results/reference/` capsules; none reads `results/local/`, loads a fit
checkpoint, or invokes a sampler.

## Targets and launch entry points

The root `_targets.R` dispatches through `SBLR_BENCH_STUDY` to
`studies/<study>/targets.R`. Current target files exist for the contract smoke
study and Studies 01--07; Study 04 also has `validation_targets.R`. Study launch
scripts remain under `scripts/`, including the current refresh, Studies 05 and
06, the Study 06 v2 validation, and the paused Study 07 scaffold.

Targets remains the execution engine. This reorganization does not run a target
pipeline, rewrite `_targets.R`, or introduce `run_benchmark()`.
