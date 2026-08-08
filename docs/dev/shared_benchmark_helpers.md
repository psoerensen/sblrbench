# Shared benchmark helpers after Phase 3

## Study 06 annotation-model extension

The runner supports `annotation_models` with explicit `validate_only`,
`qualification`, and `final` modes. Shared spec, method, extraction,
convergence, reporting, provenance, and semantic-checkpoint functions are
reused. `R/metrics-annotation.R` owns stable annotation recovery metrics;
deterministic scientific construction remains in Study 06
`annotation-design.R`. Comparable BED/CSR marker priors require true retained
alpha traces and are reconstructed draw by draw. Missing traces are reported
as unavailable rather than replaced by posterior means or final states.

## Study 05 operator-validation extension

`R/metrics-operator.R` owns reusable matrix and deterministic-action
comparisons. `benchmark-spec.R` owns the `ld_operator` spec, coordinate grids,
and historical seed arithmetic. `benchmark-execution.R` supplies the explicit
task/validation boundary, and `benchmark-reporting.R` supplies operator design
summaries and named error, rank, spectral, recovery, and runtime plots.

Block construction, retained-eigen policy, projected residual/SSE identities,
and integrated SBayesR audits remain Study 05-specific. No operator registry or
class was introduced. Paused Study 08 sources the one authoritative operator
design file.

Study 01 adds fine-mapping task validation and execution, shared PIP/marker
extraction, fine-mapping metrics, and tidy plotting helpers. Separated causal
locus selection remains in `studies/01_finemapping/locus-design.R` because it
is a Study 01 scientific choice.

All helpers are ordinary R functions. The small design-table and plotting
surface used by downloadable workflows is exported; lower-level mechanics stay
internal. They introduce no class system, registry, plugin interface, workflow
language, or replacement execution engine.

## Study 02 prediction framework

- `R/benchmark-spec.R`: reads and validates ordinary-list specifications,
  resolves profiles, and constructs coordinates and exact seeds.
- `R/benchmark-data.R`: validates/aligned example data, fixed splits, training
  scaling, LD, summary statistics, and `Glist` structure.
- `R/benchmark-simulation.R`: deterministic Study 02-compatible sparse and
  mixture effects, heritability scaling, phenotypes, truth, and oracle checks.
- `R/benchmark-methods.R`: translates the four Study 02 method lists into exact
  BED/CSR calls and controls; it is not a registry.
- `R/benchmark-extraction.R`: extracts effects, probabilities, variance
  components, runtime, metadata, and honest chain information.
- `R/metrics-prediction.R`: prediction/effect metrics, runtime, summaries, and
  paired comparisons as data frames.
- `R/benchmark-execution.R`: authoritative `run_benchmark()` for fine-mapping,
  prediction, parameter estimation, convergence, and LD-operator validation;
  validation, standard output paths, science-identity checkpoints, dispatch,
  extraction, metrics, and provenance. Unsupported task types fail explicitly.

Current callers are the Study 01--05 analyses and CLI, all five reusable
templates, focused tests, and narrowly updated development-study callers.
Annotation and multitrait interfaces remain deferred. There are no Study 02
compatibility wrappers; removal is complete because all callers are internal.

## Study 03 parameter estimation

`benchmark-spec.R` and `benchmark-execution.R` explicitly support the
`parameter_estimation` task. Full-sample data and simulations reuse the Study 02
mechanics. `benchmark-extraction.R` owns scalar trace extraction and draw-wise
nonlinear transformations. `metrics-parameter-estimation.R` owns the frozen
bias, error, coverage, RMSE, MAE, and paired-comparison definitions.

## Study 04 convergence

`benchmark-spec.R` derives the two-stage matched grid from the authoritative
Study 03 spec. `benchmark-extraction.R` accepts only identifiable native chain
traces. `benchmark-convergence.R` owns exact window slicing, R-hat, bulk/tail
ESS, MCSE, relative MCSE, burn-in stability, recommendation selection, and
five-replicate support summaries. `benchmark-execution.R` runs selection then
validation through semantic checkpoints; validation-only resolves 4 workshop
or 24 complete benchmark coordinates without data preparation or fit dispatch.

## `R/benchmark-reporting.R`

- Responsibilities: method/architecture labels and order, factors, plot
  scales/theme, number/runtime/interval/status formatting, replicate summaries,
  compact design/coordinate/method/estimand/output tables, tidy prediction and
  parameter plots, and display of compact capsule scripts.
- Authoritative functions: `sblrbench_method_*`,
  `sblrbench_architecture_*`, `format_sblrbench_*`, `theme_sblrbench()`,
  `prepare_sblrbench_replicates()`, `benchmark_*_table()`,
  `benchmark_data_summary()`, `benchmark_output_inventory()`, the focused
  `plot_*()` workflow helpers, and `display_capsule_script*()`.
- Compatibility: none; the former report loader was removed after all callers
  moved to the authoritative implementation.
- Current callers: Studies 01--06 standard reports, the Study 02 and Study 03
  exact workflows and templates, and reporting/workflow tests.
- Deferred callers: future reports after individual study migration.
- Limitation: plotting consumes extracted tidy tables only; native fits are
  intentionally not accepted.
- Wrapper removal: complete.

## `R/benchmark-provenance.R`

- Responsibilities: Git SHA discovery, installed-package version/path/SHA,
  strict expected-SHA checks, canonical MD5, file SHA-256, session information.
- Authoritative functions: `benchmark_git_sha()`,
  `benchmark_package_provenance()`, `benchmark_assert_package_sha()`,
  `benchmark_canonical_md5()`, `benchmark_file_sha256()`, and
  `benchmark_session_information()`.
- Compatibility: public `sblrbench_git_commit()` delegates to the Git helper;
  no repository-only provenance loader remains.
- Current callers: package provenance, active Studies 01--06 promotion, and
  focused tests.
- Deferred callers: future Study 06 and Study 08 migrations.
- Limitation: installed metadata only; no network lookup.
- Wrapper removal: at the corresponding study migration, after callers use the
  authoritative names.

## `R/benchmark-checkpoints.R`

- Responsibilities: serialized input hashing, atomic RDS replacement, safe RDS
  load, strict hash comparison, and caller-supplied semantic validation.
- Authoritative functions: `benchmark_hash_object()`,
  `benchmark_atomic_save_rds()`, `benchmark_load_checkpoint()`,
  `benchmark_semantic_checkpoint_identity()`,
  `benchmark_semantic_checkpoint_hash()`, and
  `benchmark_load_semantic_checkpoint()`.
- Compatibility: no LD-operator checkpoint wrapper remains. New Study 05 work uses
  semantic-v2 identities; frozen capsules remain the numerical authority.
- Current callers: migrated Studies 01--05, diagnostics, and focused tests.
- Current callers also include the Study 03 SBayesR comparison and integrated
  Study 05 SBayesR evidence. Their semantic schema
  is `sblrbench-semantic-v2`; legacy source-hashed caches are rejected.
- Deferred callers: paused Study 08 and Study 06 annotation CSV status checkpoints.
- Limitation: diagnostic IDs and scientific payload construction remain local.
- Wrapper removal: complete for migrated Studies 01--05.

## `R/benchmark-convergence.R`

- Responsibilities: retained chain windows, trace-array-to-long extraction,
  rank-normalized R-hat, bulk/tail ESS, mean MCSE, relative MCSE, and compact
  threshold flags.
- Authoritative functions: `benchmark_chain_window()`,
  `benchmark_trace_array_long()`, and `benchmark_scalar_diagnostics()`.
- Compatibility: `.study06_one_diagnostic()` adapts annotation output details;
  Study 05 operator diagnostics use the integrated shared/operator helpers.
- Current callers: Study 04 execution and exact workflow, Studies 05--06
  diagnostic wrappers, and focused tests.
- Deferred callers: paused Study 08 and later study-specific trace extraction.
- Limitation: accepts true numeric chain traces only; it does not reconstruct
  traces from posterior means, compact summaries, or final states.
- Wrapper removal: after each study’s trace schema is migrated and compared to
  frozen diagnostic fixtures.

## `R/benchmark-capsules.R`

- Responsibilities: canonical checksum inventory, unique staging directory,
  selected-file copying, validator-first atomic promotion.
- Authoritative functions: `benchmark_capsule_checksums()`,
  `benchmark_capsule_staging_directory()`,
  `benchmark_copy_capsule_files()`, and `benchmark_promote_capsule()`.
- Compatibility: study promotion functions retain destinations, required-file
  lists, manifests, README text, and semantic checks; Study 05--06 checksum
  functions delegate to the shared inventory.
- Current callers: promotion checksum wrappers and focused tests.
- Deferred callers: paused Study 08 and full staging/finalization adoption
  during each study migration, because current cleanup and overwrite policies
  differ.
- Limitation: never chooses scientific outputs or promotes on its own.
- Wrapper removal: once a migrated study supplies destination, selected files,
  semantic validator, title, and description directly.

## `R/benchmark-validation.R`

- Responsibilities: required-file presence, checksum schema/path safety,
  canonical MD5 validation, and mutation detection.
- Authoritative function: `benchmark_validate_capsule_checksums()`.
- Compatibility: Studies 05, 06, and 06 v2 retain named wrapper errors and
  semantic validation after this structural check.
- Current callers: those promotion validators and focused tests.
- Deferred callers: Studies 01--04 and paused Study 08, whose existing combined
  validators are retained until migration.
- Limitation: structural integrity only; it cannot decide scientific
  completeness.
- Wrapper removal: after study-specific semantic validators accept a shared
  structural-validation result directly.
