# Shared benchmark helpers after Phase 3

All helpers are ordinary internal R functions. They introduce no class system,
registry, plugin interface, workflow language, or replacement execution engine.

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
- `R/benchmark-execution.R`: authoritative prediction-only `run_benchmark()`;
  validation, standard output paths, science-identity checkpoints, dispatch,
  extraction, metrics, and provenance. Unsupported task types fail explicitly.

Current callers are Study 02's analysis and CLI, the reusable prediction
template, focused tests, and narrowly updated Study 05--07 data callers.
Parameter-estimation, fine-mapping, annotation, operator, and MTBLR interfaces
remain deferred until those studies migrate. There are no Study 02
compatibility wrappers; removal is complete because all callers are internal.

## Study 03 parameter estimation

`benchmark-spec.R` and `benchmark-execution.R` explicitly support the
`parameter_estimation` task. Full-sample data and simulations reuse the Study 02
mechanics. `benchmark-extraction.R` owns scalar trace extraction and draw-wise
nonlinear transformations. `metrics-parameter-estimation.R` owns the frozen
bias, error, coverage, RMSE, MAE, and paired-comparison definitions.

## `R/benchmark-reporting.R`

- Responsibilities: method/architecture labels and order, factors, plot
  scales/theme, number/runtime/interval/status formatting, replicate summaries,
  and display of compact capsule scripts.
- Authoritative functions: `sblrbench_method_*`,
  `sblrbench_architecture_*`, `format_sblrbench_*`, `theme_sblrbench()`,
  `prepare_sblrbench_replicates()`, and `display_capsule_script*()`.
- Compatibility: `studies/reporting_helpers.R` sources this file and contains
  no implementation.
- Current callers: Studies 01--06 standard reports and reporting tests.
- Deferred callers: future reports after individual study migration.
- Limitation: table-only presentation; native fits are intentionally rejected
  by convention rather than accepted as inputs.
- Wrapper removal: after every report loads package/shared helpers directly.

## `R/benchmark-provenance.R`

- Responsibilities: Git SHA discovery, installed-package version/path/SHA,
  strict expected-SHA checks, canonical MD5, file SHA-256, session information.
- Authoritative functions: `benchmark_git_sha()`,
  `benchmark_package_provenance()`, `benchmark_assert_package_sha()`,
  `benchmark_canonical_md5()`, `benchmark_file_sha256()`, and
  `benchmark_session_information()`.
- Compatibility: public `sblrbench_git_commit()` delegates to the Git helper;
  `.five_replicate_sblr_provenance()` adapts the generic record to its existing
  field names.
- Current callers: package provenance, active Studies 01--06 promotion, and
  focused tests.
- Deferred callers: current-refresh and Study 06 v2 detailed provenance where
  compiler/archive fields remain study-specific.
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
- Compatibility: `.study06_atomic_rds()` and `.study06v2_atomic_rds()` preserve
  current signatures; Study 06 v2 preserves its exact identity payload and
  delegates only digest/save/load mechanics.
- Current callers: Study 06, Study 06 v2, and focused tests.
- Current callers also include the Study 03 SBayesR comparison and Study 06
  supplemental scheduler and exact/sparse diagnostics. Their semantic schema
  is `sblrbench-semantic-v2`; legacy source-hashed caches are rejected.
- Deferred callers: paused Study 07 and Study 05 CSV status checkpoints.
- Limitation: diagnostic IDs and scientific payload construction remain local.
- Wrapper removal: after each study migration demonstrates byte-compatible
  checkpoint names and identities.

## `R/benchmark-convergence.R`

- Responsibilities: retained chain windows, trace-array-to-long extraction,
  rank-normalized R-hat, bulk/tail ESS, mean MCSE, relative MCSE, and compact
  threshold flags.
- Authoritative functions: `benchmark_chain_window()`,
  `benchmark_trace_array_long()`, and `benchmark_scalar_diagnostics()`.
- Compatibility: `.study05_one_diagnostic()` and `.study06_diagnostic_one()`
  adapt output details while retaining study thresholds and historical labels.
- Current callers: Studies 05--06 diagnostics and focused tests.
- Deferred callers: Study 04 registry diagnostics, paused Study 07, and full
  study-specific trace extraction.
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
- Deferred callers: paused Study 07 and full staging/finalization adoption
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
- Deferred callers: Studies 01--04 and paused Study 07, whose existing combined
  validators are retained until migration.
- Limitation: structural integrity only; it cannot decide scientific
  completeness.
- Wrapper removal: after study-specific semantic validators accept a shared
  structural-validation result directly.
