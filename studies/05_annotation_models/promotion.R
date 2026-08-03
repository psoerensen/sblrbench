.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-provenance.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-capsules.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-validation.R"), local = TRUE)

.study05_required <- function(type = c("convergence", "benchmark")) {
  type <- match.arg(type)
  common <- c("README.md", "benchmark_manifest.json", "config.R",
    "example_data_manifest.csv", "source_files.csv",
    "interface_audit_sources.csv", "session_info.txt",
    "reproduce.R", "contract_smoke_test.R", "checksums.csv")
  if (type == "convergence") c(common,
    "annotation_design_summary.csv", "true_alpha.csv", "fit_status.csv",
    "scalar_chain_draws.csv", "convergence_diagnostics.csv",
    "candidate_settings.csv", "method_recommendations.csv",
    "computational_summary.csv", "seed_registry.csv")
  else c(common,
    "annotation_design.csv", "annotation_design_summary.csv",
    "annotation_column_metadata.csv", "sample_split.csv", "true_alpha.csv",
    "simulation_summary.csv", "simulation_marker_truth.csv",
    "fit_status.csv", "prediction_metrics.csv",
    "prediction_summary.csv", "parameter_estimates.csv",
    "parameter_recovery_summary.csv", "annotation_coefficient_estimates.csv",
    "annotation_recovery_summary.csv", "annotation_variance_summary.csv",
    "marker_effect_metrics.csv", "probability_recovery_metrics.csv",
    "probability_recovery_summary.csv", "probability_output_inventory.csv",
    "enrichment_metrics.csv", "enrichment_summary.csv",
    "paired_replicate_differences.csv", "paired_comparison_summary.csv",
    "annotation_value_interactions.csv",
    "annotation_value_interaction_summary.csv",
    "convergence_diagnostics.csv", "convergence_validation_summary.csv",
    "computational_summary.csv", "seed_registry.csv",
    "method_output_availability.csv")
}

.study05_checksums <- function(path) {
  benchmark_capsule_checksums(path)
}

.study05_validate_checksums <- function(path, required) {
  tryCatch(benchmark_validate_capsule_checksums(path, required),
    error = function(error) stop("Study 05 capsule checksum validation failed: ",
      conditionMessage(error), call. = FALSE))
}

.study05_validate_convergence_capsule <- function(path) {
  required <- .study05_required("convergence")
  if (!dir.exists(path) || length(setdiff(required, list.files(path))))
    stop("Study 05 convergence capsule is incomplete.", call. = FALSE)
  .study05_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  status <- read.csv(file.path(path, "fit_status.csv"))
  rec <- read.csv(file.path(path, "method_recommendations.csv"))
  diagnostics <- read.csv(file.path(path, "convergence_diagnostics.csv"))
  draws <- read.csv(file.path(path, "scalar_chain_draws.csv"))
  seeds <- read.csv(file.path(path, "seed_registry.csv"))
  expected <- c("st_bed_bayesrc", "st_csr_sbayesrc")
  if (!identical(manifest$task, "annotation_model_convergence_selection") ||
      !identical(manifest$benchmark_status, "complete") ||
      !identical(manifest$installed_sblr_commit,
        "02e8c74baa906e83c4a08d42a9cc6339b4e81072") ||
      manifest$expected_fit_count != 2L ||
      nrow(status) != 2L || !setequal(status$method, expected) ||
      any(status$status != "ok") || any(status$chain_count != 4L) ||
      nrow(rec) != 2L || !setequal(rec$method, expected) ||
      any(rec$recommendation_status != "available") ||
      any(rec$nchains != 4L) || any(rec$ncores != 4L) ||
      any(rec$nthin != 1L) || any(rec$nburn + rec$nit > 3000L) ||
      nrow(seeds) != 8L || anyDuplicated(seeds[c("method", "chain")]) ||
      anyNA(seeds$effective_chain_seed) ||
      any(vapply(split(seeds$chain_seed, seeds$method),
        function(x) length(unique(x)) != 4L, logical(1))) ||
      any(vapply(split(seeds$effective_chain_seed, seeds$method),
        function(x) length(unique(x)) != 4L, logical(1))) ||
      any(!is.finite(diagnostics$rhat)) ||
      any(!is.finite(diagnostics$ess_bulk)) ||
      any(!is.finite(diagnostics$ess_tail)) ||
      any(!is.finite(diagnostics$relative_mcse)) ||
      !all(c("alpha", "sigmaSqAlpha") %in% unique(draws$parameter_name)) ||
      !identical(sort(unique(draws$chain)), 1:4))
    stop("Study 05 convergence capsule semantic validation failed.",
      call. = FALSE)
  invisible(TRUE)
}

.study05_validate_stop_capsule <- function(path) {
  required <- .study05_required("convergence")
  if (!dir.exists(path) || length(setdiff(required, list.files(path))))
    stop("Study 05 stop capsule is incomplete.", call. = FALSE)
  .study05_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  status <- read.csv(file.path(path, "fit_status.csv"), stringsAsFactors = FALSE)
  rec <- read.csv(file.path(path, "method_recommendations.csv"), stringsAsFactors = FALSE)
  diagnostics <- read.csv(file.path(path, "convergence_diagnostics.csv"),
    stringsAsFactors = FALSE)
  seeds <- read.csv(file.path(path, "seed_registry.csv"), stringsAsFactors = FALSE)
  expected <- c("st_bed_bayesrc", "st_csr_sbayesrc")
  if (!identical(manifest$benchmark_status, "stopped") ||
      !identical(manifest$validation_status,
        "prespecified_convergence_stop_triggered") ||
      !identical(manifest$installed_sblr_commit,
        "02e8c74baa906e83c4a08d42a9cc6339b4e81072") ||
      nrow(status) != 2L || !setequal(status$method, expected) ||
      any(status$status != "ok") || any(status$chain_count != 4L) ||
      nrow(rec) != 2L || !setequal(rec$method, expected) ||
      all(rec$recommendation_status == "available") ||
      nrow(seeds) != 8L || anyDuplicated(seeds[c("method", "chain")]) ||
      any(!is.finite(diagnostics$rhat)) ||
      any(!is.finite(diagnostics$ess_bulk)) ||
      any(!is.finite(diagnostics$ess_tail)) ||
      any(!is.finite(diagnostics$relative_mcse)))
    stop("Study 05 stop capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}

.study05_validate_benchmark_capsule <- function(path) {
  required <- .study05_required("benchmark")
  if (!dir.exists(path) || length(setdiff(required, list.files(path))))
    stop("Study 05 benchmark capsule is incomplete.", call. = FALSE)
  .study05_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  status <- read.csv(file.path(path, "fit_status.csv"))
  annotation <- read.csv(file.path(path, "annotation_design.csv"))
  marker_truth <- read.csv(file.path(path, "simulation_marker_truth.csv"))
  inventory <- read.csv(file.path(path, "probability_output_inventory.csv"))
  alpha <- read.csv(file.path(path, "annotation_coefficient_estimates.csv"))
  diagnostics <- read.csv(file.path(path, "convergence_diagnostics.csv"))
  paired <- read.csv(file.path(path, "paired_replicate_differences.csv"))
  methods <- c("st_bed_bayesr", "st_bed_bayesrc",
    "st_csr_sbayesr", "st_csr_sbayesrc")
  scenarios <- c("informative_annotations", "uninformative_annotations")
  expected <- expand.grid(scenario = scenarios, replicate = 1:5,
    method = methods, stringsAsFactors = FALSE)
  key <- function(x) paste(x$scenario, x$replicate, x$method, sep = "|")
  required_inventory <- inventory$required
  if (!identical(manifest$task, "single_trait_annotation_informed_models") ||
      !identical(manifest$benchmark_status, "complete") ||
      manifest$replicate_count != 5L || manifest$expected_fit_count != 40L ||
      manifest$successful_fit_count != 40L || manifest$failed_fit_count != 0L ||
      manifest$chains_per_fit != 4L || manifest$expected_chain_count != 160L ||
      !identical(as.numeric(manifest$mixture_var), c(0, .01, .1, 1)) ||
      nrow(status) != 40L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || any(status$status != "ok") ||
      ncol(annotation) != 5L ||
      nrow(annotation) != manifest$marker_count ||
      !identical(names(annotation), c("marker_id", "Intercept",
        "enriched_binary", "continuous_signal", "null_annotation")) ||
      anyDuplicated(annotation$marker_id) || any(!is.finite(as.matrix(annotation[-1L]))) ||
      nrow(marker_truth) != manifest$marker_count * 10L ||
      anyDuplicated(marker_truth[c("scenario", "replicate", "marker_id")]) ||
      any(!is.finite(as.matrix(marker_truth[paste0(
        "true_prior_component_", 0:3)]))) ||
      any(abs(rowSums(marker_truth[paste0(
        "true_prior_component_", 0:3)]) - 1) > 1e-8) ||
      nrow(inventory) != 80L ||
      any(inventory$row_count[required_inventory] != manifest$marker_count) ||
      any(inventory$column_count[required_inventory] != 4L) ||
      any(inventory$minimum_probability[required_inventory] < 0) ||
      any(inventory$maximum_probability[required_inventory] > 1) ||
      any(inventory$maximum_row_sum_error[required_inventory] > 1e-8) ||
      nrow(alpha) != 240L ||
      any(!is.finite(alpha$posterior_mean)) ||
      length(unique(interaction(diagnostics$scenario, diagnostics$replicate,
        diagnostics$method, drop = TRUE))) != 40L ||
      any(diagnostics$chain_count != 4L) ||
      any(!is.finite(diagnostics$rhat)) ||
      any(!is.finite(diagnostics$ess_bulk)) ||
      any(!is.finite(diagnostics$ess_tail)) ||
      any(!is.finite(diagnostics$relative_mcse)) ||
      any(!paired$complete_pair) ||
      any(table(interaction(paired$comparison_id, paired$scenario,
        paired$metric, drop = TRUE)) != 5L) ||
      !file.exists(file.path(path, "source_files.csv")) ||
      !file.exists(file.path(path, "session_info.txt")))
    stop("Study 05 benchmark capsule semantic validation failed.",
      call. = FALSE)
  invisible(TRUE)
}

.study05_interface_audit_table <- function(config) data.frame(
  installed_commit = packageDescription("sblr")$RemoteSha,
  consulted_path = config$source_audit,
  access = "read_only_git_show_or_read_only_worktree",
  sibling_modified = FALSE, stringsAsFactors = FALSE)

.study05_promote <- function(type = c("convergence", "benchmark"),
                             source_dir, destination, config,
                             validator_override = NULL) {
  type <- match.arg(type)
  validator <- if (!is.null(validator_override)) validator_override else if (type == "convergence")
    .study05_validate_convergence_capsule else
      .study05_validate_benchmark_capsule
  if (dir.exists(destination)) {
    validator(destination)
    return(invisible(destination))
  }
  required <- .study05_required(type)
  generated <- setdiff(required, c("README.md", "config.R",
    "example_data_manifest.csv", "source_files.csv",
    "interface_audit_sources.csv", "session_info.txt", "reproduce.R",
    "contract_smoke_test.R", "checksums.csv"))
  missing <- generated[!file.exists(file.path(source_dir, generated))]
  if (length(missing)) stop("Study 05 promotion source is incomplete: ",
    paste(missing, collapse = ", "), call. = FALSE)
  staging <- file.path(Sys.getenv("SBLR_BENCH_STUDY05_LOCAL",
    file.path("results", "local", "study05_annotation_models")),
    "promotion_staging", paste0(basename(destination), "-", Sys.getpid()))
  if (dir.exists(staging)) stop("Study 05 promotion staging already exists.",
    call. = FALSE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(file.path(source_dir, generated), staging,
    overwrite = FALSE))) stop("Study 05 output copy failed.", call. = FALSE)
  file.copy("studies/05_annotation_models/config.R",
    file.path(staging, "config.R"))
  reproduce_source <- if (nzchar(Sys.getenv("SBLR_BENCH_REFRESH_ROOT", "")))
    "scripts/run_current_benchmark_refresh.R" else
      "scripts/run_study05_annotation_models.R"
  file.copy(reproduce_source,
    file.path(staging, "reproduce.R"))
  file.copy("studies/05_annotation_models/contract_smoke_test.R",
    file.path(staging, "contract_smoke_test.R"))
  file.copy(file.path("results", "reference", "02_prediction",
    "current",
    "example_data_manifest.csv"), file.path(staging,
      "example_data_manifest.csv"))
  source_files <- c("studies/05_annotation_models/config.R",
    "studies/05_annotation_models/annotation_design.R",
    "studies/05_annotation_models/simulation.R",
    "studies/05_annotation_models/methods.R",
    "studies/05_annotation_models/chain_extraction.R",
    "studies/05_annotation_models/diagnostics.R",
    "studies/05_annotation_models/metrics.R",
    "studies/05_annotation_models/pilot.R",
    "studies/05_annotation_models/targets.R",
    "studies/05_annotation_models/promotion.R",
    "studies/05_annotation_models/interface_fit_smoke_test.R",
    reproduce_source)
  write.csv(data.frame(file = source_files,
    md5 = unname(benchmark_canonical_md5(source_files))),
    file.path(staging, "source_files.csv"), row.names = FALSE)
  write.csv(.study05_interface_audit_table(config),
    file.path(staging, "interface_audit_sources.csv"), row.names = FALSE)
  writeLines(sub("[[:space:]]+$", "", capture.output(utils::sessionInfo())),
    file.path(staging, "session_info.txt"))
  writeLines(c(
    paste0("# Study 05: ", if (type == "convergence")
      "annotation-model convergence development" else
        "five-replicate annotation-model development benchmark"),
    "", "Package-specific development evidence for the current installed sblr implementation.",
    "The explicit component grid is c(0, 0.01, 0.1, 1).",
    "Five simulations are descriptive evidence, not universal validation.",
    "Native fit objects and target caches are deliberately excluded."),
    file.path(staging, "README.md"))
  write.csv(.study05_checksums(staging),
    file.path(staging, "checksums.csv"), row.names = FALSE)
  validator(staging)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging, destination))
    stop("Atomic Study 05 capsule promotion failed.", call. = FALSE)
  validator(destination)
  invisible(destination)
}

.study05_promote_convergence <- function(source_dir, config) .study05_promote(
  "convergence", source_dir, config$convergence_capsule, config)

.study05_promote_benchmark <- function(source_dir, config) .study05_promote(
  "benchmark", source_dir, config$benchmark_capsule, config)

.study05_promote_current_decision <- function(source_dir, config) {
  rec <- utils::read.csv(file.path(source_dir, "method_recommendations.csv"),
    stringsAsFactors = FALSE)
  supported <- nrow(rec) == 2L && all(rec$recommendation_status == "available")
  destination <- file.path("results", "reference", "05_annotation_models",
    if (supported) "current-convergence" else "current-stop")
  manifest_path <- file.path(source_dir, "benchmark_manifest.json")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  manifest$benchmark_status <- if (supported) "complete" else "stopped"
  manifest$validation_status <- if (supported) "complete_grid_validated" else
    "prespecified_convergence_stop_triggered"
  manifest$full_benchmark_started <- FALSE
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE,
    auto_unbox = TRUE, null = "null")
  config$convergence_capsule <- destination
  .study05_promote("convergence", source_dir, destination, config,
    validator_override = if (supported) .study05_validate_convergence_capsule else
      .study05_validate_stop_capsule)
  if (supported) .study05_validate_convergence_capsule(destination) else
    .study05_validate_stop_capsule(destination)
  list(destination = destination, supported = supported)
}
