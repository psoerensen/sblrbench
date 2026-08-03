.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-provenance.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-capsules.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-validation.R"), local = TRUE)

.study06v2_capsule_required <- function(type = c("convergence", "benchmark")) {
  type <- match.arg(type)
  common <- c("README.md", "benchmark_manifest.json", "config.R",
    "example_data_manifest.csv", "source_files.csv", "session_info.txt",
    "reproduce.R", "contract_smoke_test.R", "checksums.csv",
    "v1_to_v2_design_crosswalk.csv", "deterministic_identity_summary.csv",
    "low_rank_block_diagnostics.csv", "deterministic_validation_manifest.json")
  if (type == "convergence") c(common,
    "historical_v1_inventory.md", "optimized_hash_inventory.csv",
    "one_replicate_operator_pilot.csv", "pilot_gate.csv",
    "pilot_paired_comparisons.csv", "pilot_runtime_comparison.csv",
    "projected_complete_grid_runtime.csv", "historical_runtime_context.csv",
    "fit_status.csv",
    "candidate_settings.csv", "convergence_diagnostics.csv",
    "method_recommendations.csv")
  else c(common, "block_definitions.csv", "seed_registry.csv",
    "checkpoint_validation.csv",
    "simulation_summary.csv", "fit_status.csv", "prediction_metrics.csv",
    "parameter_estimates.csv", "parameter_recovery_summary.csv",
    "marker_effect_metrics.csv", "variable_selection_metrics.csv",
    "marker_effect_agreement.csv", "sbayesr_diagnostic_evidence.csv",
    "paired_replicate_differences.csv", "paired_comparison_summary.csv",
    "convergence_diagnostics.csv", "convergence_validation_summary.csv",
    "computational_summary.csv", "method_output_availability.csv",
    "low_rank_fit_block_diagnostics.csv", "method_recommendations.csv",
    "readiness_decision.txt")
}

.study06v2_checksums <- function(path) {
  benchmark_capsule_checksums(path)
}

.study06v2_validate_checksums <- function(path, required) {
  tryCatch(benchmark_validate_capsule_checksums(path, required),
    error = function(error) stop("Study 06 v2 canonical checksum validation failed: ",
      conditionMessage(error), call. = FALSE))
}

.study06v2_validate_capsule <- function(path,
                                        type = c("convergence", "benchmark"),
                                        config) {
  type <- match.arg(type)
  required <- .study06v2_capsule_required(type)
  if (!dir.exists(path) || length(setdiff(required, list.files(path))))
    stop("Study 06 v2 capsule is incomplete.", call. = FALSE)
  .study06v2_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  status <- read.csv(file.path(path, "fit_status.csv"),
    stringsAsFactors = FALSE)
  expected <- if (type == "convergence") config$convergence_fit_count else
    config$expected_fit_count
  if (!identical(manifest$study, config$study) ||
      !identical(manifest$sblr_source_sha, config$required_sblr_sha) ||
      !identical(manifest$sblr_version, config$required_sblr_version) ||
      !identical(manifest$operator_contract, config$operator_contract) ||
      !identical(manifest$low_rank_representation, config$representation) ||
      manifest$expected_fit_count != expected || nrow(status) != expected ||
      any(status$status != "ok") || any(status$chain_count != 4L))
    stop("Study 06 v2 capsule provenance or fit-grid validation failed.",
      call. = FALSE)
  status_key <- paste(status$architecture, status$replicate,
    status$configuration, sep = "::")
  expected_grid <- if (type == "convergence")
    expand.grid(architecture = config$architectures, replicate = 1L,
      configuration = setdiff(config$configurations, "bed"),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) else
    expand.grid(architecture = config$architectures,
      replicate = seq_len(config$replicate_count),
      configuration = config$configurations, KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE)
  expected_key <- paste(expected_grid$architecture, expected_grid$replicate,
    expected_grid$configuration, sep = "::")
  if (anyDuplicated(status_key) || !setequal(status_key, expected_key))
    stop("Study 06 v2 capsule fit coordinates are incomplete or duplicated.",
      call. = FALSE)
  if (type == "benchmark") {
    paired <- read.csv(file.path(path, "paired_comparison_summary.csv"),
      stringsAsFactors = FALSE)
    low <- read.csv(file.path(path, "low_rank_fit_block_diagnostics.csv"),
      stringsAsFactors = FALSE)
    checkpoint <- read.csv(file.path(path, "checkpoint_validation.csv"),
      stringsAsFactors = FALSE)
    recommendations <- read.csv(file.path(path, "method_recommendations.csv"),
      stringsAsFactors = FALSE)
    checkpoint_key <- paste(checkpoint$architecture, checkpoint$replicate,
      checkpoint$configuration, sep = "::")
    prediction <- read.csv(file.path(path, "prediction_metrics.csv"),
      stringsAsFactors = FALSE)
    low_checkpoint <- checkpoint$operator_representation == "low_rank"
    recommended_checkpoint <- checkpoint[checkpoint$configuration != "bed",
      , drop = FALSE]
    recommended_checkpoint <- merge(recommended_checkpoint,
      recommendations[c("architecture", "configuration", "nburn", "nit",
        "nthin", "nchains", "ncores")],
      by = c("architecture", "configuration"), all.x = TRUE,
      suffixes = c("_observed", "_recommended"))
    controls_match <- nrow(recommended_checkpoint) ==
      config$expected_fit_count - config$replicate_count *
        length(config$architectures) &&
      all(with(recommended_checkpoint,
        nburn_observed == nburn_recommended &
        retained_iterations == nit & nthin_observed == nthin_recommended &
        chain_count == nchains & ncores_observed == ncores_recommended))
    if (manifest$replicate_count != 5L ||
        any(paired$complete_paired_replicates != 5L) ||
        nrow(checkpoint) != config$expected_fit_count ||
        anyDuplicated(checkpoint_key) || !setequal(checkpoint_key, expected_key) ||
        any(checkpoint$status != "ok") ||
        any(checkpoint$package_version != config$required_sblr_version) ||
        any(checkpoint$package_commit != config$required_sblr_sha) ||
        !controls_match ||
        any(!is.finite(prediction$value)) ||
        !setequal(unique(low$configuration),
          c("low_rank_full", "low_rank_0999", "low_rank_0995")) ||
        any(low$representation != "low_rank") ||
        any(checkpoint$operator_representation == "dense_reconstructed") ||
        any(!is.finite(checkpoint$maximum_residual_drift[low_checkpoint])) ||
        any(checkpoint$maximum_residual_drift[low_checkpoint] >
          config$operator_tolerance$product_absolute) ||
        any(checkpoint$retained_rank[low_checkpoint] >
          checkpoint$positive_rank[low_checkpoint]) ||
        !readLines(file.path(path, "readiness_decision.txt"), warn = FALSE)[1L] %in%
          c("READY FOR MTBLR LOW-RANK EXTENSION",
            "NOT READY FOR MTBLR LOW-RANK EXTENSION"))
      stop("Study 06 v2 benchmark capsule semantic validation failed.",
        call. = FALSE)
  }
  text_files <- list.files(path, pattern = "\\.(csv|json|md|txt|R)$",
    full.names = TRUE)
  leaked <- vapply(text_files, function(file) any(grepl(
    "[A-Za-z]:[/\\\\](Users|Documents|AppData)[/\\\\]",
    readLines(file, warn = FALSE), ignore.case = TRUE)), logical(1L))
  if (any(leaked))
    stop("Study 06 v2 capsule contains an absolute local path: ",
      paste(basename(text_files[leaked]), collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.study06v2_manifest <- function(type, status, config, readiness = NULL) {
  expected <- if (type == "convergence") config$convergence_fit_count else
    config$expected_fit_count
  list(study = config$study,
    task = if (type == "convergence")
      "single_trait_retained_low_rank_convergence_selection" else config$task,
    benchmark_scope = if (type == "convergence")
      "retained_low_rank_convergence_development_v2" else
        "retained_low_rank_five_replicate_development_v2",
    benchmark_status = "complete", expected_fit_count = expected,
    successful_fit_count = sum(status$status == "ok"),
    failed_fit_count = sum(status$status != "ok"),
    replicate_count = if (type == "benchmark") config$replicate_count else 1L,
    architectures = config$architectures,
    configurations = if (type == "convergence")
      setdiff(config$configurations, "bed") else config$configurations,
    chains_per_fit = 4L, sblr_source_sha = config$required_sblr_sha,
    sblr_version = config$required_sblr_version,
    sblrbench_version = as.character(packageVersion("sblrbench")),
    repository_commit = trimws(system2("git", c("rev-parse", "HEAD"),
      stdout = TRUE)), qgdata_commit = config$example_data$commit,
    operator_contract = config$operator_contract,
    low_rank_representation = config$representation,
    eigen_props = list(near_full = config$eigen_prop_full,
      high_retention = unname(config$eigen_props[["low_rank_0999"]]),
      canonical = unname(config$eigen_props[["low_rank_0995"]])),
    projected_variance_contract =
      "yy - sum_b crossprod(w_b) + sum_b crossprod(r_b)",
    historical_v1_label = "Historical reconstructed-dense Study 06 v1",
    readiness_decision = readiness,
    thread_controls = list(OMP_NUM_THREADS = 4L, OMP_THREAD_LIMIT = 4L,
      OPENBLAS_NUM_THREADS = 1L, MKL_NUM_THREADS = 1L,
      VECLIB_MAXIMUM_THREADS = 1L),
    source_status = "pinned_local_installed_package_offline",
    validation_status = "complete_grid_validated",
    failures = unname(status$error_message[status$status != "ok"]))
}

.study06v2_copy_unique <- function(name, source_dirs, destination) {
  candidates <- file.path(source_dirs, name)
  existing <- candidates[file.exists(candidates)]
  if (!length(existing))
    stop("Study 06 v2 promotion source missing or ambiguous: ", name,
      call. = FALSE)
  if (length(existing) > 1L) {
    hashes <- unname(benchmark_canonical_md5(existing))
    if (length(unique(hashes)) != 1L)
      stop("Study 06 v2 promotion sources disagree: ", name,
        call. = FALSE)
  }
  if (!file.copy(existing[1L], file.path(destination, name),
      overwrite = FALSE))
    stop("Study 06 v2 promotion source copy failed: ", name,
      call. = FALSE)
}

.study06v2_readiness_decision <- function(config) {
  deterministic <- read.csv(file.path(config$local_dir, "deterministic",
    "deterministic_identity_summary.csv"), stringsAsFactors = FALSE)
  pilot <- read.csv(file.path(config$local_dir, "operator_pilot",
    "pilot_gate.csv"), stringsAsFactors = FALSE)
  recommendations <- read.csv(file.path(config$local_dir, "convergence",
    "method_recommendations.csv"), stringsAsFactors = FALSE)
  aggregate <- file.path(config$local_dir, "aggregate")
  paired <- read.csv(file.path(aggregate, "paired_replicate_differences.csv"),
    stringsAsFactors = FALSE)
  marker <- read.csv(file.path(aggregate, "marker_effect_agreement.csv"),
    stringsAsFactors = FALSE)
  sensitivity <- paired[paired$comparison_id ==
    "low_rank_0995_minus_low_rank_0999", , drop = FALSE]
  pred <- sensitivity[sensitivity$metric ==
    "phenotype_prediction_correlation", , drop = FALSE]
  parameter <- read.csv(file.path(aggregate, "parameter_estimates.csv"),
    stringsAsFactors = FALSE)
  h2 <- parameter[parameter$estimand == "heritability",
    c("architecture", "replicate", "configuration", "posterior_mean")]
  a <- h2[h2$configuration == "low_rank_0995", ]
  b <- h2[h2$configuration == "low_rank_0999", ]
  h2_pair <- merge(a, b, by = c("architecture", "replicate"),
    suffixes = c("_0995", "_0999"))
  marker_pair <- marker[marker$comparison_id ==
    "low_rank_0995_minus_low_rank_0999", , drop = FALSE]
  pass <- all(deterministic$pass) && isTRUE(pilot$pass[1L]) &&
    nrow(recommendations) == config$convergence_fit_count &&
    all(recommendations$recommendation_status == "available") &&
    nrow(pred) == 2L * config$replicate_count && all(pred$complete_pair) &&
    max(abs(pred$difference)) <= config$readiness_gate$
      maximum_0995_vs_0999_prediction_correlation_difference &&
    nrow(h2_pair) == 2L * config$replicate_count &&
    max(abs(h2_pair$posterior_mean_0995 - h2_pair$posterior_mean_0999)) <=
      config$readiness_gate$maximum_0995_vs_0999_heritability_difference &&
    nrow(marker_pair) == 2L * config$replicate_count &&
    min(marker_pair$posterior_mean_effect_correlation) >=
      config$readiness_gate$minimum_0995_vs_0999_marker_effect_correlation
  if (isTRUE(pass)) "READY FOR MTBLR LOW-RANK EXTENSION" else
    "NOT READY FOR MTBLR LOW-RANK EXTENSION"
}

.study06v2_promote <- function(type = c("convergence", "benchmark"), config,
                               readiness = NULL) {
  type <- match.arg(type)
  destination <- if (type == "convergence") config$convergence_capsule else
    config$benchmark_capsule
  if (dir.exists(destination)) {
    .study06v2_validate_capsule(destination, type, config)
    return(invisible(destination))
  }
  local <- config$local_dir
  source_dirs <- if (type == "convergence") c(
    file.path(local, "deterministic"), file.path(local, "operator_pilot"),
    file.path(local, "convergence"),
    file.path(local, "historical_v1_inventory")) else c(
      file.path(local, "deterministic"), file.path(local, "aggregate"))
  required <- .study06v2_capsule_required(type)
  staging <- file.path(local, "promotion_staging",
    paste0(basename(destination), "-", Sys.getpid()))
  if (dir.exists(staging))
    stop("Study 06 v2 promotion staging already exists.", call. = FALSE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  generated <- setdiff(required, c("README.md", "benchmark_manifest.json",
    "config.R", "example_data_manifest.csv", "source_files.csv",
    "session_info.txt", "reproduce.R", "contract_smoke_test.R",
    "checksums.csv", "historical_v1_inventory.md", "readiness_decision.txt"))
  for (name in generated) {
    dirs <- if (type == "benchmark" &&
        identical(name, "method_recommendations.csv"))
      file.path(local, "convergence") else source_dirs
    .study06v2_copy_unique(name, dirs, staging)
  }
  file.copy(file.path("studies", "06_ld_operator", "v2", "config.R"),
    file.path(staging, "config.R"))
  file.copy("scripts/run_study06_low_rank_operator_v2.R",
    file.path(staging, "reproduce.R"))
  file.copy(file.path("studies", "06_ld_operator", "v2",
    "contract_smoke_test.R"), file.path(staging, "contract_smoke_test.R"))
  file.copy(file.path("results", "reference", "06_ld_operator",
    "current",
    "example_data_manifest.csv"), file.path(staging,
      "example_data_manifest.csv"))
  if (type == "convergence") file.copy(file.path(local,
    "historical_v1_inventory", "README.md"), file.path(staging,
      "historical_v1_inventory.md"))
  if (type == "benchmark") writeLines(readiness,
    file.path(staging, "readiness_decision.txt"))
  sources <- c("scripts/run_study06_low_rank_operator_v2.R",
    "scripts/run_study06_low_rank_operator_v2.ps1",
    list.files(file.path("studies", "06_ld_operator", "v2"),
      full.names = TRUE),
    "studies/06_ld_operator/low-rank-operator.qmd")
  sources <- sort(sources[file.exists(sources)])
  write.csv(data.frame(file = sources,
    md5 = unname(benchmark_canonical_md5(sources))),
    file.path(staging, "source_files.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(staging,
    "session_info.txt"))
  writeLines(c("# Study 06 v2: retained low-rank LD operator validation", "",
    "Five paired development simulations using the pinned sblr 0.2.0",
    "retained low-rank implementation. Historical v1 reconstructed-dense",
    "results remain separate. Native fits, Q matrices, and target stores are excluded."),
    file.path(staging, "README.md"))
  status <- read.csv(file.path(if (type == "convergence")
    file.path(local, "convergence") else file.path(local, "aggregate"),
    "fit_status.csv"), stringsAsFactors = FALSE)
  manifest <- .study06v2_manifest(type, status, config, readiness)
  jsonlite::write_json(manifest, file.path(staging, "benchmark_manifest.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  write.csv(.study06v2_checksums(staging), file.path(staging,
    "checksums.csv"), row.names = FALSE)
  .study06v2_validate_capsule(staging, type, config)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging, destination))
    stop("Atomic Study 06 v2 capsule promotion failed.", call. = FALSE)
  .study06v2_validate_capsule(destination, type, config)
  invisible(destination)
}
