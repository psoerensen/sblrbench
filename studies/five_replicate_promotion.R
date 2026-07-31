source(file.path("studies", "02_prediction", "promotion.R"), local = TRUE)

.five_capsule_checksums <- function(path) {
  files <- sort(setdiff(list.files(path, recursive = FALSE), "checksums.csv"))
  info <- file.info(file.path(path, files))
  data.frame(file = files, size_bytes = info$size,
    md5 = unname(.study02_canonical_md5(file.path(path, files))),
    stringsAsFactors = FALSE)
}

.five_validate_checksum_inventory <- function(path, required) {
  checks <- utils::read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  if (!identical(names(checks), c("file", "size_bytes", "md5")) ||
      anyNA(checks$file) || anyDuplicated(checks$file) ||
      any(checks$file != basename(checks$file)) ||
      any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])\\.\\.([/\\\\]|$))", checks$file)) ||
      !all(setdiff(required, "checksums.csv") %in% checks$file) ||
      any(!grepl("^[0-9a-f]{32}$", checks$md5)))
    stop("Invalid capsule checksum inventory.", call. = FALSE)
  paths <- file.path(path, checks$file)
  if (any(!file.exists(paths))) stop("Checksum inventory lists missing files.", call. = FALSE)
  actual <- unname(.study02_canonical_md5(paths))
  if (any(actual != checks$md5))
    stop("Capsule checksum validation failed: ",
      paste(checks$file[actual != checks$md5], collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.five_required_files <- function(study) switch(study,
  study02 = c("README.md", "benchmark_manifest.json", "prediction_metrics.csv",
    "benchmark_summary.csv", "paired_method_differences.csv",
    "paired_comparison_summary.csv", "computational_summary.csv",
    "replicate_status.csv", "simulation_summary.csv", "seed_registry.csv",
    "example_data_manifest.csv", "source_files.csv", "session_info.txt", "checksums.csv"),
  study03 = c("README.md", "benchmark_manifest.json", "simulation_truth.csv",
    "parameter_estimates.csv", "parameter_recovery_summary.csv",
    "paired_parameter_differences.csv", "paired_comparison_summary.csv",
    "computational_summary.csv", "replicate_status.csv", "seed_registry.csv",
    "estimand_registry.csv", "example_data_manifest.csv", "source_files.csv",
    "session_info.txt", "checksums.csv"),
  study04 = c("README.md", "benchmark_manifest.json", "scalar_chain_draws.csv",
    "convergence_diagnostics.csv", "replicate_diagnostic_summary.csv",
    "method_validation_summary.csv", "chain_summaries.csv", "chain_agreement.csv",
    "computational_summary.csv", "chain_status.csv", "seed_registry.csv",
    "example_data_manifest.csv", "source_files.csv", "session_info.txt", "checksums.csv"))

.five_validate_study02_capsule <- function(path) {
  required <- .five_required_files("study02")
  if (length(setdiff(required, list.files(path)))) stop("Study 02 capsule files are incomplete.", call. = FALSE)
  .five_validate_checksum_inventory(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- read.csv(file.path(path, "replicate_status.csv"), stringsAsFactors = FALSE)
  metrics <- read.csv(file.path(path, "prediction_metrics.csv"), stringsAsFactors = FALSE)
  paired <- read.csv(file.path(path, "paired_method_differences.csv"), stringsAsFactors = FALSE)
  summary <- read.csv(file.path(path, "benchmark_summary.csv"), stringsAsFactors = FALSE)
  expected <- expand.grid(architecture = c("sparse_homogeneous", "sparse_mixture"),
    replicate = 1:5, method = .study02_expected_methods(), stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate, x$method, sep = "|")
  if (manifest$task != "single_trait_prediction" ||
      manifest$benchmark_scope != "five_replicate_development" ||
      manifest$replicate_count != 5L || manifest$expected_fit_count != 40L ||
      manifest$successful_fit_count != 40L || manifest$failed_fit_count != 0L ||
      nrow(status) != 40L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || any(status$status != "ok") ||
      any(!is.finite(metrics$value)) || any(!paired$complete_pair) ||
      any(!is.finite(paired$advantage)) ||
      !all(c("replicate_count", "successful_replicates", "mean", "sd", "median",
        "minimum", "maximum") %in% names(summary)) ||
      any(summary$replicate_count != 5L) || any(summary$successful_replicates != 5L) ||
      any(!is.finite(summary$sd)))
    stop("Study 02 five-replicate capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}

.five_validate_study03_capsule <- function(path) {
  required <- .five_required_files("study03")
  if (length(setdiff(required, list.files(path)))) stop("Study 03 capsule files are incomplete.", call. = FALSE)
  .five_validate_checksum_inventory(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- read.csv(file.path(path, "replicate_status.csv"), stringsAsFactors = FALSE)
  estimates <- read.csv(file.path(path, "parameter_estimates.csv"), stringsAsFactors = FALSE)
  summary <- read.csv(file.path(path, "parameter_recovery_summary.csv"), stringsAsFactors = FALSE)
  expected <- expand.grid(architecture = c("sparse_homogeneous", "sparse_mixture"),
    replicate = 1:5, method = .study02_expected_methods(), stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate, x$method, sep = "|")
  available <- estimates$status == "ok"
  if (manifest$task != "single_trait_parameter_estimation" ||
      manifest$benchmark_scope != "five_replicate_development" ||
      manifest$replicate_count != 5L || manifest$expected_fit_count != 40L ||
      manifest$successful_fit_count != 40L || manifest$failed_fit_count != 0L ||
      nrow(status) != 40L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || any(status$status != "ok") ||
      any(estimates$chain_count[available] != 4L) ||
      any(estimates$draws_per_chain[available] <= 0L) ||
      any(!is.finite(estimates$posterior_mean[available])) ||
      !all(c("mean_bias", "rmse", "mae", "observed_coverage_count",
        "observed_coverage_proportion", "mean_interval_width_95") %in% names(summary)))
    stop("Study 03 five-replicate capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}

.five_validate_study04_capsule <- function(path) {
  required <- .five_required_files("study04")
  if (length(setdiff(required, list.files(path)))) stop("Study 04 capsule files are incomplete.", call. = FALSE)
  .five_validate_checksum_inventory(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- read.csv(file.path(path, "chain_status.csv"), stringsAsFactors = FALSE)
  diagnostics <- read.csv(file.path(path, "convergence_diagnostics.csv"), stringsAsFactors = FALSE)
  summary <- read.csv(file.path(path, "method_validation_summary.csv"), stringsAsFactors = FALSE)
  key <- paste(status$architecture, status$replicate, status$method, status$chain, sep = "|")
  if (manifest$task != "single_trait_multichain_convergence_validation" ||
      manifest$benchmark_scope != "five_replicate_fixed_setting_validation" ||
      manifest$replicate_count != 5L || manifest$expected_method_fit_count != 20L ||
      manifest$expected_chain_count != 80L || manifest$successful_chain_count != 80L ||
      nrow(status) != 80L || anyDuplicated(key) || any(status$status != "ok") ||
      !all(table(interaction(status$architecture, status$replicate, status$method,
        drop = TRUE)) == 4L) ||
      any(!is.finite(diagnostics$rhat)) || any(!is.finite(diagnostics$ess_bulk)) ||
      any(!is.finite(diagnostics$ess_tail)) || any(!is.finite(diagnostics$relative_mcse)) ||
      nrow(summary) != 4L || any(summary$replicate_count != 5L) ||
      any(!summary$recommendation_validation_status %in% c(
        "supported_in_all_replicates", "supported_in_most_replicates", "mixed",
        "not_supported", "indeterminate")))
    stop("Study 04 validation capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}

.five_promote <- function(study, source_dir, destination, source_files, title, validator) {
  if (dir.exists(destination)) {
    validator(destination)
    return(invisible(destination))
  }
  manifest_name <- if (study == "study02") "prediction_manifest.json" else "benchmark_manifest.json"
  output_files <- setdiff(.five_required_files(study),
    c("README.md", "benchmark_manifest.json", "example_data_manifest.csv",
      "source_files.csv", "session_info.txt", "checksums.csv"))
  missing <- output_files[!file.exists(file.path(source_dir, output_files))]
  if (length(missing) || !file.exists(file.path(source_dir, manifest_name)))
    stop("Promotion source outputs are incomplete: ", paste(missing, collapse = ", "), call. = FALSE)
  staging <- file.path("results", "local", "five_replicate_overnight", "promotion_staging",
    paste0(basename(destination), "-", Sys.getpid()))
  if (dir.exists(staging)) stop("Promotion staging directory already exists.", call. = FALSE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(file.path(source_dir, output_files), staging, overwrite = FALSE)
  if (!all(copied)) stop("Failed to copy promotion outputs.", call. = FALSE)
  file.copy(file.path(source_dir, manifest_name),
    file.path(staging, "benchmark_manifest.json"), overwrite = FALSE)
  old_manifest <- file.path("results", "reference", switch(study,
    study02 = "02_prediction/st-bayesc-bayesr-one-replicate-development-v1",
    study03 = "03_parameter_estimation/st-parameter-estimation-one-replicate-development-v1",
    study04 = "04_convergence/st-multichain-convergence-development-v1"),
    "example_data_manifest.csv")
  file.copy(old_manifest, file.path(staging, "example_data_manifest.csv"), overwrite = FALSE)
  src <- data.frame(file = source_files,
    md5 = unname(.study02_canonical_md5(source_files)), stringsAsFactors = FALSE)
  write.csv(src, file.path(staging, "source_files.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(staging, "session_info.txt"))
  writeLines(c(paste0("# ", title), "", "Complete five-replicate development evidence.",
    "", "This capsule contains compact tabular outputs and provenance only; native fits and target objects are excluded.",
    "The corresponding one-replicate/development capsule remains an immutable historical reference."),
    file.path(staging, "README.md"))
  checks <- .five_capsule_checksums(staging)
  write.csv(checks, file.path(staging, "checksums.csv"), row.names = FALSE)
  validator(staging)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging, destination)) stop("Atomic capsule promotion rename failed.", call. = FALSE)
  validator(destination)
  invisible(destination)
}

.five_promote_study02 <- function(source_dir) .five_promote("study02", source_dir,
  file.path("results", "reference", "02_prediction",
    "st-bayesc-bayesr-five-replicate-development-v1"),
  c("studies/five_replicate_helpers.R", "studies/five_replicate_promotion.R",
    "studies/02_prediction/config.R", "studies/02_prediction/pilot.R",
    "studies/02_prediction/targets.R", "studies/02_prediction/promotion.R"),
  "Single-trait prediction five-replicate development benchmark",
  .five_validate_study02_capsule)

.five_promote_study03 <- function(source_dir) .five_promote("study03", source_dir,
  file.path("results", "reference", "03_parameter_estimation",
    "st-parameter-estimation-five-replicate-development-v1"),
  c("studies/five_replicate_helpers.R", "studies/five_replicate_promotion.R",
    "studies/03_parameter_estimation/config.R", "studies/03_parameter_estimation/simulation.R",
    "studies/03_parameter_estimation/methods.R", "studies/03_parameter_estimation/metrics.R",
    "studies/03_parameter_estimation/pilot.R", "studies/03_parameter_estimation/targets.R"),
  "Single-trait parameter estimation five-replicate development benchmark",
  .five_validate_study03_capsule)

.five_promote_study04 <- function(source_dir) .five_promote("study04", source_dir,
  file.path("results", "reference", "04_convergence",
    "st-multichain-convergence-validation-five-replicate-v1"),
  c("studies/five_replicate_helpers.R", "studies/five_replicate_promotion.R",
    "studies/04_convergence/config.R", "studies/04_convergence/methods.R",
    "studies/04_convergence/chain_extraction.R", "studies/04_convergence/diagnostics.R",
    "studies/04_convergence/pilot.R", "studies/04_convergence/validation_targets.R"),
  "Five-replicate validation of recommended MCMC settings",
  .five_validate_study04_capsule)
