# Reusable structural checks for compact reference capsules.

benchmark_validate_capsule_checksums <- function(path, required,
                                                 allow_extra = FALSE) {
  if (!dir.exists(path)) stop("Capsule directory does not exist: ", path,
    call. = FALSE)
  missing <- setdiff(required, list.files(path, recursive = FALSE))
  if (length(missing))
    stop("Capsule is missing required files: ", paste(missing,
      collapse = ", "), call. = FALSE)
  inventory <- utils::read.csv(file.path(path, "checksums.csv"),
    stringsAsFactors = FALSE)
  expected_names <- c("file", "size_bytes", "md5")
  expected_files <- setdiff(required, "checksums.csv")
  invalid <- !identical(names(inventory), expected_names) ||
    anyNA(inventory$file) || anyDuplicated(inventory$file) ||
    any(inventory$file != basename(inventory$file)) ||
    any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])[.][.]([/\\\\]|$))",
      inventory$file)) || any(!grepl("^[0-9a-f]{32}$", inventory$md5))
  if (!allow_extra) invalid <- invalid ||
    !setequal(inventory$file, expected_files)
  if (invalid) stop("Invalid capsule checksum inventory.", call. = FALSE)
  paths <- file.path(path, inventory$file)
  if (any(!file.exists(paths)) ||
      any(unname(benchmark_canonical_md5(paths)) != inventory$md5))
    stop("Capsule checksum validation failed.", call. = FALSE)
  invisible(TRUE)
}

validate_prediction_capsule <- function(path, spec) {
  validate_benchmark_spec(spec)
  required <- c("README.md", "benchmark_manifest.json", "benchmark_summary.csv",
    "paired_comparison_summary.csv", "prediction_metrics.csv",
    "paired_method_differences.csv", "computational_summary.csv",
    "replicate_status.csv", "simulation_summary.csv", "seed_registry.csv",
    "example_data_manifest.csv", "source_files.csv", "session_info.txt",
    "checksums.csv")
  benchmark_validate_capsule_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  metrics <- utils::read.csv(file.path(path, "prediction_metrics.csv"),
    stringsAsFactors = FALSE)
  status <- utils::read.csv(file.path(path, "replicate_status.csv"),
    stringsAsFactors = FALSE)
  simulations <- utils::read.csv(file.path(path, "simulation_summary.csv"),
    stringsAsFactors = FALSE)
  expected_methods <- names(spec$methods)
  expected_scenarios <- names(spec$scenarios)
  expected <- expand.grid(architecture = expected_scenarios,
    replicate = seq_len(spec$replicate_count), method = expected_methods,
    stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate, x$method, sep = "::")
  if (nrow(status) != nrow(expected) || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || any(status$status != "ok"))
    stop("Prediction capsule does not contain the exact complete fit grid.",
      call. = FALSE)
  if (!identical(unname(manifest$active_methods), expected_methods) ||
      !identical(unname(manifest$architectures), expected_scenarios) ||
      !identical(as.integer(manifest$replicate_count), spec$replicate_count) ||
      !identical(as.integer(manifest$successful_fit_count), nrow(expected)) ||
      !identical(as.integer(manifest$failed_fit_count), 0L))
    stop("Prediction capsule manifest disagrees with the Study 02 specification.",
      call. = FALSE)
  required_metrics <- c("prediction_correlation", "prediction_mse",
    "prediction_nmse", "phenotype_prediction_correlation",
    "prediction_calibration_intercept", "prediction_calibration_slope",
    "effect_rmse")
  observed <- unique(metrics[c("architecture", "replicate", "method", "metric")])
  if (nrow(observed) != nrow(expected) * length(required_metrics) ||
      !setequal(metrics$metric, required_metrics) || any(metrics$status != "ok") ||
      any(!is.finite(metrics$value)))
    stop("Prediction capsule metric grid is incomplete or invalid.",
      call. = FALSE)
  simulation_keys <- paste(simulations$architecture, simulations$replicate,
    sep = "::")
  if (nrow(simulations) != length(expected_scenarios) * spec$replicate_count ||
      anyDuplicated(simulation_keys) ||
      any(simulations$causal_count != spec$controls$simulation$n_causal) ||
      any(abs(simulations$realized_h2 - spec$controls$simulation$h2) >
        spec$validation$oracle_tolerance) || any(!simulations$oracle_ok))
    stop("Prediction capsule simulation validation is incomplete.",
      call. = FALSE)
  invisible(TRUE)
}
