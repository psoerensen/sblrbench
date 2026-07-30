.study02_expected_methods <- function() c(
  "st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"
)

.study02_required_metrics <- function() c(
  "prediction_correlation", "prediction_mse", "prediction_nmse",
  "phenotype_prediction_correlation", "prediction_calibration_intercept",
  "prediction_calibration_slope", "effect_rmse"
)

.study02_required_comparisons <- function() c(
  "bayesr_vs_bayesc_bed", "sbayesr_vs_sbayesc_csr",
  "csr_vs_bed_bayesc", "csr_vs_bed_bayesr"
)

.study02_method_labels <- function() c(
  st_bed_bayesc = "ST-BED BayesC", st_bed_bayesr = "ST-BED BayesR",
  st_csr_sbayesc = "ST-CSR SBayesC", st_csr_sbayesr = "ST-CSR SBayesR"
)

.study02_comparison_labels <- function() c(
  bayesr_vs_bayesc_bed = "BED BayesR versus BayesC",
  sbayesr_vs_sbayesc_csr = "CSR SBayesR versus SBayesC",
  csr_vs_bed_bayesc = "CSR versus BED BayesC class",
  csr_vs_bed_bayesr = "CSR versus BED BayesR class"
)

.study02_benchmark_summary <- function(metrics) {
  key <- interaction(metrics$architecture, metrics$method, metrics$metric,
    drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(metrics, key), function(x) {
    values <- x$value[x$status == "ok"]
    data.frame(architecture = x$architecture[[1]], method = x$method[[1]],
      method_label = unname(.study02_method_labels()[x$method[[1]]]),
      metric = x$metric[[1]], value = if (length(values) == 1L) values else mean(values),
      replicate_count = length(values), mean = mean(values),
      sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median = stats::median(values), minimum = min(values), maximum = max(values))
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

.study02_paired_summary <- function(paired) {
  key <- interaction(paired$architecture, paired$comparison_id,
    paired$paired_metric, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(paired, key), function(x) {
    values <- x$advantage[x$complete_pair]
    data.frame(architecture = x$architecture[[1]],
      comparison_id = x$comparison_id[[1]],
      comparison_label = unname(.study02_comparison_labels()[x$comparison_id[[1]]]),
      focal_method = x$focal_method[[1]], comparison_method = x$comparison_method[[1]],
      paired_metric = x$paired_metric[[1]], advantage = if (length(values) == 1L) values else mean(values),
      replicate_count = length(unique(x$replicate)), mean_advantage = mean(values),
      sd_advantage = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median_advantage = stats::median(values), complete_pairs = length(values))
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

.study02_validate_promotion_tables <- function(config, metrics, paired,
                                                computation, status,
                                                simulations, manifest,
                                                benchmark_summary = NULL,
                                                paired_summary = NULL) {
  methods <- .study02_expected_methods()
  architectures <- c("sparse_homogeneous", "sparse_mixture")
  if (!identical(config$methods, methods)) stop("Invalid active method set.", call. = FALSE)
  all_methods <- unique(c(metrics$method, computation$method, status$method))
  if (!setequal(all_methods, methods) || any(grepl("^mt_", all_methods)))
    stop("Promotion refuses MT rows or methods outside the exact four-method set.", call. = FALSE)
  if (!identical(as.integer(manifest$replicate_count), 1L))
    stop("This capsule requires exactly one replicate.", call. = FALSE)
  expected <- expand.grid(architecture = architectures, replicate = 1L,
    method = methods, stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate, x$method, sep = "|")
  if (nrow(status) != 8L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || any(status$status != "ok"))
    stop("Promotion requires exactly eight complete successful fits.", call. = FALSE)
  if (nrow(computation) != 8L || anyDuplicated(key(computation)) ||
      !setequal(key(computation), key(expected)) || any(computation$status != "ok") ||
      any(!is.finite(computation$runtime))) stop("Computational summaries are invalid.", call. = FALSE)
  observed <- unique(metrics[, c("architecture", "replicate", "method", "metric")])
  if (nrow(observed) != 56L || !setequal(unique(metrics$metric), .study02_required_metrics()) ||
      any(metrics$status != "ok") || any(!is.finite(metrics$value)))
    stop("Required prediction metrics are incomplete or invalid.", call. = FALSE)
  if (nrow(simulations) != 2L || anyDuplicated(simulations$architecture) ||
      !setequal(simulations$architecture, architectures) ||
      any(simulations$causal_count != 50L) || any(!simulations$oracle_ok) ||
      any(abs(simulations$realized_h2 - 0.30) > 1e-10))
    stop("Simulation validation is incomplete.", call. = FALSE)
  if (!setequal(unique(paired$comparison_id), .study02_required_comparisons()) ||
      nrow(paired) != 56L || any(!paired$complete_pair) ||
      any(!is.finite(paired$advantage))) stop("Paired comparisons are incomplete.", call. = FALSE)
  if (!identical(manifest$task, "single_trait_prediction") ||
      !identical(manifest$benchmark_scope, "one_replicate_development") ||
      !identical(manifest$benchmark_status, "complete") ||
      !identical(unname(manifest$active_methods), methods) ||
      manifest$training_sample_count != 3500L || manifest$test_sample_count != 1500L ||
      manifest$canonical_marker_count != 37991L || manifest$expected_fit_count != 8L ||
      manifest$successful_fit_count != 8L || manifest$failed_fit_count != 0L)
    stop("Manifest does not describe the completed benchmark.", call. = FALSE)
  if (is.null(manifest$qgdata$commit) || is.null(manifest$qgdata$files))
    stop("qgdata provenance is missing.", call. = FALSE)
  if (!is.null(benchmark_summary)) {
    expected_summary <- .study02_benchmark_summary(metrics)
    comparable <- setdiff(names(expected_summary), "sd")
    if (!isTRUE(all.equal(benchmark_summary[comparable], expected_summary[comparable],
        check.attributes = FALSE)) || !all(is.na(benchmark_summary$sd)))
      stop("Benchmark summary disagrees with metric rows.", call. = FALSE)
  }
  if (!is.null(paired_summary)) {
    expected_paired <- .study02_paired_summary(paired)
    comparable <- setdiff(names(expected_paired), "sd_advantage")
    if (!isTRUE(all.equal(paired_summary[comparable], expected_paired[comparable],
        check.attributes = FALSE)) || !all(is.na(paired_summary$sd_advantage)))
      stop("Paired summary disagrees with paired rows.", call. = FALSE)
  }
  invisible(TRUE)
}

.study02_validate_capsule <- function(path) {
  required <- c("README.md", "benchmark_manifest.json", "benchmark_summary.csv",
    "paired_comparison_summary.csv", "prediction_metrics.csv",
    "paired_method_differences.csv", "computational_summary.csv",
    "replicate_status.csv", "simulation_summary.csv", "config.R",
    "run_prediction_benchmark.R", "worked_prediction_example.R",
    "prediction_contract_smoke_test.R", "example_data_manifest.csv",
    "source_files.csv", "session_info.txt", "checksums.csv", "pilot.R", "targets.R")
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Incomplete prediction capsule: ", paste(missing, collapse = ", "), call. = FALSE)
  checks <- read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  actual <- unname(tools::md5sum(file.path(path, checks$file)))
  if (anyNA(actual) || !identical(actual, checks$md5)) stop("Capsule checksum validation failed.", call. = FALSE)
  metrics <- read.csv(file.path(path, "prediction_metrics.csv"))
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  if (manifest$replicate_count != 1L || manifest$successful_fit_count != 8L ||
      any(grepl("^mt_", metrics$method)) || any(!is.finite(metrics$value)))
    stop("Capsule report inputs are invalid.", call. = FALSE)
  invisible(TRUE)
}
