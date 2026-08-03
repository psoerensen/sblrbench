.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-provenance.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-capsules.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-validation.R"), local = TRUE)

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
      metric = x$metric[[1]], value = if (length(values) == 1L) values else if (length(values)) mean(values) else NA_real_,
      replicate_count = length(unique(x$replicate)), successful_replicates = length(values),
      mean = if (length(values)) mean(values) else NA_real_,
      sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median = if (length(values)) stats::median(values) else NA_real_,
      minimum = if (length(values)) min(values) else NA_real_,
      maximum = if (length(values)) max(values) else NA_real_)
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
      replicate_count = length(unique(x$replicate)), successful_replicates = length(values),
      mean_advantage = if (length(values)) mean(values) else NA_real_,
      sd_advantage = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median_advantage = if (length(values)) stats::median(values) else NA_real_,
      minimum_advantage = if (length(values)) min(values) else NA_real_,
      maximum_advantage = if (length(values)) max(values) else NA_real_,
      complete_pairs = length(values))
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

# Compatibility alias for callers not yet migrated from Study 02 promotion.
# Remove after Studies 01 and 03--07 call benchmark_canonical_md5() directly.
.study02_canonical_md5 <- benchmark_canonical_md5

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
  if (!identical(as.integer(manifest$replicate_count), 5L))
    stop("The current capsule requires exactly five replicates.", call. = FALSE)
  expected <- expand.grid(architecture = architectures, replicate = 1:5,
    method = methods, stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate, x$method, sep = "|")
  if (nrow(status) != 40L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || any(status$status != "ok"))
    stop("Promotion requires exactly 40 complete successful fits.", call. = FALSE)
  if (nrow(computation) != 40L || anyDuplicated(key(computation)) ||
      !setequal(key(computation), key(expected)) || any(computation$status != "ok") ||
      any(!is.finite(computation$runtime))) stop("Computational summaries are invalid.", call. = FALSE)
  observed <- unique(metrics[, c("architecture", "replicate", "method", "metric")])
  if (nrow(observed) != 280L || !setequal(unique(metrics$metric), .study02_required_metrics()) ||
      any(metrics$status != "ok") || any(!is.finite(metrics$value)))
    stop("Required prediction metrics are incomplete or invalid.", call. = FALSE)
  simulation_key <- paste(simulations$architecture, simulations$replicate, sep = "::")
  if (nrow(simulations) != 10L || anyDuplicated(simulation_key) ||
      !setequal(simulations$architecture, architectures) ||
      any(simulations$causal_count != 50L) || any(!simulations$oracle_ok) ||
      any(abs(simulations$realized_h2 - 0.30) > 1e-10))
    stop("Simulation validation is incomplete.", call. = FALSE)
  if (!setequal(unique(paired$comparison_id), .study02_required_comparisons()) ||
      nrow(paired) != 280L || any(!paired$complete_pair) ||
      any(!is.finite(paired$advantage))) stop("Paired comparisons are incomplete.", call. = FALSE)
  if (!identical(manifest$task, "single_trait_prediction") ||
      !identical(manifest$benchmark_scope, "current") ||
      !identical(manifest$benchmark_status, "complete") ||
      !identical(unname(manifest$active_methods), methods) ||
      manifest$training_sample_count != 3500L || manifest$test_sample_count != 1500L ||
      manifest$canonical_marker_count != 37991L || manifest$expected_fit_count != 40L ||
      manifest$successful_fit_count != 40L || manifest$failed_fit_count != 0L)
    stop("Manifest does not describe the completed benchmark.", call. = FALSE)
  if (is.null(manifest$qgdata$commit) || is.null(manifest$qgdata$files))
    stop("qgdata provenance is missing.", call. = FALSE)
  if (!is.null(benchmark_summary)) {
    expected_summary <- .study02_benchmark_summary(metrics)
    if (!isTRUE(all.equal(benchmark_summary, expected_summary,
        check.attributes = FALSE)))
      stop("Benchmark summary disagrees with metric rows.", call. = FALSE)
  }
  if (!is.null(paired_summary)) {
    expected_paired <- .study02_paired_summary(paired)
    if (!isTRUE(all.equal(paired_summary, expected_paired,
        check.attributes = FALSE)))
      stop("Paired summary disagrees with paired rows.", call. = FALSE)
  }
  invisible(TRUE)
}

.study02_validate_capsule <- function(path) {
  required <- c("README.md", "benchmark_manifest.json", "benchmark_summary.csv",
    "paired_comparison_summary.csv", "prediction_metrics.csv",
    "paired_method_differences.csv", "computational_summary.csv",
    "replicate_status.csv", "simulation_summary.csv", "seed_registry.csv",
    "example_data_manifest.csv", "source_files.csv", "session_info.txt",
    "checksums.csv")
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Incomplete prediction capsule: ", paste(missing, collapse = ", "), call. = FALSE)
  checks <- read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  if (!identical(names(checks), c("file", "size_bytes", "md5")))
    stop("checksums.csv has an invalid schema.", call. = FALSE)
  if (anyNA(checks$file) || anyDuplicated(checks$file))
    stop("checksums.csv contains duplicate or missing filenames.", call. = FALSE)
  if (any(checks$file != basename(checks$file)) ||
      any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])\\.\\.([/\\\\]|$))", checks$file)))
    stop("checksums.csv contains a path outside the capsule directory.", call. = FALSE)
  required_listed <- setdiff(required, "checksums.csv")
  absent_checksums <- setdiff(required_listed, checks$file)
  if (length(absent_checksums))
    stop("Required capsule files are absent from checksums.csv: ",
      paste(absent_checksums, collapse = ", "), call. = FALSE)
  listed_paths <- file.path(path, checks$file)
  missing_listed <- checks$file[!file.exists(listed_paths)]
  if (length(missing_listed))
    stop("Files listed in checksums.csv are missing: ",
      paste(missing_listed, collapse = ", "), call. = FALSE)
  if (anyNA(checks$md5) || any(!grepl("^[0-9a-f]{32}$", checks$md5)))
    stop("checksums.csv contains malformed MD5 values.", call. = FALSE)
  actual <- unname(.study02_canonical_md5(listed_paths))
  mismatched <- which(actual != checks$md5)
  if (length(mismatched)) {
    details <- paste0(checks$file[mismatched], ": expected ",
      checks$md5[mismatched], ", observed ", actual[mismatched])
    stop("Capsule checksum validation failed:\n", paste(details, collapse = "\n"),
      call. = FALSE)
  }
  metrics <- read.csv(file.path(path, "prediction_metrics.csv"))
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  methods <- .study02_expected_methods()
  status <- read.csv(file.path(path, "replicate_status.csv"),
    stringsAsFactors = FALSE)
  keys <- paste(status$architecture, status$replicate, status$method, sep = "::")
  if (manifest$replicate_count != 5L || manifest$successful_fit_count != 40L ||
      nrow(status) != 40L || anyDuplicated(keys) || any(status$status != "ok") ||
      !identical(unname(manifest$active_methods), methods) ||
      !setequal(unique(metrics$method), methods) || any(grepl("^mt_", metrics$method)) ||
      any(!is.finite(metrics$value)))
    stop("Capsule report inputs are invalid.", call. = FALSE)
  invisible(TRUE)
}
