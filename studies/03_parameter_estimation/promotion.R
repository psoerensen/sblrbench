study02_promotion <- c(file.path("studies", "02_prediction", "promotion.R"),
  file.path("..", "..", "studies", "02_prediction", "promotion.R"))
study02_promotion <- study02_promotion[file.exists(study02_promotion)][1L]
if (is.na(study02_promotion)) stop("Cannot locate shared canonical checksum helper.", call. = FALSE)
source(study02_promotion, local = TRUE)
.study03_canonical_md5 <- .study02_canonical_md5

.study03_required_files <- function() c("README.md", "benchmark_manifest.json",
  "estimand_registry.csv", "simulation_truth.csv", "parameter_estimates.csv",
  "parameter_recovery_summary.csv", "paired_parameter_differences.csv",
  "paired_comparison_summary.csv", "computational_summary.csv", "replicate_status.csv",
  "seed_registry.csv", "example_data_manifest.csv", "source_files.csv", "session_info.txt",
  "config.R", "estimands.R", "simulation.R", "methods.R", "metrics.R",
  "run_parameter_estimation_benchmark.R", "parameter_estimation_contract_smoke_test.R",
  "worked_parameter_estimation_example.R", "pilot.R", "targets.R")

.study03_validate_capsule <- function(path) {
  required <- .study03_required_files()
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Missing capsule files: ", paste(missing, collapse = ", "), call. = FALSE)
  checks <- utils::read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  if (anyDuplicated(checks$file) || any(grepl("(^|[\\/])[.][.]([\\/]|$)|^[A-Za-z]:|^/", checks$file)))
    stop("Invalid checksum paths.", call. = FALSE)
  if (!all(required %in% checks$file) || any(!grepl("^[0-9a-f]{32}$", checks$md5)))
    stop("Incomplete or malformed checksum inventory.", call. = FALSE)
  observed <- .study03_canonical_md5(file.path(path, checks$file))
  bad <- checks$file[unname(observed) != checks$md5]
  if (length(bad)) stop("Capsule checksum validation failed: ", paste(bad, collapse = ", "), call. = FALSE)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- utils::read.csv(file.path(path, "replicate_status.csv"), stringsAsFactors = FALSE)
  estimates <- utils::read.csv(file.path(path, "parameter_estimates.csv"), stringsAsFactors = FALSE)
  if (!identical(sort(unique(status$method)), sort(c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"))) ||
      nrow(status) != 8L || any(status$status != "ok") || manifest$replicate_count != 1L ||
      manifest$successful_fit_count != 8L || any(!is.finite(estimates$posterior_mean)) ||
      any(!is.finite(estimates$truth))) stop("Study 03 capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}
