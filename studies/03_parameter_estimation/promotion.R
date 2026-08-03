.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-provenance.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-capsules.R"), local = TRUE)
source(file.path(.sblrbench_root, "R", "benchmark-validation.R"), local = TRUE)
.study03_canonical_md5 <- benchmark_canonical_md5

.study03_required_files <- function() c("README.md", "benchmark_manifest.json",
  "estimand_registry.csv", "simulation_truth.csv", "parameter_estimates.csv",
  "parameter_recovery_summary.csv", "paired_parameter_differences.csv",
  "paired_comparison_summary.csv", "computational_summary.csv", "replicate_status.csv",
  "seed_registry.csv", "example_data_manifest.csv", "source_files.csv", "session_info.txt",
  "checksums.csv")

.study03_validate_capsule <- function(path) {
  required <- .study03_required_files()
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Missing capsule files: ", paste(missing, collapse = ", "), call. = FALSE)
  checks <- utils::read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  if (anyDuplicated(checks$file) || any(grepl("(^|[\\/])[.][.]([\\/]|$)|^[A-Za-z]:|^/", checks$file)))
    stop("Invalid checksum paths.", call. = FALSE)
  if (!all(setdiff(required, "checksums.csv") %in% checks$file) ||
      any(!grepl("^[0-9a-f]{32}$", checks$md5)))
    stop("Incomplete or malformed checksum inventory.", call. = FALSE)
  observed <- .study03_canonical_md5(file.path(path, checks$file))
  bad <- checks$file[unname(observed) != checks$md5]
  if (length(bad)) stop("Capsule checksum validation failed: ", paste(bad, collapse = ", "), call. = FALSE)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- utils::read.csv(file.path(path, "replicate_status.csv"), stringsAsFactors = FALSE)
  estimates <- utils::read.csv(file.path(path, "parameter_estimates.csv"), stringsAsFactors = FALSE)
  available <- estimates$status == "ok"
  keys <- paste(status$architecture, status$replicate, status$method, sep = "::")
  if (!identical(sort(unique(status$method)), sort(c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"))) ||
      nrow(status) != 40L || anyDuplicated(keys) || any(status$status != "ok") ||
      manifest$replicate_count != 5L || manifest$successful_fit_count != 40L ||
      any(!is.finite(estimates$posterior_mean[available])) ||
      any(!is.finite(estimates$truth))) stop("Study 03 capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}
