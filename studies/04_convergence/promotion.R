study02_promotion <- c(file.path("studies", "02_prediction", "promotion.R"),
  file.path("..", "..", "studies", "02_prediction", "promotion.R"))
study02_promotion <- study02_promotion[file.exists(study02_promotion)][1L]
if (is.na(study02_promotion)) stop("Cannot locate canonical checksum helper.", call. = FALSE)
source(study02_promotion, local = TRUE)
.study04_canonical_md5 <- .study02_canonical_md5

.study04_required_files <- function() c(
  "README.md", "benchmark_manifest.json", "diagnostic_registry.csv",
  "diagnostic_thresholds.csv", "scalar_chain_draws.csv", "convergence_diagnostics.csv",
  "checkpoint_summary.csv", "burnin_stability.csv", "chain_summaries.csv",
  "chain_agreement.csv", "method_recommendations.csv", "computational_summary.csv",
  "chain_status.csv", "seed_registry.csv", "simulation_truth.csv",
  "example_data_manifest.csv", "source_files.csv", "session_info.txt", "checksums.csv",
  "config.R", "diagnostic_registry.R", "methods.R", "chain_extraction.R",
  "diagnostics.R", "recommendations.R", "run_convergence_benchmark.R",
  "convergence_contract_smoke_test.R", "worked_convergence_example.R", "pilot.R", "targets.R")

.study04_validate_capsule <- function(path) {
  required <- .study04_required_files()
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing)) stop("Missing capsule files: ", paste(missing, collapse = ", "), call. = FALSE)
  checks <- utils::read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  if (anyDuplicated(checks$file) || any(grepl("(^|[\\/])[.][.]([\\/]|$)|^[A-Za-z]:|^/", checks$file)))
    stop("Invalid checksum paths.", call. = FALSE)
  if (!all(setdiff(required, "checksums.csv") %in% checks$file) ||
      any(!grepl("^[0-9a-f]{32}$", checks$md5)))
    stop("Incomplete or malformed checksum inventory.", call. = FALSE)
  bad <- checks$file[unname(.study04_canonical_md5(file.path(path, checks$file))) != checks$md5]
  if (length(bad)) stop("Capsule checksum validation failed: ", paste(bad, collapse = ", "), call. = FALSE)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- utils::read.csv(file.path(path, "chain_status.csv"), stringsAsFactors = FALSE)
  draws <- utils::read.csv(file.path(path, "scalar_chain_draws.csv"), stringsAsFactors = FALSE)
  diagnostics <- utils::read.csv(file.path(path, "convergence_diagnostics.csv"), stringsAsFactors = FALSE)
  rec <- utils::read.csv(file.path(path, "method_recommendations.csv"), stringsAsFactors = FALSE)
  expected <- c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr")
  numeric_draws <- as.matrix(draws[c("effect_variance", "genetic_variance", "residual_variance", "heritability")])
  if (!identical(sort(unique(status$method)), sort(expected)) || nrow(status) != 16L ||
      anyDuplicated(status[c("architecture", "method", "chain")]) || any(status$status != "ok") ||
      !all(table(interaction(status$architecture, status$method, drop = TRUE)) == 4L) ||
      anyDuplicated(draws[c("architecture", "method", "chain", "raw_iteration")]) ||
      any(!is.finite(numeric_draws)) || any(numeric_draws[, 1:3] < 0) ||
      any(draws$heritability < 0 | draws$heritability > 1) || nrow(rec) != 4L ||
      any(!diagnostics$status %in% c("pass", "fail", "indeterminate", "unavailable")) ||
      manifest$expected_chain_count != 16L || manifest$successful_chain_count != 16L)
    stop("Study 04 capsule semantic validation failed.", call. = FALSE)
  invisible(TRUE)
}
