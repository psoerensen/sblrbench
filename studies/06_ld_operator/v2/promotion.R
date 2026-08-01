.study06v2_promotion_root <- if (file.exists(file.path("studies",
  "02_prediction", "promotion.R"))) "." else file.path("..", "..", "..")
source(file.path(.study06v2_promotion_root, "studies", "02_prediction",
  "promotion.R"), local = TRUE)

.study06v2_capsule_required <- function(type = c("convergence", "benchmark")) {
  type <- match.arg(type)
  common <- c("README.md", "benchmark_manifest.json", "config.R",
    "example_data_manifest.csv", "source_files.csv", "session_info.txt",
    "reproduce.R", "contract_smoke_test.R", "checksums.csv",
    "v1_to_v2_design_crosswalk.csv", "deterministic_identity_summary.csv",
    "low_rank_block_diagnostics.csv", "deterministic_validation_manifest.json")
  if (type == "convergence") c(common,
    "historical_v1_inventory.md", "repin_hash_comparison.csv",
    "one_replicate_operator_pilot.csv", "pilot_gate.csv", "fit_status.csv",
    "candidate_settings.csv", "convergence_diagnostics.csv",
    "method_recommendations.csv")
  else c(common, "block_definitions.csv", "seed_registry.csv",
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
  files <- sort(setdiff(list.files(path, recursive = FALSE), "checksums.csv"))
  info <- file.info(file.path(path, files))
  data.frame(file = files, size_bytes = info$size,
    md5 = unname(.study02_canonical_md5(file.path(path, files))),
    stringsAsFactors = FALSE)
}

.study06v2_validate_checksums <- function(path, required) {
  inventory <- read.csv(file.path(path, "checksums.csv"),
    stringsAsFactors = FALSE)
  if (!identical(names(inventory), c("file", "size_bytes", "md5")) ||
      anyNA(inventory$file) || anyDuplicated(inventory$file) ||
      any(inventory$file != basename(inventory$file)) ||
      any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])\\.\\.([/\\\\]|$))",
        inventory$file)) ||
      !setequal(inventory$file, setdiff(required, "checksums.csv")) ||
      any(!grepl("^[0-9a-f]{32}$", inventory$md5)))
    stop("Invalid Study 06 v2 checksum inventory.", call. = FALSE)
  paths <- file.path(path, inventory$file)
  if (any(!file.exists(paths)) ||
      any(unname(.study02_canonical_md5(paths)) != inventory$md5))
    stop("Study 06 v2 canonical checksum validation failed.", call. = FALSE)
  invisible(TRUE)
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
  if (type == "benchmark") {
    paired <- read.csv(file.path(path, "paired_comparison_summary.csv"),
      stringsAsFactors = FALSE)
    low <- read.csv(file.path(path, "low_rank_fit_block_diagnostics.csv"),
      stringsAsFactors = FALSE)
    if (manifest$replicate_count != 5L ||
        any(paired$complete_paired_replicates != 5L) ||
        !setequal(unique(low$configuration),
          c("low_rank_full", "low_rank_0999", "low_rank_0995")) ||
        any(low$representation != "low_rank") ||
        !readLines(file.path(path, "readiness_decision.txt"), warn = FALSE)[1L] %in%
          c("READY FOR MTBLR LOW-RANK EXTENSION",
            "NOT READY FOR MTBLR LOW-RANK EXTENSION"))
      stop("Study 06 v2 benchmark capsule semantic validation failed.",
        call. = FALSE)
  }
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
    source_status = "pinned_local_source_snapshot_offline",
    validation_status = "complete_grid_validated",
    failures = unname(status$error_message[status$status != "ok"]))
}

.study06v2_copy_unique <- function(name, source_dirs, destination) {
  candidates <- file.path(source_dirs, name)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) != 1L ||
      !file.copy(existing, file.path(destination, name), overwrite = FALSE))
    stop("Study 06 v2 promotion source missing or ambiguous: ", name,
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
      file.path(local, "deterministic"), file.path(local, "aggregate"),
      file.path(local, "convergence"))
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
  for (name in generated) .study06v2_copy_unique(name, source_dirs, staging)
  file.copy(file.path("studies", "06_ld_operator", "v2", "config.R"),
    file.path(staging, "config.R"))
  file.copy("scripts/run_study06_low_rank_operator_v2.R",
    file.path(staging, "reproduce.R"))
  file.copy(file.path("studies", "06_ld_operator", "v2",
    "contract_smoke_test.R"), file.path(staging, "contract_smoke_test.R"))
  file.copy(file.path("results", "reference", "06_ld_operator",
    "st-ld-operator-five-replicate-development-v1",
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
    "studies/06_ld_operator/retained-low-rank-operator-development-v2.qmd")
  sources <- sort(sources[file.exists(sources)])
  write.csv(data.frame(file = sources,
    md5 = unname(.study02_canonical_md5(sources))),
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
