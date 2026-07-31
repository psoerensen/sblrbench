.study06_promotion_root <- if (file.exists(file.path("studies",
  "02_prediction", "promotion.R"))) "." else file.path("..", "..")
source(file.path(.study06_promotion_root, "studies",
  "02_prediction", "promotion.R"), local = TRUE)

.study06_required <- function(type = c("convergence", "benchmark")) {
  type <- match.arg(type)
  common <- c("README.md", "benchmark_manifest.json", "config.R",
    "example_data_manifest.csv", "source_files.csv",
    "interface_audit_sources.csv", "session_info.txt",
    "reproduce.R", "contract_smoke_test.R", "checksums.csv")
  if (type == "convergence") c(common,
    "block_design_candidates.csv", "selected_block_design.csv",
    "selected_block_definitions.csv",
    "operator_equivalence_summary.csv",
    "block_level_operator_diagnostics.csv",
    "hard_filter_candidates.csv", "selected_hard_filter.csv",
    "synthetic_filter_validation.csv", "ridge_lw_diagnostics.csv",
    "ridge_lw_audit.csv", "fixed_ridge_candidates.csv",
    "selected_ridge_filter.csv", "one_replicate_operator_pilot.csv",
    "operator_perturbation_metrics.csv",
    "full_csr_operator_action_metrics.csv",
    "fit_status.csv",
    "convergence_diagnostics.csv", "candidate_settings.csv",
    "method_recommendations.csv", "seed_registry.csv",
    "computational_summary.csv")
  else c(common, "selected_block_definitions.csv",
    "selected_block_design.csv", "filter_specification.csv",
    "seed_registry.csv", "simulation_summary.csv",
    "operator_summaries.csv",
    "operator_perturbation_metrics.csv",
    "full_csr_operator_action_metrics.csv",
    "synthetic_filter_validation.csv",
    "block_level_operator_diagnostics.csv", "fit_status.csv",
    "prediction_metrics.csv", "parameter_estimates.csv",
    "parameter_recovery_summary.csv", "marker_effect_metrics.csv",
    "marker_effect_agreement.csv",
    "sbayesr_diagnostic_evidence.csv",
    "paired_replicate_differences.csv",
    "paired_comparison_summary.csv",
    "convergence_diagnostics.csv",
    "convergence_validation_summary.csv",
    "computational_summary.csv",
    "method_output_availability.csv")
}

.study06_checksums <- function(path) {
  files <- sort(setdiff(list.files(path, recursive = FALSE),
    "checksums.csv"))
  info <- file.info(file.path(path, files))
  data.frame(file = files, size_bytes = info$size,
    md5 = unname(.study02_canonical_md5(file.path(path, files))),
    stringsAsFactors = FALSE)
}

.study06_validate_checksums <- function(path, required) {
  x <- utils::read.csv(file.path(path, "checksums.csv"),
    stringsAsFactors = FALSE)
  if (!identical(names(x), c("file", "size_bytes", "md5")) ||
      anyNA(x$file) || anyDuplicated(x$file) ||
      any(x$file != basename(x$file)) ||
      any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])\\.\\.([/\\\\]|$))",
        x$file)) ||
      !setequal(x$file, setdiff(required, "checksums.csv")) ||
      any(!grepl("^[0-9a-f]{32}$", x$md5)))
    stop("Invalid Study 06 checksum inventory.", call. = FALSE)
  paths <- file.path(path, x$file)
  if (any(!file.exists(paths)) ||
      any(unname(.study02_canonical_md5(paths)) != x$md5))
    stop("Study 06 canonical checksum validation failed.",
      call. = FALSE)
  invisible(TRUE)
}

.study06_validate_convergence_capsule <- function(path) {
  required <- .study06_required("convergence")
  if (!dir.exists(path) ||
      length(setdiff(required, list.files(path))))
    stop("Study 06 convergence capsule is incomplete.", call. = FALSE)
  .study06_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path,
    "benchmark_manifest.json"), simplifyVector = TRUE)
  gate <- read.csv(file.path(path,
    "operator_equivalence_summary.csv"))
  blocks <- read.csv(file.path(path,
    "selected_block_definitions.csv"))
  filter <- read.csv(file.path(path, "selected_hard_filter.csv"))
  synthetic <- read.csv(file.path(path, "synthetic_filter_validation.csv"))
  ridge <- read.csv(file.path(path, "selected_ridge_filter.csv"))
  status <- read.csv(file.path(path, "fit_status.csv"))
  rec <- read.csv(file.path(path, "method_recommendations.csv"))
  candidates <- read.csv(file.path(path, "candidate_settings.csv"))
  expected <- expand.grid(
    architecture = c("sparse_homogeneous", "sparse_mixture"),
    configuration = c("block_csr", "block_eigen_unfiltered",
      "block_eigen_hard", "block_eigen_ridge_fixed"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$configuration,
    sep = "|")
  if (!identical(manifest$task,
      "single_trait_ld_operator_convergence_selection") ||
      !identical(manifest$benchmark_status, "complete") ||
      manifest$expected_fit_count != 8L ||
      nrow(gate) != 1L || !isTRUE(gate$pass[1L]) ||
      nrow(blocks) < 1L || blocks$start[1L] != 1L ||
      any(blocks$start[-1L] != blocks$end[-nrow(blocks)] + 1L) ||
      nrow(filter) != 1L || !isTRUE(filter$pass[1L]) ||
      !is.finite(filter$eigen_tau[1L]) || filter$eigen_tau[1L] <= 0 ||
      !filter$filter_activity_status[1L] %in% c("active", "effective_no_op") ||
      !nrow(synthetic) || any(!synthetic$pass) ||
      nrow(ridge) != 1L || !isTRUE(ridge$pass[1L]) ||
      !is.finite(ridge$eigen_eta[1L]) || ridge$eigen_eta[1L] <= 0 ||
      nrow(status) != 8L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) ||
      any(status$status != "ok") || any(status$chain_count != 4L) ||
      nrow(rec) != 8L || anyDuplicated(key(rec)) ||
      !setequal(key(rec), key(expected)) ||
      any(rec$recommendation_status != "available") ||
      any(rec$nchains != 4L) || any(rec$ncores != 4L) ||
      any(rec$nthin != 1L) ||
      !all(vapply(split(candidates$pass,
        key(candidates)), any, logical(1))))
    stop("Study 06 convergence capsule semantic validation failed.",
      call. = FALSE)
  invisible(TRUE)
}

.study06_validate_benchmark_capsule <- function(path) {
  required <- .study06_required("benchmark")
  if (!dir.exists(path) ||
      length(setdiff(required, list.files(path))))
    stop("Study 06 benchmark capsule is incomplete.", call. = FALSE)
  .study06_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path,
    "benchmark_manifest.json"), simplifyVector = TRUE)
  status <- read.csv(file.path(path, "fit_status.csv"))
  prediction <- read.csv(file.path(path, "prediction_metrics.csv"))
  diagnostics <- read.csv(file.path(path,
    "convergence_diagnostics.csv"))
  paired <- read.csv(file.path(path,
    "paired_replicate_differences.csv"))
  blocks <- read.csv(file.path(path,
    "selected_block_definitions.csv"))
  filters <- read.csv(file.path(path, "filter_specification.csv"))
  expected <- expand.grid(
    architecture = c("sparse_homogeneous", "sparse_mixture"),
    replicate = 1:5,
    configuration = c("bed", "full_csr", "block_csr",
      "block_eigen_unfiltered", "block_eigen_hard",
      "block_eigen_ridge_fixed"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate,
    x$configuration, sep = "|")
  if (!identical(manifest$task,
      "single_trait_ld_operator_benchmark") ||
      !identical(manifest$benchmark_status, "complete") ||
      manifest$replicate_count != 5L ||
      manifest$expected_fit_count != 60L ||
      manifest$successful_fit_count != 60L ||
      manifest$failed_fit_count != 0L ||
      manifest$chains_per_fit != 4L ||
      nrow(status) != 60L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) ||
      any(status$status != "ok") || any(status$chain_count != 4L) ||
      length(unique(key(prediction))) != 60L ||
      any(!is.finite(prediction$value)) ||
      length(unique(key(diagnostics))) != 60L ||
      any(diagnostics$chain_count != 4L) ||
      any(!is.finite(diagnostics$rhat)) ||
      any(!is.finite(diagnostics$ess_bulk)) ||
      any(!is.finite(diagnostics$ess_tail)) ||
      any(!is.finite(diagnostics$relative_mcse)) ||
      any(table(interaction(paired$comparison_id,
        paired$architecture, paired$metric, drop = TRUE)) != 5L) ||
      any(!paired$complete_pair) ||
      blocks$start[1L] != 1L ||
      any(blocks$start[-1L] != blocks$end[-nrow(blocks)] + 1L) ||
      !setequal(filters$configuration,
        c("block_eigen_unfiltered", "block_eigen_hard",
          "block_eigen_ridge_fixed")) ||
      !file.exists(file.path(path, "source_files.csv")) ||
      !file.exists(file.path(path, "session_info.txt")))
    stop("Study 06 benchmark capsule semantic validation failed.",
      call. = FALSE)
  invisible(TRUE)
}

.study06_manifest <- function(type, source_dir, config) {
  status <- read.csv(file.path(source_dir, "fit_status.csv"),
    stringsAsFactors = FALSE)
  if (type == "convergence") list(
    study = "06_ld_operator",
    task = "single_trait_ld_operator_convergence_selection",
    benchmark_scope = "block_operator_maximum_history_development",
    benchmark_status = if (nrow(status) == 8L &&
      all(status$status == "ok")) "complete" else "incomplete",
    expected_fit_count = 8L, successful_fit_count =
      sum(status$status == "ok"),
    failed_fit_count = sum(status$status != "ok"),
    architectures = config$architectures,
    configurations = c("block_csr", "block_eigen_unfiltered",
      "block_eigen_hard", "block_eigen_ridge_fixed"),
    chains_per_fit = 4L, expected_chain_count = 32L)
  else list(
    study = "06_ld_operator",
    task = "single_trait_ld_operator_benchmark",
    benchmark_scope = "five_replicate_development",
    benchmark_status = if (nrow(status) == 60L &&
      all(status$status == "ok")) "complete" else "incomplete",
    architectures = config$architectures,
    configurations = config$configurations,
    replicate_count = 5L, expected_fit_count = 60L,
    successful_fit_count = sum(status$status == "ok"),
    failed_fit_count = sum(status$status != "ok"),
    chains_per_fit = 4L, expected_chain_count = 240L)
}

.study06_promote <- function(type = c("convergence", "benchmark"),
                             source_dirs, destination, config) {
  type <- match.arg(type)
  validator <- if (type == "convergence")
    .study06_validate_convergence_capsule else
      .study06_validate_benchmark_capsule
  if (dir.exists(destination)) {
    validator(destination)
    return(invisible(destination))
  }
  required <- .study06_required(type)
  status_source <- file.path(source_dirs, "fit_status.csv")
  status_source <- status_source[file.exists(status_source)]
  if (length(status_source) != 1L)
    stop("Study 06 promotion requires one fit-status source.",
      call. = FALSE)
  status <- read.csv(status_source, stringsAsFactors = FALSE)
  staging <- file.path(config$local_dir, "promotion_staging",
    paste0(basename(destination), "-", Sys.getpid()))
  if (dir.exists(staging))
    stop("Study 06 promotion staging already exists.", call. = FALSE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  generated <- setdiff(required, c("README.md",
    "benchmark_manifest.json", "config.R",
    "example_data_manifest.csv", "source_files.csv",
    "interface_audit_sources.csv", "session_info.txt",
    "reproduce.R", "contract_smoke_test.R", "checksums.csv"))
  for (name in generated) {
    candidates <- file.path(source_dirs, name)
    existing <- candidates[file.exists(candidates)]
    if (length(existing) != 1L ||
        !file.copy(existing, file.path(staging, name),
          overwrite = FALSE))
      stop("Study 06 promotion source missing or ambiguous: ",
        name, call. = FALSE)
  }
  file.copy("studies/06_ld_operator/config.R",
    file.path(staging, "config.R"))
  file.copy("scripts/run_study06_ld_operator.R",
    file.path(staging, "reproduce.R"))
  file.copy("studies/06_ld_operator/contract_smoke_test.R",
    file.path(staging, "contract_smoke_test.R"))
  file.copy(file.path("results", "reference", "02_prediction",
    "st-bayesc-bayesr-one-replicate-development-v1",
    "example_data_manifest.csv"),
    file.path(staging, "example_data_manifest.csv"))
  source_files <- c(
    "studies/06_ld_operator/config.R",
    "studies/06_ld_operator/blocks.R",
    "studies/06_ld_operator/operators.R",
    "studies/06_ld_operator/operator_validation.R",
    "studies/06_ld_operator/simulation.R",
    "studies/06_ld_operator/methods.R",
    "studies/06_ld_operator/chain_extraction.R",
    "studies/06_ld_operator/diagnostics.R",
    "studies/06_ld_operator/metrics.R",
    "studies/06_ld_operator/pilot.R",
    "studies/06_ld_operator/targets.R",
    "studies/06_ld_operator/promotion.R",
    "scripts/run_study06_ld_operator.R")
  write.csv(data.frame(file = source_files,
    md5 = unname(.study02_canonical_md5(source_files))),
    file.path(staging, "source_files.csv"), row.names = FALSE)
  write.csv(data.frame(installed_commit =
      packageDescription("sblr")$RemoteSha,
    consulted_path = config$source_audit,
    access = "read_only_git_show_or_read_only_worktree",
    sibling_modified = FALSE, stringsAsFactors = FALSE),
    file.path(staging, "interface_audit_sources.csv"),
    row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()),
    file.path(staging, "session_info.txt"))
  writeLines(c("# Study 06: single-trait LD operator benchmark", "",
    "Package-specific development evidence for the installed sblr implementation.",
    "Block eigen denotes an operator representation, not a distinct statistical model.",
    "Five paired simulations provide descriptive evidence, not universal validation.",
    "Native fit objects and target caches are excluded."),
    file.path(staging, "README.md"))
  manifest <- c(.study06_manifest(type, dirname(status_source), config),
    list(installed_sblr_version =
        as.character(packageVersion("sblr")),
      installed_sblr_commit = packageDescription("sblr")$RemoteSha,
      installed_sblrbench_version =
        as.character(packageVersion("sblrbench")),
      repository_commit = system2("git", c("rev-parse", "HEAD"),
        stdout = TRUE),
      qgdata_commit = config$example_data$commit,
      selected_block_design =
        utils::read.csv(file.path(staging,
          "selected_block_design.csv"),
          stringsAsFactors = FALSE),
      selected_hard_filter = if (file.exists(file.path(staging,
        "selected_hard_filter.csv")))
          utils::read.csv(file.path(staging,
            "selected_hard_filter.csv"),
            stringsAsFactors = FALSE) else
          utils::read.csv(file.path(staging,
            "filter_specification.csv"),
            stringsAsFactors = FALSE),
      thread_controls = list(OMP_NUM_THREADS = 4L,
        OMP_THREAD_LIMIT = 4L, OPENBLAS_NUM_THREADS = 1L,
        MKL_NUM_THREADS = 1L, VECLIB_MAXIMUM_THREADS = 1L),
      source_status = "cached_checksum_validated_offline",
      provenance = list(installed_package_only = TRUE,
        sibling_source_access = "read_only_audit"),
      validation_status = "complete_grid_validated",
      failures = unname(status$error_message[
        status$status != "ok"])))
  jsonlite::write_json(manifest,
    file.path(staging, "benchmark_manifest.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
  write.csv(.study06_checksums(staging),
    file.path(staging, "checksums.csv"), row.names = FALSE)
  validator(staging)
  dir.create(dirname(destination), recursive = TRUE,
    showWarnings = FALSE)
  if (!file.rename(staging, destination))
    stop("Atomic Study 06 capsule promotion failed.", call. = FALSE)
  validator(destination)
  invisible(destination)
}

.study06_promote_convergence <- function(operator_dir,
                                          convergence_dir, config)
  .study06_promote("convergence",
    c(convergence_dir, operator_dir),
    config$convergence_capsule, config)

.study06_promote_benchmark <- function(benchmark_dir, config)
  .study06_promote("benchmark", benchmark_dir,
    config$benchmark_capsule, config)
