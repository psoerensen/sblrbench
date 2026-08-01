.study07_promotion_root <- if (file.exists(file.path("studies",
  "02_prediction", "promotion.R"))) "." else file.path("..", "..")
source(file.path(.study07_promotion_root, "studies", "02_prediction",
  "promotion.R"), local = TRUE)

.study07_optional_flag <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) isTRUE(default) else
    isTRUE(x[[1L]])
}

.study07_required <- function(type = c("contract", "convergence", "benchmark")) {
  type <- match.arg(type)
  common <- c("README.md", "benchmark_manifest.json", "config.R",
    "example_data_manifest.csv", "source_files.csv",
    "interface_audit_sources.csv", "session_info.txt", "reproduce.R",
    "contract_smoke_test.R", "checksums.csv")
  switch(type,
    contract = c(common, "interface_audit.csv", "output_semantics.csv",
      "joint_state_map.csv", "permutation_contracts.csv",
      "deterministic_simulation_summary.csv", "tiny_fit_status.csv",
      "operator_equivalence.csv", "runtime_scaling.csv",
      "runtime_projections.csv", "selected_marker_count.csv",
      "computational_limits.csv"),
    convergence = c(common, "joint_state_map.csv", "selected_marker_count.csv", "fit_status.csv",
      "convergence_diagnostics.csv", "candidate_settings.csv",
      "method_recommendations.csv"),
    benchmark = c(common, "selected_marker_count.csv", "fit_status.csv",
      "simulation_truth.csv", "prediction_metrics.csv",
      "trait_marker_metrics.csv", "parameter_estimates.csv",
      "internal_consistency.csv", "paired_replicate_differences.csv",
      "paired_comparison_summary.csv", "computational_summary.csv",
      "joint_state_map.csv", "method_output_availability.csv",
      "seed_registry.csv"))
}

.study07_checksums <- function(path) {
  files <- sort(setdiff(list.files(path, recursive = FALSE), "checksums.csv"))
  info <- file.info(file.path(path, files))
  data.frame(file = files, size_bytes = info$size,
    md5 = unname(.study02_canonical_md5(file.path(path, files))),
    stringsAsFactors = FALSE)
}

.study07_validate_checksums <- function(path, required) {
  x <- utils::read.csv(file.path(path, "checksums.csv"),
    stringsAsFactors = FALSE)
  if (anyDuplicated(x$file) || any(x$file != basename(x$file)) ||
      any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])\\.\\.([/\\\\]|$))",
        x$file)) || !setequal(x$file, setdiff(required, "checksums.csv")) ||
      any(unname(.study02_canonical_md5(file.path(path, x$file))) != x$md5))
    stop("Study 07 canonical checksum validation failed.", call. = FALSE)
  invisible(TRUE)
}

.study07_validate_capsule <- function(path,
                                      type = c("contract", "convergence", "benchmark")) {
  type <- match.arg(type); required <- .study07_required(type)
  if (!dir.exists(path) || length(setdiff(required, list.files(path))))
    stop("Study 07 capsule is incomplete: ", type, call. = FALSE)
  .study07_validate_checksums(path, required)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  if (!identical(manifest$benchmark_status, "complete") ||
      !identical(manifest$study, "07_mtblr_validation"))
    stop("Study 07 manifest is incomplete.", call. = FALSE)
  if (type == "contract") {
    states <- read.csv(file.path(path, "joint_state_map.csv"))
    tiny <- read.csv(file.path(path, "tiny_fit_status.csv"))
    gate <- read.csv(file.path(path, "operator_equivalence.csv"))
    selected <- read.csv(file.path(path, "selected_marker_count.csv"))
    if (!identical(states$internal_key, c("0_0", "1_0", "0_1", "1_1")) ||
        nrow(tiny) != 4L || any(tiny$status != "ok") ||
        any(tiny$chain_count != 4L) || !isTRUE(gate$pass[[1L]]) ||
        selected$marker_count[[1L]] < 2000L)
      stop("Study 07 contract/runtime capsule validation failed.")
  } else if (type == "convergence") {
    status <- read.csv(file.path(path, "fit_status.csv"))
    rec <- read.csv(file.path(path, "method_recommendations.csv"))
    if (nrow(status) != 4L || any(status$status != "ok") ||
        nrow(rec) != 4L || any(rec$recommendation_status != "available") ||
        any(rec$nchains != 4L) || any(rec$nthin != 1L))
      stop("Study 07 convergence capsule validation failed.")
  } else {
    status <- read.csv(file.path(path, "fit_status.csv"))
    identity <- read.csv(file.path(path, "internal_consistency.csv"))
    paired <- read.csv(file.path(path, "paired_replicate_differences.csv"))
    expected <- expand.grid(architecture = c("independent", "partially_shared"),
      replicate = 1:5, implementation = c("mt_bed_bayesc",
        "mt_csr_sbayesc", "mt_block_eigen_sbayesc"),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    key <- function(x) paste(x$architecture, x$replicate,
      x$implementation, sep = "|")
    if (nrow(status) != 30L || anyDuplicated(key(status)) ||
        !setequal(key(status), key(expected)) || any(status$status != "ok") ||
        any(status$chain_count != 4L) || any(!identity$passed) ||
        any(!paired$complete_pair))
      stop("Study 07 benchmark capsule validation failed.")
  }
  invisible(TRUE)
}

.study07_seed_registry <- function(config) {
  grid <- expand.grid(architecture = config$main_architectures,
    replicate = seq_len(config$replicate_count),
    implementation = config$implementations, chain = 1:4,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  values <- Map(function(a, r, i) .study07_fit_seeds(config, a, r, i, 4L),
    grid$architecture, grid$replicate, grid$implementation)
  grid$sample_selection_seed <- config$seeds$sample_selection
  grid$marker_selection_seed <- config$seeds$marker_selection
  grid$simulation_seed <- mapply(.study07_simulation_seed,
    grid$architecture, grid$replicate, MoreArgs = list(config = config))
  grid$fit_seed <- vapply(values, `[[`, integer(1), "fit_seed")
  grid$chain_seed <- mapply(function(x, chain) x$chain_seeds[[chain]],
    values, grid$chain)
  grid
}

.study07_promote <- function(type = c("contract", "convergence", "benchmark"),
                             source_dir, destination, config) {
  type <- match.arg(type); validator <- function(path)
    .study07_validate_capsule(path, type)
  if (dir.exists(destination)) { validator(destination); return(destination) }
  required <- .study07_required(type)
  staging <- file.path(config$local_dir, "promotion_staging",
    paste0(basename(destination), "-", Sys.getpid()))
  if (dir.exists(staging)) stop("Study 07 promotion staging already exists.")
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  generated <- setdiff(required, c("README.md", "benchmark_manifest.json",
    "config.R", "example_data_manifest.csv", "source_files.csv",
    "interface_audit_sources.csv", "session_info.txt", "reproduce.R",
    "contract_smoke_test.R", "checksums.csv",
    "method_output_availability.csv", "seed_registry.csv"))
  if (type != "contract") generated <- setdiff(generated, "joint_state_map.csv")
  for (name in generated) if (!file.copy(file.path(source_dir, name),
      file.path(staging, name), overwrite = FALSE))
    stop("Missing Study 07 promotion source: ", name)
  file.copy("studies/07_mtblr_validation/config.R",
    file.path(staging, "config.R"))
  file.copy("scripts/run_study07_mtblr_validation.R",
    file.path(staging, "reproduce.R"))
  file.copy("studies/07_mtblr_validation/contract_smoke_test.R",
    file.path(staging, "contract_smoke_test.R"))
  file.copy(file.path("results", "reference", "06_ld_operator",
    "st-ld-operator-five-replicate-development-v1",
    "example_data_manifest.csv"), file.path(staging,
      "example_data_manifest.csv"))
  if (type != "contract") .study07_write_csv(
    .study07_state_table(config$trait_names),
    file.path(staging, "joint_state_map.csv"))
  if (type == "benchmark") {
    .study07_write_csv(data.frame(
      implementation = config$implementations,
      marker_joint_state_probabilities = FALSE,
      marker_joint_state_reason =
        "not returned for BayesC by installed implementation",
      trait_specific_pip = TRUE, global_state_probabilities = TRUE,
      covariance_draws = TRUE, stringsAsFactors = FALSE),
      file.path(staging, "method_output_availability.csv"))
    .study07_write_csv(.study07_seed_registry(config),
      file.path(staging, "seed_registry.csv"))
  }
  sources <- c(list.files("studies/07_mtblr_validation", full.names = TRUE),
    "scripts/run_study07_mtblr_validation.R",
    "scripts/run_study07_mtblr_validation.ps1")
  write.csv(data.frame(file = sources,
    md5 = unname(.study02_canonical_md5(sources))),
    file.path(staging, "source_files.csv"), row.names = FALSE)
  write.csv(data.frame(installed_commit = packageDescription("sblr")$RemoteSha,
    consulted_path = config$consulted_sblr_sources,
    access = "read_only_sibling_source_audit", sibling_modified = FALSE,
    stringsAsFactors = FALSE), file.path(staging,
      "interface_audit_sources.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()),
    file.path(staging, "session_info.txt"))
  writeLines(c("# Study 07 multivariate implementation validation", "",
    "Reduced-marker, two-trait BayesC/SBayesC development evidence.",
    "Residual covariance is diagonal and generating residual covariance is zero.",
    "Both traits use identical individuals; sample-overlap modelling is not assessed.",
    "Five replicates provide descriptive evidence, not universal validation."),
    file.path(staging, "README.md"))
  status <- if (file.exists(file.path(staging, "fit_status.csv")))
    read.csv(file.path(staging, "fit_status.csv")) else NULL
  manifest <- list(study = "07_mtblr_validation",
    task = switch(type, contract = "mtblr_contract_runtime_validation",
      convergence = "mtblr_convergence_selection",
      benchmark = "mtblr_five_replicate_validation"),
    benchmark_scope = paste0(type, "_development"),
    benchmark_status = "complete", trait_count = 2L,
    model_family = "BayesC/SBayesC", residual_covariance = "diagonal",
    identical_samples = TRUE, sample_overlap_modeling = FALSE,
    marker_count = if (file.exists(file.path(staging,
      "selected_marker_count.csv"))) read.csv(file.path(staging,
        "selected_marker_count.csv"))$marker_count[[1L]] else NA_integer_,
    replicate_count = if (type == "benchmark") 5L else 1L,
    expected_fit_count = switch(type, contract = 4L,
      convergence = 4L, benchmark = 30L),
    successful_fit_count = if (is.null(status)) 4L else sum(status$status == "ok"),
    failed_fit_count = if (is.null(status)) 0L else sum(status$status != "ok"),
    chains_per_fit = 4L,
    installed_sblr_version = as.character(packageVersion("sblr")),
    installed_sblr_commit = packageDescription("sblr")$RemoteSha,
    installed_sblrbench_version = as.character(packageVersion("sblrbench")),
    repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    qgdata_commit = config$example_data$commit,
    source_status = "cached_checksum_validated_offline",
    validation_status = "complete_grid_validated")
  jsonlite::write_json(manifest, file.path(staging,
    "benchmark_manifest.json"), pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null")
  write.csv(.study07_checksums(staging), file.path(staging,
    "checksums.csv"), row.names = FALSE)
  validator(staging)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging, destination))
    stop("Atomic Study 07 capsule promotion failed.")
  validator(destination)
  invisible(destination)
}

.study07_promote_contract <- function(config)
  .study07_promote("contract", .study07_paths(config)$contract_output,
    config$contract_capsule, config)
.study07_promote_convergence <- function(config)
  .study07_promote("convergence", .study07_paths(config)$convergence_output,
    config$convergence_capsule, config)
.study07_promote_benchmark <- function(config)
  .study07_promote("benchmark", .study07_paths(config)$benchmark_output,
    config$benchmark_capsule, config)
