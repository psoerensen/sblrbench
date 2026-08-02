.study06v2_sha256 <- function(path) digest::digest(file = path,
  algo = "sha256", serialize = FALSE)

.study06v2_summary_statistics_hash <- function(store = file.path("results",
    "local", "study06_ld_operator", "_targets")) {
  x <- targets::tar_read_raw("study06_operator_stats", store = store)
  digest::digest(list(marker_names = x$marker_names, rows = x$rows,
    af = x$af, wy = x$wy, ww = x$ww, yy = x$yy, n = x$n),
    algo = "sha256")
}

.study06v2_installed_provenance <- function(config) {
  info <- .study06v2_load_pinned_sblr(config, recompile = FALSE)
  description <- packageDescription("sblr")
  list(version = as.character(packageVersion("sblr")),
    sha = description[["RemoteSha"]], path = info$path,
    library = info$library,
    installation_command = config$installation_command,
    r_executable = normalizePath(file.path(R.home("bin"), "Rscript.exe"),
      winslash = "/", mustWork = TRUE),
    r_version = R.version.string,
    compiler = "g++.exe (GCC) 13.2.0",
    openmp = "SHLIB_OPENMP_CXXFLAGS=-fopenmp",
    blas = "Rblas", lapack = "Rlapack")
}

.study06v2_validate_preoptimization_archive <- function(config) {
  archive <- normalizePath(config$preoptimization_archive, winslash = "/",
    mustWork = TRUE)
  metadata <- file.path(archive, "ARCHIVE_METADATA.md")
  checkpoints <- list.files(file.path(archive, "fit_checkpoints",
    "convergence"), pattern = "[.]rds$", full.names = TRUE)
  if (!file.exists(metadata) || length(checkpoints) < 2L)
    stop("Pre-optimization Study 06 v2 archive is incomplete.", call. = FALSE)
  provenance <- vapply(checkpoints, function(path) {
    x <- readRDS(path)
    identical(x$source_sha, config$required_sblr_sha)
  }, logical(1))
  if (any(provenance))
    stop("Pre-optimization archive contains a checkpoint labelled with the optimized SHA.",
      call. = FALSE)
  fresh <- list.files(file.path(config$local_dir, "fit_checkpoints"),
    pattern = "[.]rds$", recursive = TRUE, full.names = TRUE)
  if (length(fresh)) {
    valid <- vapply(fresh, function(path) {
      x <- try(readRDS(path), silent = TRUE)
      !inherits(x, "try-error") &&
        identical(x$source_sha, config$required_sblr_sha) &&
        identical(x$package_version, config$required_sblr_version)
    }, logical(1))
    if (!all(valid))
      stop("Fresh checkpoint tree contains pre-optimization or invalid provenance.",
        call. = FALSE)
  }
  invisible(list(path = archive, checkpoint_count = length(checkpoints)))
}

.study06v2_write_installation_provenance <- function(config) {
  output <- file.path(config$local_dir, "preflight")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  info <- .study06v2_installed_provenance(config)
  writeLines(c(
    paste("installation_command:", info$installation_command),
    paste("R_executable:", info$r_executable),
    paste("library:", info$library), paste("package_path:", info$path),
    paste("package_version:", info$version), paste("package_commit:", info$sha),
    paste("R_version:", info$r_version), paste("compiler:", info$compiler),
    paste("OpenMP:", info$openmp), paste("BLAS:", info$blas),
    paste("LAPACK:", info$lapack)), file.path(output,
      "installed_sblr_provenance.txt"))
  invisible(info)
}

.study06v2_write_deterministic_manifest <- function(config) {
  output <- file.path(config$local_dir, "deterministic")
  low_rank_sources <- c("R/stblr-block-eigen.R", "R/RcppExports.R",
    "src/st_block_eigen.cpp", "src/st_block_low_rank.cpp",
    "src/st_cpg_omp_csr.cpp", "src/st_cpg_omp_csr_bayesr.cpp",
    "src/st_ld_operator.h", "src/blr_block_low_rank.h")
  sibling <- normalizePath(file.path("..", "sblr"), winslash = "/",
    mustWork = TRUE)
  inventory <- data.frame(category = "optimized_low_rank_source",
    path = low_rank_sources,
    sha256 = vapply(file.path(sibling, low_rank_sources),
      .study06v2_sha256, ""), stringsAsFactors = FALSE)
  input_paths <- c(file.path("results", "local", "02_prediction", "data",
    config$example_data$files), file.path(output, "block_definitions.csv"))
  input_inventory <- data.frame(category = "deterministic_input",
    path = input_paths, sha256 = vapply(input_paths, .study06v2_sha256, ""),
    stringsAsFactors = FALSE)
  stats <- data.frame(category = "summary_statistics",
    path = "targets:study06_operator_stats",
    sha256 = .study06v2_summary_statistics_hash(), stringsAsFactors = FALSE)
  .study06v2_write_csv(rbind(inventory, input_inventory, stats),
    file.path(output, "optimized_hash_inventory.csv"))
  summary <- read.csv(file.path(output, "deterministic_identity_summary.csv"),
    stringsAsFactors = FALSE)
  manifest <- list(validation_status = if (all(summary$pass)) "passed" else
      "failed", source_sha = config$required_sblr_sha,
    package_version = config$required_sblr_version,
    commit_subject = "Optimize retained low-rank scalar sampling",
    deterministic_evidence_reused = FALSE,
    low_rank_residual_rebuild_every = config$low_rank_residual_rebuild_every,
    numerical_summary_sha256 = .study06v2_sha256(file.path(output,
      "deterministic_identity_summary.csv")),
    hash_inventory_sha256 = .study06v2_sha256(file.path(output,
      "optimized_hash_inventory.csv")),
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(manifest, file.path(output,
    "deterministic_validation_manifest.json"), auto_unbox = TRUE,
    pretty = TRUE)
  if (!identical(manifest$validation_status, "passed"))
    stop("Optimized deterministic manifest validation failed.", call. = FALSE)
  invisible(manifest)
}
