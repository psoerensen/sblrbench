source(file.path("studies", "02_prediction", "promotion.R"), local = TRUE)

.study01_current_required <- c(
  "README.md", "benchmark_manifest.json", "benchmark_summary.csv",
  "marker_metrics.csv", "credible_set_metrics.csv", "credible_set_summary.csv",
  "computational_summary.csv", "replicate_status.csv", "simulation_summary.csv",
  "seed_registry.csv", "target_warnings.csv", "ld_warning_validation.csv",
  "example_data_manifest.csv", "source_files.csv",
  "session_info.txt", "checksums.csv"
)

.study01_current_checksums <- function(path) {
  files <- sort(setdiff(list.files(path, recursive = FALSE), "checksums.csv"))
  info <- file.info(file.path(path, files))
  data.frame(file = files, size_bytes = info$size,
    md5 = unname(.study02_canonical_md5(file.path(path, files))),
    stringsAsFactors = FALSE)
}

.study01_validate_current_capsule <- function(path) {
  missing <- setdiff(.study01_current_required, list.files(path))
  if (length(missing)) stop("Study 01 current capsule is missing: ",
    paste(missing, collapse = ", "), call. = FALSE)
  status <- utils::read.csv(file.path(path, "replicate_status.csv"),
    stringsAsFactors = FALSE)
  expected <- expand.grid(replicate = seq_len(10L),
    method = c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"),
    stringsAsFactors = FALSE)
  key <- function(x) paste(x$replicate, x$method, sep = "::")
  if (nrow(status) != 40L || anyDuplicated(key(status)) ||
      !setequal(key(status), key(expected)) || !all(status$status == "ok"))
    stop("Study 01 current capsule does not contain the exact successful 40-fit grid.", call. = FALSE)
  manifest <- jsonlite::read_json(file.path(path, "benchmark_manifest.json"),
    simplifyVector = TRUE)
  if (!identical(manifest$sblr_source_commit,
        "02e8c74baa906e83c4a08d42a9cc6339b4e81072") ||
      !identical(as.integer(manifest$expected_fit_count), 40L) ||
      !identical(manifest$validation_status, "complete_grid_validated"))
    stop("Study 01 current capsule provenance or completion contract is invalid.", call. = FALSE)
  ld_validation <- utils::read.csv(file.path(path, "ld_warning_validation.csv"),
    stringsAsFactors = FALSE)
  if (!nrow(ld_validation) || any(!ld_validation$threshold_consistent) ||
      any(ld_validation$direct_edges_at_r2_threshold != 0L))
    stop("Study 01 sparse-LD warning validation is incomplete or inconsistent.", call. = FALSE)
  checks <- utils::read.csv(file.path(path, "checksums.csv"), stringsAsFactors = FALSE)
  paths <- file.path(path, checks$file)
  if (any(!file.exists(paths)) || any(unname(.study02_canonical_md5(paths)) != checks$md5))
    stop("Study 01 current capsule checksum validation failed.", call. = FALSE)
  invisible(TRUE)
}

.study01_write_target_warnings <- function(store, path) {
  meta <- targets::tar_meta(store = store, fields = c(name, warnings, error))
  meta <- meta[grepl("^method_run_", meta$name) & !is.na(meta$warnings), , drop = FALSE]
  rows <- lapply(seq_len(nrow(meta)), function(i) {
    object <- targets::tar_read_raw(meta$name[[i]], store = store)
    data.frame(target = meta$name[[i]],
      architecture = object$computational$architecture,
      replicate = object$computational$replicate,
      method = object$computational$method,
      warning = meta$warnings[[i]], stringsAsFactors = FALSE)
  })
  out <- if (length(rows)) do.call(rbind, rows) else data.frame(
    target = character(), architecture = character(), replicate = integer(),
    method = character(), warning = character(), stringsAsFactors = FALSE)
  out <- out[order(out$replicate, out$method, out$target), , drop = FALSE]
  utils::write.csv(out, path, row.names = FALSE)
  path
}

.study01_validate_sparse_warning_loci <- function(store, path, r2_threshold = 0.001) {
  meta <- targets::tar_meta(store = store, fields = c(name, warnings))
  meta <- meta[grepl("^method_run_", meta$name) & !is.na(meta$warnings), , drop = FALSE]
  Glist <- targets::tar_read(sparse_ld_glist, store = store)
  Z <- targets::tar_read(standardized_genotypes, store = store)
  csr <- sblr::sparseLD_read_CSR(Glist$sparseLD$prefix, one_based = TRUE)
  rows <- list()
  for (i in seq_len(nrow(meta))) {
    object <- targets::tar_read_raw(meta$name[[i]], store = store)
    fit <- object$fit$result$native_fit
    map <- sblr:::.stblr_marker_map_from_Glist(Glist, fit)
    sets <- object$credible_sets$native$locus_sets
    for (locus in names(sets)) {
      ids <- sets[[locus]]
      idx <- map$index[match(ids, map$marker)]
      sparse <- suppressWarnings(sblr:::.extract_sparseLD_region_dense(csr, idx, ids))
      sparse_off <- abs(sparse[upper.tri(sparse)])
      if (!length(sparse_off) || max(sparse_off) != 0) next
      direct <- stats::cor(Z[, ids, drop = FALSE])
      direct_off <- abs(direct[upper.tri(direct)])
      rows[[length(rows) + 1L]] <- data.frame(target = meta$name[[i]],
        replicate = object$computational$replicate,
        method = object$computational$method, locus = locus,
        marker_count = length(ids), sparse_max_abs_r = max(sparse_off),
        direct_max_abs_r = max(direct_off),
        direct_edges_at_r2_threshold = sum(direct_off^2 >= r2_threshold),
        r2_threshold = r2_threshold,
        threshold_consistent = sum(direct_off^2 >= r2_threshold) == 0L,
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$replicate, out$method, out$locus), , drop = FALSE]
  if (any(!out$threshold_consistent))
    stop("A Study 01 empty sparse-LD region has a retained direct LD edge.", call. = FALSE)
  utils::write.csv(out, path, row.names = FALSE)
  path
}

.study01_promote_current <- function(source_dir,
    destination = file.path("results", "reference", "01_finemapping", "current")) {
  if (dir.exists(destination)) stop("Study 01 current capsule already exists.", call. = FALSE)
  staging <- file.path(Sys.getenv("SBLR_BENCH_REFRESH_ROOT",
    file.path("results", "local", "current_benchmark_refresh")),
    "promotion_staging", "study01-current")
  if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  output_files <- c("benchmark_summary.csv", "marker_metrics.csv",
    "credible_set_metrics.csv", "credible_set_summary.csv", "computational_summary.csv",
    "replicate_status.csv", "simulation_summary.csv", "seed_registry.csv",
    "target_warnings.csv", "ld_warning_validation.csv")
  copied <- file.copy(file.path(source_dir, output_files), file.path(staging, output_files))
  if (!all(copied)) stop("Failed to copy Study 01 current result tables.", call. = FALSE)
  manifest <- jsonlite::read_json(file.path(source_dir, "pilot_manifest.json"),
    simplifyVector = FALSE)
  warning_inventory <- utils::read.csv(file.path(source_dir, "target_warnings.csv"),
    stringsAsFactors = FALSE)
  manifest$warning_target_count <- nrow(warning_inventory)
  manifest$warnings <- unique(warning_inventory$warning)
  ld_validation <- utils::read.csv(file.path(source_dir, "ld_warning_validation.csv"),
    stringsAsFactors = FALSE)
  manifest$sparse_ld_warning_validation <- list(
    empty_region_count = nrow(ld_validation),
    threshold_consistent_count = sum(ld_validation$threshold_consistent),
    mismatch_count = sum(!ld_validation$threshold_consistent),
    maximum_direct_abs_r = max(ld_validation$direct_max_abs_r),
    r2_threshold = unique(ld_validation$r2_threshold))
  jsonlite::write_json(manifest, file.path(staging, "benchmark_manifest.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null")
  file.copy(file.path("results", "reference", "02_prediction", "current",
    "example_data_manifest.csv"), staging)
  source_files <- c("_targets.R", "studies/01_finemapping/config.R",
    "studies/01_finemapping/targets.R", "studies/01_finemapping/pilot.R",
    "studies/01_finemapping/promotion.R", "studies/01_finemapping/setup_example_data.R",
    "scripts/run_current_benchmark_refresh.R", "R/metrics.R", "R/alignment.R")
  src <- data.frame(file = source_files,
    md5 = unname(.study02_canonical_md5(source_files)), stringsAsFactors = FALSE)
  utils::write.csv(src, file.path(staging, "source_files.csv"), row.names = FALSE)
  session <- sub("[[:space:]]+$", "", capture.output(utils::sessionInfo()))
  writeLines(session, file.path(staging, "session_info.txt"), useBytes = TRUE)
  writeLines(c("# Current separated-locus fine-mapping benchmark", "",
    "Fresh 10-replicate, 40-fit evidence generated with the pinned current sblr source.",
    "The established Study 01 one-chain policy is explicit in benchmark_manifest.json.",
    "Native fits and target objects are excluded; reproduce with:", "",
    "    powershell -ExecutionPolicy Bypass -File scripts/run_current_benchmark_refresh.ps1 -Phase study01 -Resume"),
    file.path(staging, "README.md"))
  utils::write.csv(.study01_current_checksums(staging),
    file.path(staging, "checksums.csv"), row.names = FALSE)
  .study01_validate_current_capsule(staging)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging, destination)) stop("Atomic Study 01 promotion failed.", call. = FALSE)
  .study01_validate_current_capsule(destination)
  invisible(destination)
}
