.study07_paths <- function(config) list(
  local_dir = config$local_dir,
  data_dir = file.path(config$local_dir, "data"),
  ld_dir = file.path(config$local_dir, "ld"),
  operator_dir = file.path(config$local_dir, "operators"),
  fit_dir = file.path(config$local_dir, "fit_checkpoints"),
  contract_output = file.path(config$local_dir, "contract_output"),
  runtime_output = file.path(config$local_dir, "runtime_output"),
  convergence_output = file.path(config$local_dir, "convergence_output"),
  benchmark_output = file.path(config$local_dir, "benchmark_output"))

.study07_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(x) && nrow(x) && ncol(x)) {
    order_columns <- vapply(x, function(z)
      is.atomic(z) && length(z) == nrow(x), logical(1))
    if (any(order_columns)) x <- x[do.call(order,
      unname(x[order_columns])), , drop = FALSE]
  }
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

.study07_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".study07-", dirname(path), ".rds")
  saveRDS(x, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) {
    unlink(tmp); stop("Atomic Study 07 checkpoint replacement failed.")
  }
  path
}

.study07_base_resources <- function(config) {
  store <- file.path(config$local_dir, "..", "study05_ld_operator", "_targets")
  if (!dir.exists(store))
    stop("Validated Study 05 local target store is unavailable.", call. = FALSE)
  read <- function(name) targets::tar_read_raw(name, store = store)
  list(base_glist = read("study05_base_glist"),
    marker_ids = read("study05_filtered_markers")$marker_ids,
    sample_ids = read("study05_sample_ids"),
    split = read("study05_split"))
}

.study07_raw_genotypes <- function(resources, config) {
  max_m <- max(config$marker_candidates)
  marker_ids <- .study07_marker_subset(resources$marker_ids, max_m)
  qgg::getG(Glist = resources$base_glist, chr = config$chr,
    ids = resources$sample_ids, rsids = marker_ids,
    impute = FALSE, scale = FALSE)
}

.study07_scaled_data <- function(raw, marker_count, training_rows) {
  marker_ids <- colnames(raw)[seq_len(marker_count)]
  x <- raw[, marker_ids, drop = FALSE]
  training_rows <- as.integer(training_rows)
  if (!is.matrix(x) || is.null(rownames(x)) || is.null(colnames(x)) ||
      !length(training_rows) || anyNA(training_rows) ||
      anyDuplicated(training_rows) ||
      any(training_rows < 1L | training_rows > nrow(x)))
    stop("Study 07 genotype scaling inputs are invalid.", call. = FALSE)
  train <- x[training_rows, , drop = FALSE]
  means <- colMeans(train, na.rm = TRUE); af <- means / 2
  scale <- sqrt(2 * af * (1 - af))
  if (any(!is.finite(c(means, af, scale))) ||
      any(af <= 0 | af >= 1) || any(scale <= 0))
    stop("Study 07 training genotype scale is invalid.", call. = FALSE)
  filled <- x; missing <- which(is.na(filled), arr.ind = TRUE)
  if (nrow(missing)) filled[missing] <- means[missing[, 2L]]
  scaled <- sweep(sweep(filled, 2L, means, "-"), 2L, scale, "/")
  test_rows <- setdiff(seq_len(nrow(x)), training_rows)
  list(allele_frequency = setNames(af, marker_ids),
    center = setNames(means, marker_ids), scale = setNames(scale, marker_ids),
    all = scaled, train = scaled[training_rows, , drop = FALSE],
    test = scaled[test_rows, , drop = FALSE], train_rows = training_rows,
    test_rows = test_rows)
}

.study07_working_glist <- function(base_glist, marker_ids, af, config) {
  sblrbench:::benchmark_set_training_af(sblrbench:::benchmark_set_glist_marker_order(base_glist,
    config$chr, marker_ids), config$chr, marker_ids, af)
}

.study07_make_ld <- function(Glist, rows, marker_ids, config, output_dir,
                             tag) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, paste0("ld_", tag, "_m", length(marker_ids)))
  cache <- paste0(prefix, "_glist.rds")
  if (file.exists(cache)) {
    old <- try(readRDS(cache), silent = TRUE)
    if (!inherits(old, "try-error") &&
        identical(old$sparseLD$rows, as.integer(rows)) &&
        identical(old$rsids[[config$chr]][old$sparseLD$cls[[1L]]], marker_ids))
      return(old)
  }
  out <- do.call(sblr::make_sparse_ld, c(list(Glist = Glist,
    rows = as.integer(rows), out_prefix = prefix, chr = config$chr),
    config$sparse_ld))
  ids <- out$rsids[[config$chr]][out$sparseLD$cls[[1L]]]
  if (!identical(ids, marker_ids) ||
      !identical(out$sparseLD$rows, as.integer(rows)))
    stop("Study 07 sparse-LD alignment failed.", call. = FALSE)
  saveRDS(out, cache)
  out
}

.study07_make_stats <- function(simulation, Glist, rows, config) {
  stats <- sblr::make_summary_stats(Glist = Glist,
    y = simulation$phenotype[rows, , drop = FALSE], chr = config$chr,
    rows = as.integer(rows), scale = TRUE, nthreads = 1L)
  if (is.list(stats$af) && is.list(stats$cls) &&
      length(stats$af) == length(stats$cls)) names(stats$af) <- names(stats$cls)
  .study07_validate_stats(stats, simulation$marker_ids,
    config$trait_names, length(rows))
  stats
}

.study07_cached_fit <- function(implementation, simulation, stats, Glist,
                                runtime_glist, rows, blocks, config,
                                controls, phase) {
  paths <- .study07_paths(config)
  path <- file.path(paths$fit_dir, phase, paste0(simulation$architecture,
    "__", simulation$replicate, "__", implementation, "__m",
    length(simulation$marker_ids), ".rds"))
  if (file.exists(path)) {
    old <- try(readRDS(path), silent = TRUE)
    if (!inherits(old, "try-error") && identical(old$status, "ok") &&
        identical(old$implementation$id, implementation) &&
        nrow(old$fit$bm) == length(simulation$marker_ids) &&
        length(old$fit$chains) == controls$nchains) return(old)
  }
  run <- .study07_fit(implementation, simulation, stats, Glist,
    runtime_glist, rows, blocks, config, controls)
  status_path <- file.path(paths$local_dir, "fit_status", phase,
    paste0(simulation$architecture, "__", simulation$replicate, "__",
      implementation, ".csv"))
  .study07_write_csv(data.frame(phase = phase,
    architecture = simulation$architecture, replicate = simulation$replicate,
    implementation = implementation, marker_count = length(simulation$marker_ids),
    chain_count = controls$nchains, state = run$status,
    start_time = format(run$started_at, tz = "UTC", usetz = TRUE),
    finish_time = format(run$finished_at, tz = "UTC", usetz = TRUE),
    elapsed_seconds = run$runtime, memory_estimate = if (run$status == "ok")
      run$fit$memory_estimate$execution_estimated_total_bytes else NA_real_,
    output_path = path, validation_status = if (run$status == "ok")
      "passed" else "failed", error_message = run$error,
    stringsAsFactors = FALSE), status_path)
  if (run$status != "ok") stop("Study 07 fit failed: ", implementation,
    " / ", simulation$architecture, " / replicate ", simulation$replicate,
    ": ", run$error, call. = FALSE)
  .study07_atomic_rds(run, path)
  run
}

.study07_contract_evidence <- function(Z, config) {
  rows <- list()
  for (architecture in config$contract_architectures) {
    sim <- .study07_simulate(Z, architecture, 1L, config)
    .study07_validate_simulation(sim, Z, config)
    perm <- .study07_permutation_contract(Z, sim$phenotype,
      sim$effects, config$seeds$permutation)
    perm$architecture <- architecture
    rows[[architecture]] <- perm
  }
  do.call(rbind, rows)
}
