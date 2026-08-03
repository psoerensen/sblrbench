.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-checkpoints.R"), local = TRUE)

.study06_paths <- function(config) list(
  glist_path = Sys.getenv("SBLR_BENCH_GLIST", ""),
  data_dir = Sys.getenv("SBLR_BENCH_DATA_DIR",
    file.path("results", "local", "02_prediction", "data")),
  genotype_output_dir = Sys.getenv("SBLR_BENCH_LD_DIR",
    file.path("results", "local", "02_prediction", "genotype_setup")),
  local_dir = Sys.getenv("SBLR_BENCH_STUDY06_LOCAL",
    config$local_dir),
  operator_dir = file.path(Sys.getenv("SBLR_BENCH_STUDY06_LOCAL",
    config$local_dir), "operators"))

.study06_load_scaled_genotypes <- function(base_glist, chr, sample_ids,
                                            marker_ids, split) {
  prior_store <- file.path("results", "local", "five_replicate_overnight",
    "_targets_02")
  if (dir.exists(prior_store)) {
    cached <- try(targets::tar_read(prediction_scaled_genotypes,
      store = prior_store), silent = TRUE)
    if (!inherits(cached, "try-error") &&
        identical(rownames(cached$all), sample_ids) &&
        identical(colnames(cached$all), marker_ids) &&
        identical(rownames(cached$train), split$train_ids) &&
        identical(rownames(cached$test), split$test_ids))
      return(cached)
  }
  raw <- sblrbench:::benchmark_extract_raw_genotypes(base_glist, chr,
    sample_ids, marker_ids)
  sblrbench::training_scaled_genotypes(raw, split$train_rows)
}

.study06_summary_stats <- function(simulation, Glist, split, config) {
  started <- proc.time()[["elapsed"]]
  y <- simulation$phenotype[split$train_ids, , drop = FALSE]
  stats <- sblr::make_summary_stats(Glist = Glist, y = y,
    chr = config$chr, rows = split$train_rows, scale = TRUE,
    nthreads = 1L)
  if (!identical(stats$marker_names, simulation$marker_ids) ||
      !identical(stats$trait_names, config$trait) ||
      !identical(as.integer(stats$n), length(split$train_ids)) ||
      !identical(stats$rows, sort(split$train_rows)) ||
      !isTRUE(all.equal(unname(stats$af[[1L]]),
        unname(Glist$sparseLD$af[[1L]]), tolerance = 0)))
    stop("Study 06 training statistic/alignment contract failed.",
      call. = FALSE)
  list(stats = stats,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
}

.study06_select_block_design <- function(candidate_table, config) {
  feasible <- candidate_table[
    candidate_table$gate_construction_status == "passed" &
      candidate_table$runtime_diagonal_minimum > 0 &
      is.finite(candidate_table$estimated_dense_storage_bytes), ,
    drop = FALSE]
  if (!nrow(feasible))
    stop("No feasible Study 06 block design met the frozen criteria.",
      call. = FALSE)
  selected <- feasible[order(
    -feasible$fraction_full_csr_squared_ld_mass_retained,
    feasible$estimated_dense_storage_bytes), ,
    drop = FALSE][1L, ]
  if (selected$target_size != config$selected_block_size)
    stop("Selected block design differs from the frozen config: ",
      selected$target_size, " versus ", config$selected_block_size,
      call. = FALSE)
  selected
}

.study06_filter_candidate <- function(unfiltered, hard, blocks, config,
                                      tau) {
  metrics <- .study06_operator_metrics(unfiltered, hard, blocks,
    config, operator_id = paste0("hard_tau_", tau),
    reference_id = "block_eigen_unfiltered")
  data.frame(
    eigen_tau = tau,
    requested_threshold = tau,
    effective_threshold = max(tau, 0.01),
    retained_rank = sum(metrics$retained_rank),
    original_rank = sum(metrics$original_rank),
    retained_rank_proportion =
      sum(metrics$retained_rank) / sum(metrics$original_rank),
    retained_positive_eigenvalue_mass =
      sum(metrics$retained_positive_eigenvalue_mass),
    positive_eigenvalue_mass = sum(metrics$positive_eigenvalue_mass),
    retained_mass_proportion =
      sum(metrics$retained_positive_eigenvalue_mass) /
        sum(metrics$positive_eigenvalue_mass),
    eigenvalue_mass_retained =
      sum(metrics$retained_positive_eigenvalue_mass) /
        sum(metrics$positive_eigenvalue_mass),
    minimum_runtime_diagonal =
      min(metrics$runtime_diagonal_minimum),
    operator_frobenius_maximum_error =
      max(metrics$reconstruction_maximum_absolute_error),
    matrix_vector_maximum_error =
      max(metrics$matrix_vector_maximum_error),
    quadratic_form_maximum_error =
      max(metrics$quadratic_form_maximum_error),
    eigenvalues_removed = sum(metrics$original_rank) -
      sum(metrics$retained_rank),
    filter_activity_status = if (sum(metrics$retained_rank) <
      sum(metrics$original_rank)) "active" else "effective_no_op",
    prediction_check_status =
      "deferred_to_selected_filter_convergence_pilot",
    pass = min(metrics$runtime_diagonal_minimum) > 1e-8 &&
      all(vapply(metrics, function(x) !is.numeric(x) ||
        all(is.finite(x)), logical(1))),
    stringsAsFactors = FALSE)
}

.study06_select_hard_filter <- function(candidates, selected_tau = 0.10) {
  selected <- candidates[candidates$pass &
    abs(candidates$eigen_tau - selected_tau) < 1e-12, , drop = FALSE]
  if (nrow(selected) != 1L)
    stop("The frozen public hard-filter threshold is unavailable or invalid.",
      call. = FALSE)
  selected
}

.study06_select_fixed_ridge <- function(candidates, selected_shrinkage = 0.01) {
  selected <- candidates[candidates$pass & abs(candidates$shrinkage_weight -
    selected_shrinkage) < 1e-12, , drop = FALSE]
  if (nrow(selected) != 1L)
    stop("The frozen scale-aware fixed-ridge policy is unavailable.",
      call. = FALSE)
  selected
}

.study06_fixed_ridge_candidate <- function(unfiltered, candidate, blocks,
                                           config, shrinkage_weight) {
  metrics <- .study06_operator_metrics(unfiltered, candidate, blocks, config,
    operator_id = paste0("ridge_fixed_", format(shrinkage_weight,
      scientific = FALSE)), reference_id = "block_eigen_unfiltered")
  data.frame(shrinkage_weight = shrinkage_weight,
    eigen_eta = shrinkage_weight / (1 - shrinkage_weight),
    minimum_runtime_diagonal = min(metrics$runtime_diagonal_minimum),
    frobenius_relative_error = sqrt(sum(
      metrics$frobenius_relative_error^2 * metrics$frobenius_norm^2)) /
      max(sqrt(sum(metrics$frobenius_norm^2)), .Machine$double.eps),
    maximum_element_error = max(metrics$reconstruction_maximum_absolute_error),
    maximum_matrix_vector_error = max(metrics$matrix_vector_maximum_error),
    maximum_quadratic_form_error = max(metrics$quadratic_form_maximum_error),
    off_diagonal_norm_retained = mean(metrics$off_diagonal_norm_retained),
    pass = min(metrics$runtime_diagonal_minimum) > 0 &&
      all(is.finite(metrics$frobenius_relative_error)) &&
      max(metrics$reconstruction_maximum_absolute_error) > 0,
    selection_basis = "frozen 1% off-diagonal shrinkage; not prediction optimized",
    stringsAsFactors = FALSE)
}

.study06_lw_audit <- function(unfiltered, lw, blocks, config, sample_size) {
  metrics <- .study06_operator_metrics(unfiltered, lw, blocks, config,
    operator_id = "block_eigen_ridge_lw_sensitivity",
    reference_id = "block_eigen_unfiltered")
  shrink <- lw$diagnostics$shrink
  reference_blocks <- .study06_dense_blocks(unfiltered)
  lw_blocks <- .study06_dense_blocks(lw)
  delta_sq <- sum(vapply(seq_along(reference_blocks), function(k)
    sum((lw_blocks[[k]] - reference_blocks[[k]])^2), numeric(1)))
  reference_sq <- sum(vapply(reference_blocks, function(A)
    sum(A^2), numeric(1)))
  data.frame(
    formula = "a=min(bbar,d2)/d2; A_tilde=(1-a)A+a*diag(diag(A))",
    shrinkage_target = "diagonal matrix diag(diag(A))",
    shrinkage_scale = "standardized-genotype cross-product converted internally to correlation scale for a",
    sample_size = sample_size,
    shrinkage_minimum = min(shrink), shrinkage_mean = mean(shrink),
    shrinkage_median = stats::median(shrink), shrinkage_maximum = max(shrink),
    coefficient_one_semantics = "complete off-diagonal shrinkage",
    diagonal_maximum_absolute_change = max(abs(lw$diagonal -
      unfiltered$diagonal)),
    trace_difference = sum(metrics$trace_difference),
    frobenius_relative_error = sqrt(delta_sq) /
      max(sqrt(reference_sq), .Machine$double.eps),
    off_diagonal_norm_retained = mean(metrics$off_diagonal_norm_retained),
    maximum_matrix_vector_error = max(metrics$matrix_vector_maximum_error),
    retained_as_main_configuration = FALSE,
    disposition = "deterministic sensitivity only; replaced by frozen 1% fixed ridge",
    stringsAsFactors = FALSE)
}

.study06_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (is.data.frame(x) && nrow(x) && ncol(x))
    x <- x[do.call(order, unname(x)), , drop = FALSE]
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

.study06_atomic_rds <- function(x, path) {
  benchmark_atomic_save_rds(x, path, compress = FALSE,
    temporary_prefix = ".study06-fit-")
  path
}

.study06_validate_cached_run <- function(run, method, replicate) {
  identical(run$status, "ok") &&
    identical(run$architecture, method$architecture) &&
    identical(run$replicate, as.integer(replicate)) &&
    identical(run$method$configuration, method$configuration) &&
    length(run$fit$chains) == 4L &&
    identical(sort(vapply(run$fit$chains, `[[`, integer(1),
      "chain_index")), 1:4)
}

.study06_run_cached_fit <- function(method, simulation, stats, Glist,
                                    split, blocks, runtime_glist, config,
                                    phase, recommendations, hard_tau,
                                    ridge_eta = config$selected_ridge_eta) {
  root <- .study06_paths(config)$local_dir
  dir <- file.path(root, "fit_checkpoints", phase)
  path <- file.path(dir, paste0(method$architecture, "__",
    simulation$replicate, "__", method$configuration, ".rds"))
  if (file.exists(path)) {
    old <- try(readRDS(path), silent = TRUE)
    if (!inherits(old, "try-error") &&
        .study06_validate_cached_run(old, method,
          simulation$replicate)) return(old)
  }
  run <- .study06_fit(method, simulation, stats, Glist, split,
    blocks, runtime_glist, config, phase = phase,
    block_recommendations = recommendations,
    selected_hard_tau = hard_tau, selected_ridge_eta = ridge_eta)
  status_dir <- file.path(root, "fit_status", phase)
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
  status_path <- file.path(status_dir, paste0(method$architecture, "__",
    simulation$replicate, "__", method$configuration, ".csv"))
  .study06_write_csv(data.frame(
    phase = phase, architecture = method$architecture,
    replicate = simulation$replicate,
    configuration = method$configuration,
    method = method$native_method,
    operator = method$operator_family,
    filter = method$filter_policy, state = run$status,
    start_time = format(run$started_at, tz = "UTC", usetz = TRUE),
    finish_time = format(run$finished_at, tz = "UTC", usetz = TRUE),
    elapsed_seconds = run$runtime, output_path = path,
    validation_status = if (run$status == "ok") "passed" else "failed",
    error_message = run$reason, stringsAsFactors = FALSE), status_path)
  if (!identical(run$status, "ok"))
    stop("Study 06 fit failed; status preserved at ", status_path,
      ": ", run$reason, call. = FALSE)
  .study06_atomic_rds(run, path)
  run
}

.study06_seed_registry <- function(config) {
  grid <- expand.grid(architecture = config$architectures,
    replicate = seq_len(config$replicate_count),
    configuration = config$configurations, chain = 1:4,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$data_selection_seed <- config$seeds$data_selection
  grid$simulation_seed <- mapply(.study06_seed,
    MoreArgs = list(config = config),
    architecture = grid$architecture,
    replicate = grid$replicate)
  grid$fit_seed <- mapply(.study06_seed,
    MoreArgs = list(config = config),
    architecture = grid$architecture,
    replicate = grid$replicate,
    configuration = grid$configuration)
  grid$chain_seed <- mapply(.study06_seed,
    MoreArgs = list(config = config),
    architecture = grid$architecture,
    replicate = grid$replicate,
    configuration = grid$configuration,
    chain = grid$chain)
  grid
}
