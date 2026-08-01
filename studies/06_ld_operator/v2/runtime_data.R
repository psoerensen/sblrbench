.study06v2_load_runtime_data <- function(config) {
  store <- file.path("results", "local", "study06_ld_operator", "_targets")
  design <- .study06v2_safe_v1_design(store)
  split <- targets::tar_read_raw("study06_split", store = store)
  blocks <- read.csv(file.path(config$local_dir, "deterministic",
    "block_definitions.csv"), stringsAsFactors = FALSE)
  near_full <- readRDS(file.path(config$local_dir, "deterministic",
    "near_full_operator_checkpoint.rds"))
  if (!identical(near_full$source_sha, config$required_sblr_sha) ||
      !identical(near_full$representation, "low_rank") ||
      !identical(near_full$marker_ids,
        design$study06_operator_stats$marker_names))
    stop("Near-full deterministic checkpoint provenance failed.",
      call. = FALSE)
  prefix <- file.path(config$local_dir, "deterministic", "block_csr",
    "runtime_matched")
  meta <- paste0(prefix, ".meta.txt")
  if (!file.exists(meta)) {
    source("studies/06_ld_operator/operators.R", local = TRUE)
    inspect <- list(
      packed_upper_triangle = lapply(near_full$inspect$factor,
        function(Q) .study06_pack_triangle(crossprod(Q))),
      block_size = near_full$inspect$block_size)
    csr <- .study06_write_runtime_csr(inspect, prefix,
      near_full$marker_ids, zero_tolerance = 0)
  } else {
    csr_read <- sblr::sparseLD_read_CSR(prefix, one_based = FALSE)
    csr <- list(prefix = prefix, row_ptr = csr_read$row_ptr,
      col_idx = csr_read$col_idx, values = csr_read$values,
      nnz = length(csr_read$values))
  }
  source("studies/06_ld_operator/operators.R", local = TRUE)
  block_glist <- .study06_runtime_glist(design$study06_ld_bundle$Glist,
    design$study06_operator_stats, csr)
  list(design = design, split = split, blocks = blocks,
    block_csr = csr, block_glist = block_glist)
}

.study06v2_simulation <- function(architecture, replicate, runtime, config) {
  source("studies/06_ld_operator/simulation.R", local = TRUE)
  spec <- list(architecture = architecture, replicate = as.integer(replicate),
    simulation_seed = .study06v2_seed(config, architecture, replicate))
  sim <- .study06_simulate(spec,
    runtime$design$study06_scaled_genotypes$all, config)
  .study06_validate_simulation(sim,
    runtime$design$study06_scaled_genotypes$all, config)
  sim
}

.study06v2_stats <- function(simulation, runtime, config) {
  source("studies/06_ld_operator/pilot.R", local = TRUE)
  .study06_summary_stats(simulation,
    runtime$design$study06_ld_bundle$Glist, runtime$split, config)$stats
}

.study06v2_checkpoint_path <- function(config, phase, method, replicate) {
  file.path(config$local_dir, "fit_checkpoints", phase,
    paste0(method$architecture, "__", replicate, "__",
      method$configuration, ".rds"))
}

.study06v2_validate_checkpoint <- function(run, method, simulation,
                                           input_hash, config) {
  identical(run$status, "ok") &&
    identical(run$source_sha, config$required_sblr_sha) &&
    identical(run$package_version, config$required_sblr_version) &&
    identical(run$architecture, method$architecture) &&
    identical(run$replicate, simulation$replicate) &&
    identical(run$method$configuration, method$configuration) &&
    identical(run$input_hash, input_hash) &&
    length(run$fit$chains) == 4L &&
    identical(sort(vapply(run$fit$chains, `[[`, integer(1),
      "chain_index")), 1:4) &&
    (!identical(method$operator_family, "retained_low_rank") ||
      (identical(run$operator_contract, config$operator_contract) &&
       identical(run$representation, "low_rank") &&
       isTRUE(all.equal(run$eigen_prop, method$eigen_prop))))
}

.study06v2_cached_fit <- function(method, simulation, stats, runtime, config,
                                  phase, recommendations = NULL) {
  controls <- .study06v2_controls(method, config, phase, recommendations)
  if (identical(method$operator_family, "retained_low_rank")) controls <-
    c(controls, list(representation = "low_rank",
      eigen_prop = method$eigen_prop))
  input_hash <- .study06v2_input_hash(simulation, stats, method, controls,
    config)
  path <- .study06v2_checkpoint_path(config, phase, method,
    simulation$replicate)
  if (file.exists(path)) {
    old <- try(readRDS(path), silent = TRUE)
    if (!inherits(old, "try-error") &&
        .study06v2_validate_checkpoint(old, method, simulation,
          input_hash, config)) return(old)
    stop("Stale or mismatched Study 06 v2 checkpoint: ", path,
      call. = FALSE)
  }
  run <- .study06v2_fit(method, simulation, stats,
    runtime$design$study06_ld_bundle$Glist, runtime$split,
    runtime$blocks, runtime$block_glist, config, phase, recommendations)
  run$input_hash <- input_hash
  status <- data.frame(phase = phase, architecture = method$architecture,
    replicate = simulation$replicate,
    configuration = method$configuration,
    method = method$native_method, operator = method$operator_family,
    representation = method$representation,
    eigen_prop = method$eigen_prop, state = run$status,
    start_time = format(run$started_at, tz = "UTC", usetz = TRUE),
    finish_time = format(run$finished_at, tz = "UTC", usetz = TRUE),
    elapsed_seconds = run$runtime, output_path = path,
    validation_status = if (run$status == "ok") "passed" else "failed",
    error_message = run$reason, stringsAsFactors = FALSE)
  .study06v2_write_csv(status, file.path(config$local_dir, "fit_status",
    phase, paste0(method$architecture, "__", simulation$replicate,
      "__", method$configuration, ".csv")))
  if (!identical(run$status, "ok"))
    stop("Study 06 v2 fit failed: ", run$reason, call. = FALSE)
  .study06v2_atomic_rds(run, path)
  run
}
