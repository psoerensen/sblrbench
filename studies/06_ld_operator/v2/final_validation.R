.study06v2_finite_numeric <- function(x) {
  is.numeric(x) && length(x) > 0L && all(is.finite(as.numeric(x)))
}

.study06v2_validate_scientific_run <- function(run, method, config) {
  fit <- run$fit
  required <- c("bm", "vbs", "vgs", "ves", "vle", "vld")
  if (!all(vapply(required, function(field)
      .study06v2_finite_numeric(fit[[field]]), logical(1))))
    stop("Non-finite required posterior field in ", method$method_id,
      call. = FALSE)
  if (nrow(fit$bm) != 37991L || ncol(fit$bm) != 1L)
    stop("Marker-effect dimensions failed in ", method$method_id,
      call. = FALSE)
  if (!is.null(fit$dm) && (!.study06v2_finite_numeric(fit$dm) ||
      any(fit$dm < -config$operator_tolerance$probability) ||
      any(fit$dm > 1 + config$operator_tolerance$probability)))
    stop("Posterior inclusion probabilities failed in ", method$method_id,
      call. = FALSE)
  if (!is.null(fit$pi_mean) && (!.study06v2_finite_numeric(fit$pi_mean) ||
      any(fit$pi_mean < -config$operator_tolerance$probability) ||
      any(fit$pi_mean > 1 + config$operator_tolerance$probability)))
    stop("Posterior mixture probabilities failed in ", method$method_id,
      call. = FALSE)
  trace_fields <- c("vbs", "vgs", "ves", "vle", "vld")
  for (chain in fit$chains) {
    if (!all(vapply(trace_fields, function(field)
        .study06v2_finite_numeric(chain[[field]]), logical(1))))
      stop("Non-finite chain trace in ", method$method_id, call. = FALSE)
  }
  if (!is.finite(run$runtime) || run$runtime < 0)
    stop("Invalid runtime in ", method$method_id, call. = FALSE)

  low_rank <- identical(method$operator_family, "retained_low_rank")
  residual <- fit$diagnostics$native$low_rank_residual
  if (low_rank) {
    blocks <- fit$input$eigen_diagnostics$blocks
    build <- fit$input$eigen_diagnostics$build
    required_blocks <- c("block_size", "positive_rank", "retained_rank",
      "positive_eigenvalue_mass", "retained_eigenvalue_mass",
      "retained_mass_fraction")
    if (!is.data.frame(blocks) || !all(required_blocks %in% names(blocks)) ||
        nrow(blocks) != 38L || any(!is.finite(as.matrix(
          blocks[required_blocks]))) || any(blocks$retained_rank <= 0L) ||
        any(blocks$retained_rank > blocks$positive_rank) ||
        any(blocks$positive_rank > blocks$block_size) ||
        any(blocks$retained_mass_fraction + 1e-12 < method$eigen_prop))
      stop("Low-rank block metadata failed in ", method$method_id,
        call. = FALSE)
    if (identical(method$configuration, "low_rank_full") &&
        any(blocks$retained_rank != blocks$positive_rank))
      stop("Near-full rank failed in ", method$method_id, call. = FALSE)
    build_values <- unlist(build[c("operator_storage_bytes",
      "chain_residual_storage_bytes", "construction_workspace_bytes",
      "construction_time", "cross_product_time",
      "eigendecomposition_time", "transformation_time")], use.names = FALSE)
    if (!length(build_values) || any(!is.finite(build_values)) ||
        any(build_values < 0))
      stop("Low-rank build diagnostics failed in ", method$method_id,
        call. = FALSE)
    if (is.null(residual) ||
        any(!is.finite(unlist(residual, use.names = FALSE))) ||
        any(residual$low_rank_residual_rebuild_every !=
          config$low_rank_residual_rebuild_every) ||
        any(residual$low_rank_residual_rebuild_count < 1) ||
        any(residual$low_rank_residual_max_abs_drift >
          config$operator_tolerance$product_absolute))
      stop("Low-rank residual drift validation failed in ", method$method_id,
        call. = FALSE)
  } else if (!is.null(residual) || identical(run$representation,
      "dense_reconstructed")) {
    stop("Canonical non-low-rank fit contains forbidden eigen metadata.",
      call. = FALSE)
  }
  invisible(TRUE)
}

.study06v2_final_checkpoint_inventory <- function(config,
                                                   validate_inputs = TRUE) {
  recommendation_path <- file.path(config$local_dir, "convergence",
    "method_recommendations.csv")
  if (!file.exists(recommendation_path))
    stop("Convergence recommendations are unavailable.", call. = FALSE)
  recommendations <- read.csv(recommendation_path, stringsAsFactors = FALSE)
  runtime <- if (validate_inputs) .study06v2_load_runtime_data(config) else NULL
  grid <- .study06v2_method_grid(config)
  expected <- expand.grid(architecture = config$architectures,
    replicate = seq_len(config$replicate_count),
    configuration = config$configurations, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  rows <- vector("list", nrow(expected)); index <- 0L
  for (architecture in config$architectures) {
    for (replicate in seq_len(config$replicate_count)) {
      simulation <- if (validate_inputs) .study06v2_simulation(architecture,
        replicate, runtime, config) else NULL
      stats <- if (validate_inputs) .study06v2_stats(simulation, runtime,
        config) else NULL
      methods <- grid[grid$architecture == architecture, , drop = FALSE]
      for (i in seq_len(nrow(methods))) {
        index <- index + 1L
        method <- as.list(methods[i, , drop = FALSE])
        path <- .study06v2_checkpoint_path(config, "benchmark", method,
          replicate)
        if (!file.exists(path) || grepl("(^|[/\\])\\.", path))
          stop("Missing or temporary benchmark checkpoint: ", path,
            call. = FALSE)
        run <- try(readRDS(path), silent = TRUE)
        if (inherits(run, "try-error"))
          stop("Unreadable benchmark checkpoint: ", path, call. = FALSE)
        controls <- .study06v2_controls(method, config, "benchmark",
          recommendations)
        if (identical(method$operator_family, "retained_low_rank"))
          controls <- c(controls, list(representation = "low_rank",
            eigen_prop = method$eigen_prop,
            low_rank_residual_rebuild_every =
              config$low_rank_residual_rebuild_every))
        input_hash <- if (validate_inputs) .study06v2_input_hash(simulation,
          stats, method, controls, config) else run$input_hash
        if (!.study06v2_validate_checkpoint(run, method,
            if (validate_inputs) simulation else list(replicate = replicate),
            input_hash, config))
          stop("Checkpoint provenance/schema failed: ", path, call. = FALSE)
        expected_controls <- .study06v2_controls(method, config, "benchmark",
          recommendations)
        fields <- c("nit", "nburn", "nthin", "nchains", "ncores")
        if (!identical(unname(unlist(run$controls[fields])),
            unname(unlist(expected_controls[fields]))))
          stop("MCMC controls differ from recommendations: ", path,
            call. = FALSE)
        .study06v2_validate_scientific_run(run, method, config)
        low <- identical(method$operator_family, "retained_low_rank")
        blocks <- if (low) run$fit$input$eigen_diagnostics$blocks else NULL
        build <- if (low) run$fit$input$eigen_diagnostics$build else list()
        residual <- if (low)
          run$fit$diagnostics$native$low_rank_residual else list()
        rows[[index]] <- data.frame(architecture = architecture,
          replicate = replicate, configuration = method$configuration,
          package_version = run$package_version,
          package_commit = run$source_sha, chain_count = length(run$fit$chains),
          nburn = run$controls$nburn, retained_iterations = run$controls$nit,
          total_iterations = run$controls$nburn + run$controls$nit,
          nthin = run$controls$nthin, ncores = run$controls$ncores,
          operator_representation = if (low) run$representation else
            method$operator_family,
          eigen_prop = if (low) run$eigen_prop else NA_real_,
          positive_rank = if (low) sum(blocks$positive_rank) else NA_real_,
          retained_rank = if (low) sum(blocks$retained_rank) else NA_real_,
          rank_fraction = if (low) sum(blocks$retained_rank) /
            sum(blocks$block_size) else NA_real_,
          retained_mass_fraction = if (low)
            sum(blocks$retained_eigenvalue_mass) /
              sum(blocks$positive_eigenvalue_mass) else NA_real_,
          residual_rebuild_interval = if (low)
            max(residual$low_rank_residual_rebuild_every) else NA_real_,
          residual_rebuild_count = if (low)
            max(residual$low_rank_residual_rebuild_count) else NA_real_,
          maximum_residual_drift = if (low)
            max(residual$low_rank_residual_max_abs_drift) else NA_real_,
          construction_seconds = as.numeric(build$construction_time %||%
            NA_real_), eigendecomposition_seconds = as.numeric(
              build$eigendecomposition_time %||% NA_real_),
          transformation_seconds = as.numeric(build$transformation_time %||%
            NA_real_), mcmc_wall_seconds = run$runtime,
          operator_storage_bytes = as.numeric(build$operator_storage_bytes %||%
            NA_real_), fit_object_bytes = as.numeric(object.size(run$fit)),
          checkpoint_bytes = file.info(path)$size,
          input_hash = run$input_hash, status = run$status,
          warning = run$warnings, stringsAsFactors = FALSE)
        rm(run); gc(FALSE)
      }
    }
  }
  inventory <- do.call(rbind, rows)
  key <- interaction(inventory$architecture, inventory$replicate,
    inventory$configuration, drop = TRUE)
  if (nrow(inventory) != config$expected_fit_count || anyDuplicated(key) ||
      !setequal(do.call(paste, expected), do.call(paste,
        inventory[c("architecture", "replicate", "configuration")])))
    stop("Final checkpoint coordinate inventory failed.", call. = FALSE)
  inventory <- inventory[order(inventory$architecture,
    inventory$replicate, match(inventory$configuration,
      config$configurations)), , drop = FALSE]
  rownames(inventory) <- NULL
  inventory
}

.study06v2_write_final_validation <- function(config) {
  inventory <- .study06v2_final_checkpoint_inventory(config,
    validate_inputs = TRUE)
  output <- file.path(config$local_dir, "verification")
  .study06v2_write_csv(inventory, file.path(output,
    "checkpoint_validation.csv"))
  summary <- data.frame(expected_fit_count = config$expected_fit_count,
    observed_fit_count = nrow(inventory), failed_fit_count = sum(
      inventory$status != "ok"), duplicate_coordinate_count = 0L,
    maximum_residual_drift = max(inventory$maximum_residual_drift,
      na.rm = TRUE), validation_status = "passed", stringsAsFactors = FALSE)
  .study06v2_write_csv(summary, file.path(output,
    "checkpoint_validation_summary.csv"))
  invisible(inventory)
}
