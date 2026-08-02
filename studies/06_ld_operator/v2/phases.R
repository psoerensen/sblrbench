.study06v2_source_analysis_helpers <- function() {
  target <- environment(.study06v2_source_analysis_helpers)
  for (file in c("chain_extraction.R", "diagnostics.R", "metrics.R",
      "simulation.R"))
    source(file.path("studies", "06_ld_operator", file), local = target)
}

.study06v2_run_phase_grid <- function(configurations, architectures,
                                      replicates, phase, config,
                                      recommendations = NULL) {
  runtime <- .study06v2_load_runtime_data(config)
  grid <- .study06v2_method_grid(config, configurations)
  grid <- grid[grid$configuration %in% configurations &
    grid$architecture %in% architectures, , drop = FALSE]
  runs <- list()
  for (architecture in architectures) for (replicate in replicates) {
    simulation <- .study06v2_simulation(architecture, replicate, runtime,
      config)
    stats <- .study06v2_stats(simulation, runtime, config)
    methods <- grid[grid$architecture == architecture, , drop = FALSE]
    for (i in seq_len(nrow(methods))) {
      method <- as.list(methods[i, , drop = FALSE])
      message(sprintf("Study 06 v2 %s: %s / replicate %d / %s",
        phase, architecture, replicate, method$configuration))
      runs[[length(runs) + 1L]] <- .study06v2_cached_fit(method,
        simulation, stats, runtime, config, phase, recommendations)
      gc()
    }
  }
  runs
}

.study06v2_low_rank_residual <- function(run) {
  if (!identical(run$method$operator_family, "retained_low_rank"))
    return(data.frame(rebuild_every = NA_integer_, rebuild_count = NA_real_,
      maximum_drift = NA_real_, available = FALSE))
  native <- run$fit$diagnostics$native$low_rank_residual
  if (is.null(native)) stop("Low-rank residual diagnostics are absent.",
    call. = FALSE)
  data.frame(
    rebuild_every = as.integer(max(native$low_rank_residual_rebuild_every)),
    rebuild_count = as.numeric(max(native$low_rank_residual_rebuild_count)),
    maximum_drift = as.numeric(max(native$low_rank_residual_max_abs_drift)),
    available = TRUE)
}

.study06v2_native_timing <- function(run) {
  native <- run$fit$diagnostics$native
  reported <- as.numeric(native$seconds_mean %||% NA_real_)
  construction <- as.numeric(run$fit$input$eigen_diagnostics$build$
    construction_time %||% 0)
  usable <- length(reported) == 1L && is.finite(reported) && reported > 0
  chain_elapsed <- if (usable) reported else
    max(as.numeric(run$runtime) - construction, 0)
  data.frame(
    chain_elapsed_seconds = chain_elapsed,
    chain_timing_source = if (usable) "native_seconds_mean" else
      "wall_minus_construction_estimate",
    chain_seconds_max = as.numeric(native$seconds_max %||% NA_real_),
    seconds_per_chain_iteration = chain_elapsed /
      (run$controls$nit + run$controls$nburn))
}

.study06v2_build_timing <- function(run) {
  build <- run$fit$input$eigen_diagnostics$build %||% list()
  data.frame(operator_construction_seconds = as.numeric(
      build$construction_time %||% NA_real_),
    cross_product_seconds = as.numeric(build$cross_product_time %||% NA_real_),
    eigendecomposition_seconds = as.numeric(
      build$eigendecomposition_time %||% NA_real_),
    transformation_seconds = as.numeric(build$transformation_time %||%
      NA_real_),
    operator_storage_bytes = as.numeric(build$operator_storage_bytes %||%
      NA_real_),
    chain_residual_storage_bytes = as.numeric(
      build$chain_residual_storage_bytes %||% NA_real_))
}

.study06v2_preoptimization_pilot <- function(config) {
  path <- file.path(config$preoptimization_archive, "operator_pilot",
    "one_replicate_operator_pilot.csv")
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE)
}

.study06v2_historical_runtime_context <- function(config) {
  v1 <- read.csv(file.path(config$historical_capsules[1L],
    "computational_summary.csv"), stringsAsFactors = FALSE)
  v1$evidence_source <- "historical_v1_maximum_history"
  old_paths <- list.files(file.path(config$preoptimization_archive,
    "fit_checkpoints", "convergence"), pattern = "[.]rds$",
    full.names = TRUE)
  old <- do.call(rbind, lapply(old_paths, function(path) {
    x <- readRDS(path)
    data.frame(architecture = x$architecture, replicate = x$replicate,
      configuration = x$method$configuration,
      method = x$method$native_method, status = x$status,
      error_message = x$reason, chain_count = length(x$fit$chains),
      elapsed_seconds = x$runtime, warnings = x$warnings,
      evidence_source = "interrupted_preoptimization_v2_convergence",
      stringsAsFactors = FALSE)
  }))
  common <- Reduce(intersect, list(names(v1), names(old)))
  rbind(v1[common], old[common])
}

.study06v2_pilot <- function(config) {
  runs <- .study06v2_run_phase_grid(config$pilot_configurations,
    config$architectures, 1L, "operator-pilot", config)
  .study06v2_source_analysis_helpers()
  runtime <- .study06v2_load_runtime_data(config)
  simulations <- setNames(lapply(config$architectures, function(a)
    .study06v2_simulation(a, 1L, runtime, config)), config$architectures)
  rows <- lapply(runs, function(run) {
    simulation <- simulations[[run$architecture]]
    prediction <- .study06_prediction_metrics(run, simulation,
      runtime$design$study06_scaled_genotypes$test, runtime$split)
    marker <- .study06_marker_metrics(run, simulation,
      runtime$design$study06_scaled_genotypes$all)
    evidence <- .study06_sbayesr_evidence(run)
    metric <- function(x, id) x$value[x$metric == id][1L]
    residual <- .study06v2_low_rank_residual(run)
    timing <- .study06v2_native_timing(run)
    build <- .study06v2_build_timing(run)
    cbind(data.frame(architecture = run$architecture,
      replicate = run$replicate,
      configuration = run$method$configuration,
      posterior_heritability = evidence$benchmark_heritability,
      vgs = evidence$vgs_posterior_mean,
      ves = evidence$ves_posterior_mean,
      vbs = evidence$vbs_posterior_mean,
      vld = evidence$vld_posterior_mean,
      vle = evidence$vle_posterior_mean,
      pis = evidence$pis_posterior_mean,
      prediction_correlation = metric(prediction,
        "phenotype_prediction_correlation"),
      prediction_rmse = metric(prediction, "phenotype_prediction_rmse"),
      prediction_slope = metric(prediction, "prediction_regression_slope"),
      genetic_value_correlation = metric(prediction,
        "genetic_value_correlation"),
      genetic_value_rmse = metric(prediction, "genetic_value_rmse"),
      marker_effect_correlation = metric(marker,
        "marker_effect_correlation"),
      marker_effect_rmse = metric(marker, "marker_effect_rmse"),
      runtime_seconds = run$runtime,
      fit_object_bytes = as.numeric(object.size(run$fit)),
      representation = run$representation,
      eigen_prop = run$eigen_prop,
      operator_contract = run$operator_contract,
      source_sha = run$source_sha,
      package_version = run$package_version,
      warnings = run$warnings, stringsAsFactors = FALSE), timing, build,
      residual)
  })
  summary <- do.call(rbind, rows)
  comparisons <- list(); gate_rows <- list()
  pairs <- list(c("full_csr", "bed"), c("block_csr", "full_csr"),
    c("low_rank_full", "block_csr"),
    c("low_rank_0999", "low_rank_full"),
    c("low_rank_0995", "low_rank_full"),
    c("dense_reconstructed_unfiltered", "low_rank_full"))
  for (architecture in config$architectures) for (pair in pairs) {
    focal <- Filter(function(x) x$architecture == architecture &&
      x$method$configuration == pair[1L], runs)[[1L]]
    reference <- Filter(function(x) x$architecture == architecture &&
      x$method$configuration == pair[2L], runs)[[1L]]
    sf <- summary[summary$architecture == architecture &
      summary$configuration == pair[1L], ]
    sr <- summary[summary$architecture == architecture &
      summary$configuration == pair[2L], ]
    comparisons[[length(comparisons) + 1L]] <- data.frame(
      architecture = architecture,
      comparison = paste(pair[1L], "minus", pair[2L]),
      heritability_difference = sf$posterior_heritability -
        sr$posterior_heritability,
      prediction_correlation_difference = sf$prediction_correlation -
        sr$prediction_correlation,
      posterior_effect_correlation = .study06_safe_cor(
        as.numeric(focal$fit$bm), as.numeric(reference$fit$bm)),
      posterior_effect_rmse = sqrt(mean((as.numeric(focal$fit$bm) -
        as.numeric(reference$fit$bm))^2)),
      pip_correlation = if (is.null(focal$fit$dm) ||
        is.null(reference$fit$dm)) NA_real_ else .study06_safe_cor(
          as.numeric(focal$fit$dm), as.numeric(reference$fit$dm)),
      stringsAsFactors = FALSE)
  }
  comparisons <- do.call(rbind, comparisons)
  low_rank <- summary[startsWith(summary$configuration, "low_rank_"), ]
  science_pairs <- comparisons[comparisons$comparison %in% c(
    "low_rank_full minus block_csr",
    "low_rank_0999 minus low_rank_full",
    "low_rank_0995 minus low_rank_full",
    "dense_reconstructed_unfiltered minus low_rank_full"), ]
  old <- .study06v2_preoptimization_pilot(config)
  old_low_rank <- old[old$configuration %in% low_rank$configuration, ]
  timing <- merge(low_rank[c("architecture", "configuration",
    "runtime_seconds", "seconds_per_chain_iteration")],
    old_low_rank[c("architecture", "configuration", "runtime_seconds")],
    by = c("architecture", "configuration"), all.x = TRUE,
    suffixes = c("_optimized", "_preoptimization"))
  timing$runtime_ratio <- timing$runtime_seconds_optimized /
    timing$runtime_seconds_preoptimization
  comparable_timing <- timing[is.finite(timing$runtime_ratio), , drop = FALSE]
  canonical <- summary[summary$configuration %in% config$configurations, ]
  canonical$target_total_iterations <- vapply(seq_len(nrow(canonical)),
    function(i) {
      method <- as.list(.study06v2_method_for(canonical$architecture[i],
        canonical$configuration[i], config))
      controls <- .study06v2_baseline_controls(method, config)
      controls$nit + controls$nburn
    }, numeric(1))
  construction <- ifelse(is.finite(canonical$operator_construction_seconds),
    canonical$operator_construction_seconds, 0)
  canonical$projected_fit_seconds <- construction +
    canonical$seconds_per_chain_iteration * canonical$target_total_iterations
  projected_grid_hours <- sum(canonical$projected_fit_seconds) * 5 / 3600
  pass <- nrow(summary) == length(config$architectures) *
      length(config$pilot_configurations) &&
    all(summary$source_sha == config$required_sblr_sha) &&
    all(is.finite(unlist(summary[c("posterior_heritability", "vgs",
      "ves", "vbs", "prediction_correlation",
      "marker_effect_correlation")]))) &&
    all(low_rank$rebuild_every == config$low_rank_residual_rebuild_every) &&
    all(is.finite(low_rank$maximum_drift)) &&
    all(low_rank$maximum_drift <=
      config$optimized_pilot_gate$maximum_residual_drift) &&
    all(abs(science_pairs$heritability_difference) <=
      config$pilot_gate$maximum_heritability_difference) &&
    all(abs(science_pairs$prediction_correlation_difference) <=
      config$pilot_gate$maximum_prediction_correlation_difference) &&
    all(science_pairs$posterior_effect_correlation >=
      config$pilot_gate$minimum_posterior_effect_correlation) &&
    nrow(comparable_timing) > 0L &&
    all(comparable_timing$runtime_ratio <=
      config$optimized_pilot_gate$maximum_preoptimization_runtime_ratio) &&
    projected_grid_hours <=
      config$optimized_pilot_gate$maximum_projected_grid_hours
  gate <- data.frame(expected_fit_count = length(config$architectures) *
      length(config$pilot_configurations), completed_fit_count = nrow(summary),
    maximum_low_rank_residual_drift = max(low_rank$maximum_drift),
    maximum_absolute_heritability_difference = max(abs(
      science_pairs$heritability_difference)),
    minimum_posterior_effect_correlation = min(
      science_pairs$posterior_effect_correlation),
    comparable_preoptimization_timing_count = nrow(comparable_timing),
    maximum_preoptimization_runtime_ratio = max(
      comparable_timing$runtime_ratio),
    projected_complete_grid_hours = projected_grid_hours,
    pass = pass, stringsAsFactors = FALSE)
  output <- file.path(config$local_dir, "operator_pilot")
  .study06v2_write_csv(summary, file.path(output,
    "one_replicate_operator_pilot.csv"))
  .study06v2_write_csv(gate, file.path(output, "pilot_gate.csv"))
  .study06v2_write_csv(comparisons, file.path(output,
    "pilot_paired_comparisons.csv"))
  .study06v2_write_csv(timing, file.path(output,
    "pilot_runtime_comparison.csv"))
  .study06v2_write_csv(canonical[c("architecture", "configuration",
    "target_total_iterations", "projected_fit_seconds")], file.path(output,
      "projected_complete_grid_runtime.csv"))
  .study06v2_write_csv(.study06v2_historical_runtime_context(config),
    file.path(output, "historical_runtime_context.csv"))
  if (!isTRUE(pass)) stop("Study 06 v2 operator pilot gate failed.",
    call. = FALSE)
  invisible(TRUE)
}

.study06v2_convergence <- function(config) {
  configurations <- setdiff(config$configurations, "bed")
  runs <- .study06v2_run_phase_grid(configurations,
    config$architectures, 1L, "convergence", config)
  .study06v2_source_analysis_helpers()
  selections <- lapply(runs, function(run)
    .study06_select_recommendation(.study06_extract_draws(run), config))
  recommendations <- do.call(rbind, lapply(selections, `[[`,
    "recommendation"))
  candidates <- do.call(rbind, lapply(selections, `[[`, "candidates"))
  diagnostics <- do.call(rbind, lapply(selections, `[[`, "diagnostics"))
  if (nrow(recommendations) != config$convergence_fit_count ||
      any(recommendations$recommendation_status != "available") ||
      !all(config$architectures %in%
        recommendations$architecture[recommendations$configuration ==
          "low_rank_0995"]))
    stop("Canonical low-rank convergence recommendation is unsupported.",
      call. = FALSE)
  status <- do.call(rbind, lapply(runs, function(run) data.frame(
    architecture = run$architecture, replicate = run$replicate,
    configuration = run$method$configuration,
    status = run$status, chain_count = length(run$fit$chains),
    elapsed_seconds = run$runtime, warnings = run$warnings,
    error_message = run$reason, stringsAsFactors = FALSE)))
  output <- file.path(config$local_dir, "convergence")
  .study06v2_write_csv(recommendations,
    file.path(output, "method_recommendations.csv"))
  .study06v2_write_csv(candidates,
    file.path(output, "candidate_settings.csv"))
  .study06v2_write_csv(diagnostics,
    file.path(output, "convergence_diagnostics.csv"))
  .study06v2_write_csv(status, file.path(output, "fit_status.csv"))
  invisible(TRUE)
}

.study06v2_benchmark <- function(config) {
  path <- file.path(config$local_dir, "convergence",
    "method_recommendations.csv")
  if (!file.exists(path))
    stop("Validated v2 convergence results are required.", call. = FALSE)
  recommendations <- read.csv(path, stringsAsFactors = FALSE)
  runs <- .study06v2_run_phase_grid(config$configurations,
    config$architectures, seq_len(config$replicate_count), "benchmark",
    config, recommendations)
  if (length(runs) != config$expected_fit_count)
    stop("Study 06 v2 benchmark fit count is incomplete.", call. = FALSE)
  invisible(TRUE)
}

.study06v2_read_benchmark_runs <- function(config) {
  grid <- .study06v2_method_grid(config)
  out <- list()
  for (replicate in seq_len(config$replicate_count))
    for (i in seq_len(nrow(grid))) {
      method <- as.list(grid[i, , drop = FALSE])
      path <- .study06v2_checkpoint_path(config, "benchmark", method,
        replicate)
      if (!file.exists(path)) stop("Missing benchmark checkpoint: ", path,
        call. = FALSE)
      out[[length(out) + 1L]] <- readRDS(path)
    }
  if (length(out) != config$expected_fit_count ||
      any(vapply(out, `[[`, "", "status") != "ok"))
    stop("Benchmark checkpoint grid is incomplete.", call. = FALSE)
  out
}

.study06v2_comparison_registry <- function() data.frame(
  comparison_id = c("full_csr_minus_bed", "block_csr_minus_full_csr",
    "low_rank_full_minus_block_csr", "low_rank_0995_minus_low_rank_full",
    "low_rank_0995_minus_low_rank_0999",
    "low_rank_0999_minus_low_rank_full",
    "low_rank_0995_minus_full_csr"),
  focal = c("full_csr", "block_csr", "low_rank_full", "low_rank_0995",
    "low_rank_0995", "low_rank_0999", "low_rank_0995"),
  reference = c("bed", "full_csr", "block_csr", "low_rank_full",
    "low_rank_0999", "low_rank_full", "full_csr"),
  interpretation_level = c(4L, 2L, 1L, 3L, 3L, 3L, 4L),
  stringsAsFactors = FALSE)

.study06v2_paired_differences <- function(table) {
  registry <- .study06v2_comparison_registry(); out <- list()
  for (i in seq_len(nrow(registry))) for (architecture in unique(
      table$architecture)) for (metric in unique(table$metric)) {
    a <- table[table$architecture == architecture &
      table$configuration == registry$focal[i] & table$metric == metric,
      c("replicate", "value"), drop = FALSE]
    b <- table[table$architecture == architecture &
      table$configuration == registry$reference[i] & table$metric == metric,
      c("replicate", "value"), drop = FALSE]
    names(a)[2L] <- "focal_value"; names(b)[2L] <- "reference_value"
    z <- merge(a, b, by = "replicate", all = TRUE)
    if (!nrow(z)) next
    z$complete_pair <- is.finite(z$focal_value) & is.finite(z$reference_value)
    z$difference <- z$focal_value - z$reference_value
    z$architecture <- architecture; z$metric <- metric
    z$comparison_id <- registry$comparison_id[i]
    z$orientation <- paste(registry$focal[i], "minus", registry$reference[i])
    z$interpretation_level <- registry$interpretation_level[i]
    out[[length(out) + 1L]] <- z
  }
  do.call(rbind, out)
}

.study06v2_paired_summary <- function(x) do.call(rbind, lapply(split(x,
  interaction(x$comparison_id, x$architecture, x$metric, drop = TRUE)),
  function(z) {
    y <- z$difference[z$complete_pair]
    data.frame(comparison_id = z$comparison_id[1L],
      architecture = z$architecture[1L], metric = z$metric[1L],
      orientation = z$orientation[1L], complete_paired_replicates = length(y),
      mean_difference = if (length(y)) mean(y) else NA_real_,
      sd_difference = if (length(y) > 1L) sd(y) else NA_real_,
      median_difference = if (length(y)) median(y) else NA_real_,
      minimum_difference = if (length(y)) min(y) else NA_real_,
      maximum_difference = if (length(y)) max(y) else NA_real_,
      stringsAsFactors = FALSE)
  }))

.study06v2_seed_registry <- function(config) {
  grid <- expand.grid(architecture = config$architectures,
    replicate = seq_len(config$replicate_count),
    configuration = config$configurations, chain = 1:4,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$data_selection_seed <- config$seeds$data_selection
  grid$simulation_seed <- mapply(.study06v2_seed,
    architecture = grid$architecture, replicate = grid$replicate,
    MoreArgs = list(config = config))
  grid$fit_seed <- mapply(.study06v2_seed,
    architecture = grid$architecture, replicate = grid$replicate,
    configuration = grid$configuration, MoreArgs = list(config = config))
  grid$chain_seed <- mapply(.study06v2_seed,
    architecture = grid$architecture, replicate = grid$replicate,
    configuration = grid$configuration, chain = grid$chain,
    MoreArgs = list(config = config))
  grid
}

.study06v2_variable_selection_metrics <- function(run, simulation) {
  pip <- if (is.null(run$fit$dm)) NULL else as.numeric(run$fit$dm[, 1L])
  truth <- as.integer(as.numeric(simulation$effects[, 1L]) != 0)
  if (is.null(pip) || length(pip) != length(truth) ||
      any(!is.finite(pip))) return(data.frame(
        architecture = run$architecture, replicate = run$replicate,
        configuration = run$method$configuration,
        metric = c("pip_causal_mean", "pip_null_mean", "pip_auprc",
          "pip_auroc", "top_marker_causal_enrichment", "total_pip_mass"),
        value = NA_real_, available = FALSE,
        unavailable_reason = "posterior non-null probabilities unavailable",
        stringsAsFactors = FALSE))
  ord <- order(pip, decreasing = TRUE)
  tp <- cumsum(truth[ord]); fp <- cumsum(1L - truth[ord])
  recall <- tp / sum(truth); precision <- tp / pmax(tp + fp, 1)
  auprc <- sum(diff(c(0, recall)) * precision)
  ranks <- rank(pip, ties.method = "average")
  n1 <- sum(truth); n0 <- length(truth) - n1
  auroc <- (sum(ranks[truth == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  top_n <- length(simulation$causal_index)
  enrichment <- mean(truth[ord[seq_len(top_n)]]) / mean(truth)
  values <- c(pip_causal_mean = mean(pip[truth == 1L]),
    pip_null_mean = mean(pip[truth == 0L]), pip_auprc = auprc,
    pip_auroc = auroc, top_marker_causal_enrichment = enrichment,
    total_pip_mass = sum(pip))
  data.frame(architecture = run$architecture, replicate = run$replicate,
    configuration = run$method$configuration, metric = names(values),
    value = unname(values), available = TRUE, unavailable_reason = "",
    stringsAsFactors = FALSE)
}

.study06v2_low_rank_fit_metadata <- function(run) {
  if (!identical(run$method$operator_family, "retained_low_rank"))
    return(NULL)
  d <- run$fit$input$eigen_diagnostics
  blocks <- d$blocks
  if (!is.data.frame(blocks))
    stop("Low-rank fit block diagnostics are absent.", call. = FALSE)
  cbind(data.frame(architecture = run$architecture,
    replicate = run$replicate,
    configuration = run$method$configuration,
    source_sha = run$source_sha,
    package_version = run$package_version,
    operator_contract = run$fit$input$operator_contract,
    representation = run$fit$input$operator_representation,
    eigen_prop = run$fit$input$eigen_prop,
    stringsAsFactors = FALSE), blocks)
}

.study06v2_marker_agreement <- function(runs) {
  registry <- .study06v2_comparison_registry()
  out <- list()
  architectures <- unique(vapply(runs, `[[`, "", "architecture"))
  replicates <- sort(unique(vapply(runs, `[[`, integer(1), "replicate")))
  for (architecture in architectures) for (replicate in replicates)
    for (i in seq_len(nrow(registry))) {
      focal <- Filter(function(x) x$architecture == architecture &&
        x$replicate == replicate && x$method$configuration ==
          registry$focal[i], runs)
      reference <- Filter(function(x) x$architecture == architecture &&
        x$replicate == replicate && x$method$configuration ==
          registry$reference[i], runs)
      if (length(focal) != 1L || length(reference) != 1L) next
      a <- as.numeric(focal[[1L]]$fit$bm[, 1L])
      b <- as.numeric(reference[[1L]]$fit$bm[, 1L])
      da <- if (is.null(focal[[1L]]$fit$dm)) NULL else
        as.numeric(focal[[1L]]$fit$dm[, 1L])
      db <- if (is.null(reference[[1L]]$fit$dm)) NULL else
        as.numeric(reference[[1L]]$fit$dm[, 1L])
      out[[length(out) + 1L]] <- data.frame(
        architecture = architecture, replicate = replicate,
        comparison_id = registry$comparison_id[i],
        orientation = paste(registry$focal[i], "minus",
          registry$reference[i]),
        posterior_mean_effect_correlation = .study06_safe_cor(a, b),
        posterior_mean_effect_rmse = sqrt(mean((a - b)^2)),
        posterior_mean_effect_maximum_absolute_difference =
          max(abs(a - b)),
        nonnull_probability_correlation =
          if (is.null(da) || is.null(db)) NA_real_ else
            .study06_safe_cor(da, db),
        high_pip_overlap = if (is.null(da) || is.null(db)) NA_real_ else {
          A <- which(da >= .5); B <- which(db >= .5)
          if (!length(union(A, B))) 1 else
            length(intersect(A, B)) / length(union(A, B))
        }, stringsAsFactors = FALSE)
    }
  if (!length(out)) return(data.frame())
  z <- do.call(rbind, out)
  z[order(z$comparison_id, z$architecture, z$replicate), , drop = FALSE]
}

.study06v2_output_availability <- function(runs) do.call(rbind,
  lapply(runs, function(run) {
    fields <- c("vgs", "ves", "vbs", "vld", "vle", "pis", "dm",
      "component_probabilities", "log_cpo", "mean_log_cpo")
    data.frame(architecture = run$architecture, replicate = run$replicate,
      configuration = run$method$configuration, output = fields,
      available = vapply(fields, function(field)
        !is.null(run$fit[[field]]), logical(1)),
      unavailable_reason = vapply(fields, function(field)
        if (is.null(run$fit[[field]])) "not returned by public fit contract"
        else "", ""), stringsAsFactors = FALSE)
  }))

.study06v2_aggregate <- function(config) {
  .study06v2_source_analysis_helpers()
  final_inventory <- .study06v2_write_final_validation(config)
  runtime <- .study06v2_load_runtime_data(config)
  runs <- .study06v2_read_benchmark_runs(config)
  simulations <- list()
  for (architecture in config$architectures)
    for (replicate in seq_len(config$replicate_count))
      simulations[[paste(architecture, replicate)]] <-
        .study06v2_simulation(architecture, replicate, runtime, config)
  find_sim <- function(run) simulations[[paste(run$architecture,
    run$replicate)]]
  draws <- lapply(runs, .study06_extract_draws)
  prediction <- do.call(rbind, Map(function(run, sim)
    .study06_prediction_metrics(run, sim,
      runtime$design$study06_scaled_genotypes$test, runtime$split),
    runs, lapply(runs, find_sim)))
  parameter <- do.call(rbind, Map(.study06_parameter_estimates, draws,
    lapply(runs, find_sim)))
  marker <- do.call(rbind, Map(function(run, sim)
    .study06_marker_metrics(run, sim,
      runtime$design$study06_scaled_genotypes$all), runs,
    lapply(runs, find_sim)))
  variable <- do.call(rbind, Map(.study06v2_variable_selection_metrics,
    runs, lapply(runs, find_sim)))
  diagnostics <- do.call(rbind, Map(function(x, run)
    .study06_diagnostics(x, 0L, run$controls$nit, config), draws, runs))
  status <- do.call(rbind, lapply(runs, function(run) data.frame(
    architecture = run$architecture, replicate = run$replicate,
    configuration = run$method$configuration,
    method = run$method$native_method, status = run$status,
    chain_count = length(run$fit$chains), input_hash = run$input_hash,
    source_sha = run$source_sha, package_version = run$package_version,
    operator_contract = run$operator_contract,
    representation = run$representation, eigen_prop = run$eigen_prop,
    warnings = run$warnings, error_message = run$reason,
    stringsAsFactors = FALSE)))
  computational <- do.call(rbind, lapply(runs, function(run) data.frame(
    architecture = run$architecture, replicate = run$replicate,
    configuration = run$method$configuration,
    method = run$method$native_method,
    mcmc_runtime_seconds = run$runtime,
    seconds_per_retained_draw = run$runtime /
      (run$controls$nit * run$controls$nchains),
    fit_object_bytes = as.numeric(object.size(run$fit)),
    checkpoint_bytes = file.info(.study06v2_checkpoint_path(config,
      "benchmark", run$method, run$replicate))$size,
    immutable_operator_bytes = if (identical(run$method$operator_family,
      "retained_low_rank")) as.numeric(run$fit$input$eigen_diagnostics$build$operator_storage_bytes %||% NA_real_) else NA_real_,
    per_chain_state_bytes = if (identical(run$method$operator_family,
      "retained_low_rank")) as.numeric(run$fit$input$eigen_diagnostics$build$chain_residual_storage_bytes %||% NA_real_) else NA_real_,
    operator_construction_seconds = if (identical(
      run$method$operator_family, "retained_low_rank"))
        as.numeric(run$fit$input$eigen_diagnostics$build$construction_time %||%
          NA_real_) else NA_real_,
    warnings = run$warnings, stringsAsFactors = FALSE)))
  low_rank_metadata <- do.call(rbind, Filter(Negate(is.null),
    lapply(runs, .study06v2_low_rank_fit_metadata)))
  evidence <- do.call(rbind, lapply(runs, .study06_sbayesr_evidence))
  marker_agreement <- .study06v2_marker_agreement(runs)
  availability <- .study06v2_output_availability(runs)
  variable_combined <- variable[c("architecture", "replicate",
    "configuration", "metric", "value")]
  variable_combined$method <- ""
  combined <- rbind(prediction[c("architecture", "replicate",
    "configuration", "method", "metric", "value")],
    marker[c("architecture", "replicate", "configuration", "method",
      "metric", "value")],
    variable_combined[c("architecture", "replicate", "configuration",
      "method", "metric", "value")])
  paired <- .study06v2_paired_differences(combined)
  paired_summary <- .study06v2_paired_summary(paired)
  recovery <- .study06_recovery_summary(parameter)
  convergence <- .study06_convergence_validation_summary(diagnostics)
  simulation_summary <- do.call(rbind, lapply(simulations,
    .study06_simulation_summary))
  outputs <- list(seed_registry.csv = .study06v2_seed_registry(config),
    checkpoint_validation.csv = final_inventory,
    simulation_summary.csv = simulation_summary,
    block_definitions.csv = runtime$blocks,
    fit_status.csv = status, prediction_metrics.csv = prediction,
    parameter_estimates.csv = parameter,
    parameter_recovery_summary.csv = recovery,
    marker_effect_metrics.csv = marker,
    variable_selection_metrics.csv = variable,
    marker_effect_agreement.csv = marker_agreement,
    sbayesr_diagnostic_evidence.csv = evidence,
    paired_replicate_differences.csv = paired,
    paired_comparison_summary.csv = paired_summary,
    convergence_diagnostics.csv = diagnostics,
    convergence_validation_summary.csv = convergence,
    computational_summary.csv = computational,
    method_output_availability.csv = availability,
    low_rank_fit_block_diagnostics.csv = low_rank_metadata)
  output <- file.path(config$local_dir, "aggregate")
  vapply(names(outputs), function(name)
    .study06v2_write_csv(outputs[[name]], file.path(output, name)), "")
  if (nrow(status) != config$expected_fit_count ||
      any(status$status != "ok") || any(status$chain_count != 4L) ||
      any(paired_summary$complete_paired_replicates != 5L))
    stop("Study 06 v2 aggregate completeness validation failed.",
      call. = FALSE)
  invisible(TRUE)
}
