.study06v2_source_analysis_helpers <- function() {
  for (file in c("chain_extraction.R", "diagnostics.R", "metrics.R"))
    source(file.path("studies", "06_ld_operator", file), local = parent.frame())
}

.study06v2_run_phase_grid <- function(configurations, architectures,
                                      replicates, phase, config,
                                      recommendations = NULL) {
  runtime <- .study06v2_load_runtime_data(config)
  grid <- .study06v2_method_grid(config)
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

.study06v2_pilot <- function(config) {
  runs <- .study06v2_run_phase_grid(config$configurations,
    config$architectures[1L], 1L, "operator-pilot", config)
  .study06v2_source_analysis_helpers()
  runtime <- .study06v2_load_runtime_data(config)
  simulation <- .study06v2_simulation(config$architectures[1L], 1L,
    runtime, config)
  rows <- lapply(runs, function(run) {
    prediction <- .study06_prediction_metrics(run, simulation,
      runtime$design$study06_scaled_genotypes$test, runtime$split)
    marker <- .study06_marker_metrics(run, simulation,
      runtime$design$study06_scaled_genotypes$all)
    evidence <- .study06_sbayesr_evidence(run)
    metric <- function(x, id) x$value[x$metric == id][1L]
    data.frame(architecture = run$architecture,
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
      marker_effect_correlation = metric(marker,
        "marker_effect_correlation"),
      runtime_seconds = run$runtime,
      warnings = run$warnings, stringsAsFactors = FALSE)
  })
  summary <- do.call(rbind, rows)
  reference <- summary[summary$configuration == "block_csr", ]
  full <- summary[summary$configuration == "low_rank_full", ]
  canonical <- summary[summary$configuration == "low_rank_0995", ]
  pass <- nrow(reference) == 1L && nrow(full) == 1L && nrow(canonical) == 1L &&
    all(is.finite(unlist(summary[c("posterior_heritability", "vgs",
      "ves", "vbs", "prediction_correlation",
      "marker_effect_correlation")]))) &&
    abs(full$posterior_heritability - reference$posterior_heritability) <=
      config$pilot_gate$maximum_heritability_difference &&
    abs(full$prediction_correlation - reference$prediction_correlation) <=
      config$pilot_gate$maximum_prediction_correlation_difference &&
    full$marker_effect_correlation >=
      config$pilot_gate$minimum_posterior_effect_correlation
  gate <- data.frame(
    block_csr_minus_full_rank_h2_absolute = abs(
      reference$posterior_heritability - full$posterior_heritability),
    block_csr_minus_full_rank_prediction_absolute = abs(
      reference$prediction_correlation - full$prediction_correlation),
    full_rank_truth_effect_correlation = full$marker_effect_correlation,
    canonical_truth_effect_correlation = canonical$marker_effect_correlation,
    pass = pass, stringsAsFactors = FALSE)
  output <- file.path(config$local_dir, "operator_pilot")
  .study06v2_write_csv(summary, file.path(output,
    "one_replicate_operator_pilot.csv"))
  .study06v2_write_csv(gate, file.path(output, "pilot_gate.csv"))
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

.study06v2_aggregate <- function(config) {
  .study06v2_source_analysis_helpers()
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
      "retained_low_rank")) as.numeric(run$fit$input$eigen_diagnostics$build$immutable_operator_bytes %||% NA_real_) else NA_real_,
    per_chain_state_bytes = if (identical(run$method$operator_family,
      "retained_low_rank")) as.numeric(run$fit$input$eigen_diagnostics$build$per_chain_residual_bytes %||% NA_real_) else NA_real_,
    warnings = run$warnings, stringsAsFactors = FALSE)))
  low_rank_metadata <- do.call(rbind, Filter(Negate(is.null),
    lapply(runs, .study06v2_low_rank_fit_metadata)))
  evidence <- do.call(rbind, lapply(runs, .study06_sbayesr_evidence))
  marker_agreement <- .study06_marker_agreement(runs)
  combined <- rbind(prediction[c("architecture", "replicate",
    "configuration", "method", "metric", "value")],
    marker[c("architecture", "replicate", "configuration", "method",
      "metric", "value")],
    transform(variable[c("architecture", "replicate", "configuration",
      "metric", "value")], method = "", stringsAsFactors = FALSE)[
        c("architecture", "replicate", "configuration", "method",
          "metric", "value")])
  paired <- .study06v2_paired_differences(combined)
  paired_summary <- .study06v2_paired_summary(paired)
  recovery <- .study06_recovery_summary(parameter)
  convergence <- .study06_convergence_validation_summary(diagnostics)
  simulation_summary <- do.call(rbind, lapply(simulations,
    .study06_simulation_summary))
  outputs <- list(seed_registry.csv = .study06v2_seed_registry(config),
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
