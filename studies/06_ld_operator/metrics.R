.study06_safe_cor <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2L ||
      any(!is.finite(x)) || any(!is.finite(y)) ||
      stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}

.study06_prediction_metrics <- function(run, simulation, Z_test, split) {
  if (!identical(run$status, "ok"))
    stop("Cannot score failed Study 06 fit.", call. = FALSE)
  effects <- as.numeric(run$fit$bm[, 1L])
  if (length(effects) != ncol(Z_test) || any(!is.finite(effects)))
    stop("Study 06 posterior effects are invalid.", call. = FALSE)
  pred <- as.numeric(Z_test %*% effects)
  phenotype <- as.numeric(simulation$phenotype[split$test_ids, 1L])
  genetic <- as.numeric(simulation$genetic_values[split$test_ids, 1L])
  slope <- if (stats::var(pred) > 0)
    unname(stats::coef(stats::lm(phenotype ~ pred))[2L]) else NA_real_
  data.frame(architecture = run$architecture,
    replicate = run$replicate,
    configuration = run$method$configuration,
    method = run$method$native_method,
    metric = c("phenotype_prediction_correlation",
      "phenotype_prediction_rmse", "prediction_regression_slope",
      "genetic_value_correlation", "genetic_value_rmse"),
    value = c(.study06_safe_cor(pred, phenotype),
      sqrt(mean((pred - phenotype)^2)), slope,
      .study06_safe_cor(pred, genetic),
      sqrt(mean((pred - genetic)^2))),
    status = "ok", reason = "", stringsAsFactors = FALSE)
}

.study06_marker_metrics <- function(run, simulation, Z_all) {
  estimated <- as.numeric(run$fit$bm[, 1L])
  truth <- as.numeric(simulation$effects[, 1L])
  causal <- truth != 0
  reconstructed <- as.numeric(Z_all %*% estimated)
  pip <- if (is.null(run$fit$dm)) NULL else
    as.numeric(run$fit$dm[, 1L])
  values <- c(
    marker_effect_correlation = .study06_safe_cor(estimated, truth),
    marker_effect_rmse = sqrt(mean((estimated - truth)^2)),
    causal_marker_effect_correlation =
      .study06_safe_cor(estimated[causal], truth[causal]),
    noncausal_marker_rmse = sqrt(mean(estimated[!causal]^2)),
    genetic_value_reconstruction_correlation =
      .study06_safe_cor(reconstructed,
        as.numeric(simulation$genetic_values)),
    nonnull_probability_mean = if (is.null(pip)) NA_real_ else mean(pip),
    high_pip_marker_count = if (is.null(pip)) NA_real_ else sum(pip >= .5))
  data.frame(architecture = run$architecture,
    replicate = run$replicate,
    configuration = run$method$configuration,
    method = run$method$native_method,
    metric = names(values), value = unname(values),
    available = is.finite(values),
    unavailable_reason = ifelse(is.finite(values), "",
      "output unavailable or zero-variance edge case"),
    stringsAsFactors = FALSE)
}

.study06_comparison_registry <- function() {
  base <- data.frame(
    comparison_id = c("unfiltered_block_eigen_minus_block_csr",
      "block_csr_minus_full_csr",
      "hard_minus_unfiltered_block_eigen",
      "ridge_fixed_minus_unfiltered_block_eigen"),
    focal = c("block_eigen_unfiltered", "block_csr",
      "block_eigen_hard", "block_eigen_ridge_fixed"),
    reference = c("block_csr", "full_csr",
      "block_eigen_unfiltered", "block_eigen_unfiltered"),
    interpretation_level = c(1L, 2L, 3L, 3L),
    stringsAsFactors = FALSE)
  bed <- setdiff(c("full_csr", "block_csr",
    "block_eigen_unfiltered", "block_eigen_hard",
    "block_eigen_ridge_fixed"), "bed")
  rbind(base, data.frame(
    comparison_id = paste0(bed, "_minus_bed"),
    focal = bed, reference = "bed",
    interpretation_level = 4L, stringsAsFactors = FALSE))
}

.study06_operator_pilot_summary <- function(runs, Z_test, split, Z_all,
                                            simulations) {
  find_sim <- function(run) Filter(function(x)
    x$architecture == run$architecture && x$replicate == run$replicate,
    simulations)[[1L]]
  rows <- lapply(runs, function(run) {
    sim <- find_sim(run)
    pred <- .study06_prediction_metrics(run, sim, Z_test, split)
    marker <- .study06_marker_metrics(run, sim, Z_all)
    evidence <- .study06_sbayesr_evidence(run)
    value <- function(table, name) {
      x <- table$value[table$metric == name]
      if (length(x) == 1L) x else NA_real_
    }
    data.frame(architecture = run$architecture, replicate = run$replicate,
      configuration = run$method$configuration,
      posterior_heritability = evidence$benchmark_heritability,
      vgs = evidence$vgs_posterior_mean, ves = evidence$ves_posterior_mean,
      vbs = evidence$vbs_posterior_mean, vld = evidence$vld_posterior_mean,
      vle = evidence$vle_posterior_mean,
      pis = evidence$pis_posterior_mean,
      prediction_correlation = value(pred,
        "phenotype_prediction_correlation"),
      marker_effect_correlation = value(marker, "marker_effect_correlation"),
      runtime_seconds = run$runtime, warnings = run$warnings,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.study06_scalar_equivalence <- function(mean_focal, mean_reference,
                                        mcse_focal, mcse_reference,
                                        lower_focal, upper_focal,
                                        lower_reference, upper_reference,
                                        numerical_tolerance = 1e-6,
                                        standardized_small = 2,
                                        standardized_material = 4) {
  values <- c(mean_focal, mean_reference, mcse_focal, mcse_reference,
    lower_focal, upper_focal, lower_reference, upper_reference)
  if (any(!is.finite(values))) return(list(classification = "indeterminate",
    difference = NA_real_, standardized_difference = NA_real_,
    relative_difference = NA_real_, interval_overlap = NA_real_))
  difference <- mean_focal - mean_reference
  scale <- sqrt(mcse_focal^2 + mcse_reference^2)
  z <- if (scale > 0) abs(difference) / scale else Inf
  relative <- abs(difference) / max(abs(mean_reference), numerical_tolerance)
  overlap <- max(0, min(upper_focal, upper_reference) -
    max(lower_focal, lower_reference)) /
    max(min(upper_focal - lower_focal, upper_reference - lower_reference),
      numerical_tolerance)
  classification <- if (abs(difference) <= numerical_tolerance)
    "numerically_equivalent" else if (z <= standardized_small)
      "consistent_with_monte_carlo_error" else if (z <= standardized_material &&
        relative <= 0.05) "small_but_detectable_difference" else
          "material_difference"
  list(classification = classification, difference = difference,
    standardized_difference = z, relative_difference = relative,
    interval_overlap = overlap)
}

.study06_paired_differences <- function(table, metric_col = "metric",
                                        value_col = "value") {
  registry <- .study06_comparison_registry()
  out <- list()
  for (i in seq_len(nrow(registry))) {
    for (architecture in unique(table$architecture)) {
      focal <- table[table$architecture == architecture &
        table$configuration == registry$focal[i], , drop = FALSE]
      ref <- table[table$architecture == architecture &
        table$configuration == registry$reference[i], , drop = FALSE]
      metrics <- intersect(unique(focal[[metric_col]]),
        unique(ref[[metric_col]]))
      for (metric in metrics) {
        a <- focal[focal[[metric_col]] == metric,
          c("replicate", value_col), drop = FALSE]
        b <- ref[ref[[metric_col]] == metric,
          c("replicate", value_col), drop = FALSE]
        names(a)[2L] <- "focal_value"
        names(b)[2L] <- "reference_value"
        z <- merge(a, b, by = "replicate", all = TRUE)
        z$complete_pair <- is.finite(z$focal_value) &
          is.finite(z$reference_value)
        z$difference <- z$focal_value - z$reference_value
        z$architecture <- architecture
        z$comparison_id <- registry$comparison_id[i]
        z$metric <- metric
        z$orientation <- paste(registry$focal[i], "minus",
          registry$reference[i])
        z$interpretation_level <- registry$interpretation_level[i]
        out[[length(out) + 1L]] <- z
      }
    }
  }
  z <- do.call(rbind, out)
  z[order(z$interpretation_level, z$comparison_id,
    z$architecture, z$metric, z$replicate), ]
}

.study06_paired_summary <- function(x) {
  groups <- split(x, interaction(x$comparison_id,
    x$architecture, x$metric, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    y <- z$difference[z$complete_pair]
    data.frame(comparison_id = z$comparison_id[1L],
      architecture = z$architecture[1L], metric = z$metric[1L],
      orientation = z$orientation[1L],
      complete_paired_replicates = length(y),
      mean_difference = if (length(y)) mean(y) else NA_real_,
      sd_difference = if (length(y) >= 2L) stats::sd(y) else NA_real_,
      median_difference = if (length(y)) stats::median(y) else NA_real_,
      minimum_difference = if (length(y)) min(y) else NA_real_,
      maximum_difference = if (length(y)) max(y) else NA_real_,
      interpretation_level = z$interpretation_level[1L],
      stringsAsFactors = FALSE)
  }))
}

.study06_recovery_summary <- function(estimates) {
  z <- estimates[estimates$available &
    is.finite(estimates$truth), , drop = FALSE]
  groups <- split(z, interaction(z$architecture,
    z$configuration, z$estimand, drop = TRUE))
  do.call(rbind, lapply(groups, function(x)
    data.frame(architecture = x$architecture[1L],
      configuration = x$configuration[1L],
      estimand = x$estimand[1L],
      replicate_count = 5L,
      successful_replicates = nrow(x),
      mean_bias = mean(x$bias),
      rmse = sqrt(mean(x$squared_error)),
      mae = mean(abs(x$bias)),
      observed_coverage_count =
        sum(x$interval_coverage, na.rm = TRUE),
      observed_coverage_proportion =
        mean(x$interval_coverage, na.rm = TRUE),
      mean_interval_width = mean(x$upper_95 - x$lower_95),
      coverage_interpretation =
        "coarse five-replicate observed coverage; each miss changes 20 percentage points",
      stringsAsFactors = FALSE)))
}

.study06_sbayesr_evidence <- function(run) {
  fit <- run$fit
  scalar <- function(x) {
    if (is.null(x)) return(NA_real_)
    x <- as.numeric(unlist(x, use.names = FALSE))
    if (!length(x) || !any(is.finite(x))) NA_real_ else
      mean(x[is.finite(x)])
  }
  vg <- scalar(fit$vgs)
  ve <- scalar(fit$ves)
  data.frame(architecture = run$architecture,
    replicate = run$replicate,
    configuration = run$method$configuration,
    method = run$method$native_method,
    vgs_posterior_mean = vg, ves_posterior_mean = ve,
    vbs_posterior_mean = scalar(fit$vbs),
    pis_posterior_mean = scalar(fit$pi_trace %||%
      lapply(fit$chains, `[[`, "pis")),
    vld_posterior_mean = scalar(fit$vld),
    vle_posterior_mean = scalar(fit$vle),
    benchmark_heritability = vg / (vg + ve),
    package_reported_heritability = scalar(fit$heritability),
    marker_effect_variance = stats::var(as.numeric(fit$bm)),
    operator_diagonal_minimum = if (!is.null(
      fit$data$operator$diagonal))
        min(fit$data$operator$diagonal) else NA_real_,
    operator_diagonal_mean = if (!is.null(
      fit$data$operator$diagonal))
        mean(fit$data$operator$diagonal) else NA_real_,
    log_cpo = scalar(fit$log_cpo),
    mean_log_cpo = scalar(fit$mean_log_cpo),
    stringsAsFactors = FALSE)
}

.study06_marker_agreement <- function(runs) {
  registry <- .study06_comparison_registry()
  out <- list()
  architectures <- unique(vapply(runs, `[[`, "", "architecture"))
  for (architecture in architectures)
    for (replicate in sort(unique(vapply(runs, `[[`, integer(1),
      "replicate"))))
      for (i in seq_len(nrow(registry))) {
        focal <- Filter(function(x) x$architecture == architecture &&
          x$replicate == replicate &&
          x$method$configuration == registry$focal[i], runs)
        reference <- Filter(function(x) x$architecture == architecture &&
          x$replicate == replicate &&
          x$method$configuration == registry$reference[i], runs)
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
          posterior_mean_effect_correlation = .study06_safe_cor(a, b),
          posterior_mean_effect_rmse = sqrt(mean((a - b)^2)),
          nonnull_probability_correlation =
            if (is.null(da) || is.null(db)) NA_real_ else
              .study06_safe_cor(da, db),
          high_pip_overlap = if (is.null(da) || is.null(db)) NA_real_ else {
            A <- which(da >= .5); B <- which(db >= .5)
            if (!length(union(A, B))) 1 else
              length(intersect(A, B)) / length(union(A, B))
          },
          component_probability_agreement =
            if (is.null(focal[[1L]]$fit$component_probabilities) ||
                is.null(reference[[1L]]$fit$component_probabilities))
              NA_real_ else .study06_safe_cor(
                as.numeric(focal[[1L]]$fit$component_probabilities[[1L]]),
                as.numeric(reference[[1L]]$fit$component_probabilities[[1L]])),
          stringsAsFactors = FALSE)
      }
  do.call(rbind, out)
}

.study06_convergence_validation_summary <- function(diagnostics) {
  fit_pass <- do.call(rbind, lapply(split(diagnostics,
    interaction(diagnostics$architecture, diagnostics$replicate,
      diagnostics$configuration, drop = TRUE)), function(x)
    data.frame(architecture = x$architecture[1L],
      replicate = x$replicate[1L],
      configuration = x$configuration[1L],
      pass = all(x$overall_pass),
      limiting_estimand = x$estimand[which.max(pmax(
        x$rhat / 1.01, 400 / x$ess_bulk, 400 / x$ess_tail,
        x$relative_mcse / .05))],
      stringsAsFactors = FALSE)))
  groups <- split(fit_pass, interaction(fit_pass$architecture,
    fit_pass$configuration, drop = TRUE))
  do.call(rbind, lapply(groups, function(x) {
    d <- diagnostics[diagnostics$architecture == x$architecture[1L] &
      diagnostics$configuration == x$configuration[1L], ,
      drop = FALSE]
    limiting <- names(sort(table(x$limiting_estimand),
      decreasing = TRUE))[1L]
    count <- sum(x$pass)
    data.frame(architecture = x$architecture[1L],
      configuration = x$configuration[1L],
      replicate_count = 5L, pass_count = count,
      pass_proportion = count / 5,
      maximum_observed_rhat = max(d$rhat),
      minimum_observed_bulk_ess = min(d$ess_bulk),
      minimum_observed_tail_ess = min(d$ess_tail),
      maximum_observed_relative_mcse = max(d$relative_mcse),
      limiting_estimand = limiting,
      limiting_diagnostic = d$limiting_diagnostic[
        which.max(pmax(d$rhat / 1.01, 400 / d$ess_bulk,
          400 / d$ess_tail, d$relative_mcse / .05))],
      limiting_replicate = d$replicate[
        which.max(pmax(d$rhat / 1.01, 400 / d$ess_bulk,
          400 / d$ess_tail, d$relative_mcse / .05))],
      recommendation_validation_status = if (count == 5L)
        "supported_in_all_replicates" else if (count > 0L)
          "partially_supported" else "not_supported",
      interpretation =
        "five-simulation support is not universal convergence validation",
      stringsAsFactors = FALSE)
  }))
}
