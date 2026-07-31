.study05_cor <- function(x, y) {
  if (length(x) != length(y) || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(NA_real_)
  stats::cor(x, y)
}

.study05_rmse <- function(x, y) sqrt(mean((x - y)^2))

.study05_auroc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || all(truth)) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[truth]) - sum(seq_len(sum(truth)))) /
    (sum(truth) * sum(!truth))
}

.study05_auprc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth)) return(NA_real_)
  o <- order(score, decreasing = TRUE)
  y <- truth[o]
  mean(cumsum(y)[y] / which(y))
}

.study05_component_probabilities <- function(fit) {
  x <- fit$component_probabilities
  if (is.list(x)) x <- x[[1L]]
  x <- as.matrix(x)
  if (!nrow(x) || any(!is.finite(x)) || any(x < 0 | x > 1) ||
      any(abs(rowSums(x) - 1) > 1e-8)) return(NULL)
  x
}

.study05_prior_probabilities <- function(run, A = NULL) {
  if (!isTRUE(run$method$annotation)) return(NULL)
  x <- run$fit$annotation_prior
  if (is.list(x)) x <- x[[1L]]
  if (is.null(x) && !is.null(A) && length(run$fit$alpha)) {
    alpha <- if (is.list(run$fit$alpha)) run$fit$alpha[[1L]] else run$fit$alpha
    x <- sblr::sbayesrc_marker_pi(A, alpha,
      gamma = run$controls$mixture_var)
  }
  if (is.null(x)) return(NULL)
  x <- as.matrix(x)
  if (!nrow(x) || any(!is.finite(x)) || any(x < 0 | x > 1) ||
      any(abs(rowSums(x) - 1) > 1e-8)) return(NULL)
  x
}

.study05_probability_inventory <- function(run, A, component_count) {
  posterior <- .study05_component_probabilities(run$fit)
  prior <- .study05_prior_probabilities(run, A)
  marker_count <- nrow(A)
  one <- function(x, object, expected) data.frame(
    scenario = run$scenario, replicate = run$replicate, method = run$method$id,
    probability_object = object, available = !is.null(x),
    row_count = if (is.null(x)) NA_integer_ else nrow(x),
    column_count = if (is.null(x)) NA_integer_ else ncol(x),
    expected_row_count = marker_count, expected_column_count = component_count,
    minimum_probability = if (is.null(x)) NA_real_ else min(x),
    maximum_probability = if (is.null(x)) NA_real_ else max(x),
    maximum_row_sum_error = if (is.null(x)) NA_real_ else
      max(abs(rowSums(x) - 1)),
    required = expected,
    reason = if (!is.null(x)) "" else if (expected)
      "required output missing" else "not returned by this model",
    stringsAsFactors = FALSE)
  rbind(one(posterior, "posterior_marker_component_allocation", TRUE),
    one(prior, "annotation_implied_marker_prior",
      isTRUE(run$method$annotation)))
}

.study05_effects <- function(fit, marker_ids) {
  x <- fit$bm
  if (is.matrix(x)) x <- x[, 1L]
  x <- as.numeric(x)
  ids <- rownames(fit$bm)
  if (is.null(ids)) ids <- marker_ids
  names(x) <- ids
  x[match(marker_ids, names(x))]
}

.study05_prediction_metrics <- function(run, simulation, Z_test, test_ids) {
  if (run$status != "ok") return(data.frame())
  effect <- .study05_effects(run$fit, names(simulation$effect))
  prediction <- as.numeric(Z_test %*% effect)
  phenotype <- simulation$phenotype[test_ids, 1L]
  genetic <- simulation$genetic_value[test_ids]
  slope <- if (stats::sd(prediction) > 0)
    unname(stats::coef(stats::lm(phenotype ~ prediction))[2L]) else NA_real_
  values <- c(
    phenotype_prediction_correlation = .study05_cor(prediction, phenotype),
    phenotype_prediction_rmse = .study05_rmse(prediction, phenotype),
    prediction_regression_slope = slope,
    genetic_value_correlation = .study05_cor(prediction, genetic),
    genetic_value_rmse = .study05_rmse(prediction, genetic))
  data.frame(scenario = simulation$scenario, replicate = simulation$replicate,
    method = run$method$id, metric = names(values), value = as.numeric(values),
    available = is.finite(values),
    reason = ifelse(is.finite(values), "", "zero variance or non-finite prediction"),
    stringsAsFactors = FALSE)
}

.study05_marker_metrics <- function(run, simulation, A, config) {
  if (run$status != "ok") return(list(marker = data.frame(),
    probability = data.frame(), enrichment = data.frame()))
  truth <- simulation$marker_truth
  effect <- .study05_effects(run$fit, truth$marker_id)
  active <- truth$true_nonnull
  posterior <- .study05_component_probabilities(run$fit)
  if (!is.null(posterior) && !identical(rownames(posterior), truth$marker_id))
    posterior <- posterior[match(truth$marker_id, rownames(posterior)), , drop = FALSE]
  marker_values <- c(
    marker_effect_correlation = .study05_cor(effect, truth$effect),
    marker_effect_rmse = .study05_rmse(effect, truth$effect),
    causal_marker_effect_correlation = .study05_cor(effect[active], truth$effect[active]),
    noncausal_marker_rmse = .study05_rmse(effect[!active], truth$effect[!active]))
  marker <- data.frame(scenario = simulation$scenario,
    replicate = simulation$replicate, method = run$method$id,
    metric = names(marker_values), value = as.numeric(marker_values),
    available = is.finite(marker_values),
    reason = ifelse(is.finite(marker_values), "", "zero variance edge case"),
    stringsAsFactors = FALSE)
  probability <- data.frame()
  if (!is.null(posterior)) {
    nonnull <- 1 - posterior[, 1L]
    expected_variance <- as.numeric(posterior %*% config$mixture_var)
    true_onehot <- matrix(0, nrow(truth), length(config$mixture_var))
    true_onehot[cbind(seq_len(nrow(truth)), truth$component_index + 1L)] <- 1
    pclip <- pmax(posterior, 1e-12)
    posterior_values <- c(
      posterior_nonnull_auroc = .study05_auroc(nonnull, active),
      posterior_nonnull_auprc = .study05_auprc(nonnull, active),
      posterior_nonnull_brier = mean((nonnull - active)^2),
      posterior_component_accuracy =
        mean(max.col(posterior, ties.method = "first") - 1L == truth$component_index),
      posterior_multiclass_brier = mean(rowSums((posterior - true_onehot)^2)),
      posterior_component_log_score =
        -mean(log(pclip[cbind(seq_len(nrow(truth)), truth$component_index + 1L)])),
      posterior_expected_variance_correlation =
        .study05_cor(expected_variance, truth$component_variance))
    probability <- rbind(probability, data.frame(
      scenario = simulation$scenario, replicate = simulation$replicate,
      method = run$method$id, probability_object = "posterior_marker_component_allocation",
      metric = names(posterior_values), value = as.numeric(posterior_values),
      available = is.finite(posterior_values), reason = "", stringsAsFactors = FALSE))
  }
  prior <- .study05_prior_probabilities(run, A)
  if (!is.null(prior)) {
    prior <- as.matrix(prior)
    if (!identical(rownames(prior), truth$marker_id))
      prior <- prior[match(truth$marker_id, rownames(prior)), , drop = FALSE]
    true_prior <- simulation$prior_probability
    prior_values <- c(
      prior_component_probability_correlation =
        .study05_cor(as.numeric(prior), as.numeric(true_prior)),
      prior_component_probability_rmse = .study05_rmse(prior, true_prior),
      prior_component_probability_mae = mean(abs(prior - true_prior)),
      prior_multiclass_brier = mean(rowSums((prior - true_prior)^2)),
      prior_component_log_score = -mean(rowSums(true_prior * log(pmax(prior, 1e-12)))),
      prior_expected_nonnull_probability_error =
        mean(1 - prior[, 1L]) - mean(1 - true_prior[, 1L]))
    probability <- rbind(probability, data.frame(
      scenario = simulation$scenario, replicate = simulation$replicate,
      method = run$method$id, probability_object = "annotation_implied_marker_prior",
      metric = names(prior_values), value = as.numeric(prior_values),
      available = is.finite(prior_values), reason = "", stringsAsFactors = FALSE))
  } else {
    metrics <- c("prior_component_probability_correlation",
      "prior_component_probability_rmse", "prior_component_probability_mae",
      "prior_multiclass_brier", "prior_component_log_score",
      "prior_expected_nonnull_probability_error")
    probability <- rbind(probability, data.frame(
      scenario = simulation$scenario, replicate = simulation$replicate,
      method = run$method$id, probability_object = "annotation_implied_marker_prior",
      metric = metrics, value = NA_real_, available = FALSE,
      reason = "method has no annotation-implied marker prior output",
      stringsAsFactors = FALSE))
  }
  enriched <- A[, "enriched_binary"] == 1L
  score <- if (!is.null(prior)) 1 - prior[, 1L] else
    if (!is.null(posterior)) 1 - posterior[, 1L] else rep(NA_real_, nrow(A))
  source <- if (!is.null(prior)) "annotation_implied_marker_prior" else
    "posterior_marker_nonnull_probability"
  enrichment_values <- c(
    true_nonnull_probability_fraction_enriched =
      simulation$truth$true_enriched_nonnull_probability_fraction,
    realized_causal_fraction_enriched =
      simulation$truth$realized_enriched_causal_fraction,
    mean_nonnull_probability_enriched = mean(score[enriched]),
    mean_nonnull_probability_unannotated = mean(score[!enriched]),
    enrichment_ratio = mean(score[enriched]) / mean(score[!enriched]),
    enrichment_difference = mean(score[enriched]) - mean(score[!enriched]),
    enriched_ranking_auroc = .study05_auroc(score, enriched),
    continuous_annotation_direction = .study05_cor(score, A[, "continuous_signal"]))
  enrichment <- data.frame(scenario = simulation$scenario,
    replicate = simulation$replicate, method = run$method$id,
    probability_object = source, metric = names(enrichment_values),
    value = as.numeric(enrichment_values),
    available = is.finite(enrichment_values),
    reason = ifelse(is.finite(enrichment_values), "", "probability output unavailable"),
    stringsAsFactors = FALSE)
  list(marker = marker, probability = probability, enrichment = enrichment)
}

.study05_scalar_estimates <- function(run, draws, simulation) {
  ids <- c(effect_variance = "vbs", genetic_variance = "vgs",
    residual_variance = "ves", heritability = "heritability",
    global_nonnull_proportion = "global_nonnull_proportion")
  rows <- lapply(names(ids), function(estimand) {
    z <- draws$value[draws$quantity == ids[[estimand]] |
      draws$quantity == estimand]
    truth <- simulation$truth[[estimand]]
    data.frame(scenario = simulation$scenario, replicate = simulation$replicate,
      method = run$method$id, estimand = estimand, available = length(z) > 0L,
      reason = if (length(z)) "" else "chain draws unavailable",
      truth = truth, posterior_mean = if (length(z)) mean(z) else NA_real_,
      posterior_sd = if (length(z)) stats::sd(z) else NA_real_,
      posterior_median = if (length(z)) stats::median(z) else NA_real_,
      lower_95 = if (length(z)) unname(stats::quantile(z, .025)) else NA_real_,
      upper_95 = if (length(z)) unname(stats::quantile(z, .975)) else NA_real_,
      stringsAsFactors = FALSE)
  })
  unavailable <- lapply(paste0("global_component_proportion_", 0:3),
    function(estimand) data.frame(scenario = simulation$scenario,
      replicate = simulation$replicate, method = run$method$id,
      estimand = estimand, available = FALSE,
      reason = "no comparable per-iteration scalar chain draws",
      truth = if (estimand == "global_nonnull_proportion")
        simulation$truth$global_nonnull_proportion else
          simulation$truth$global_component_proportions[
            as.integer(sub(".*_", "", estimand)) + 1L],
      posterior_mean = NA_real_, posterior_sd = NA_real_,
      posterior_median = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_,
      stringsAsFactors = FALSE))
  out <- do.call(rbind, c(rows, unavailable))
  out$bias <- out$posterior_mean - out$truth
  out$squared_error <- out$bias^2
  out$interval_coverage <- with(out, available & lower_95 <= truth & upper_95 >= truth)
  out
}

.study05_annotation_estimates <- function(run, draws, simulation) {
  if (!isTRUE(run$method$annotation)) return(data.frame())
  z <- draws[draws$parameter_name == "alpha", , drop = FALSE]
  groups <- split(z, z$quantity)
  do.call(rbind, lapply(groups, function(x) {
    truth <- simulation$alpha[x$annotation_name[1L], x$stick_name[1L]]
    p_correct <- if (truth > 0) mean(x$value > 0) else if (truth < 0)
      mean(x$value < 0) else NA_real_
    lo <- unname(stats::quantile(x$value, .025))
    hi <- unname(stats::quantile(x$value, .975))
    data.frame(scenario = simulation$scenario, replicate = simulation$replicate,
      method = run$method$id, annotation_name = x$annotation_name[1L],
      stick_name = x$stick_name[1L], true_alpha = truth,
      posterior_mean = mean(x$value), posterior_sd = stats::sd(x$value),
      posterior_median = stats::median(x$value), lower_95 = lo, upper_95 = hi,
      bias = mean(x$value) - truth, squared_error = (mean(x$value) - truth)^2,
      interval_coverage = lo <= truth && hi >= truth,
      sign_recovery = if (truth == 0) NA else sign(mean(x$value)) == sign(truth),
      posterior_probability_correct_sign = p_correct,
      zero_covered = lo <= 0 && hi >= 0, stringsAsFactors = FALSE)
  }))
}

.study05_sigma_alpha <- function(run, draws) {
  if (!isTRUE(run$method$annotation)) return(data.frame())
  z <- draws[draws$parameter_name == "sigmaSqAlpha", , drop = FALSE]
  do.call(rbind, lapply(split(z, z$quantity), function(x)
    data.frame(scenario = run$scenario, replicate = run$replicate,
      method = run$method$id, stick_name = x$stick_name[1L],
      interpretation = "posterior_behaviour_no_generating_truth",
      configured_prior_a = run$controls$sigmaSqAlpha_a,
      configured_prior_b = run$controls$sigmaSqAlpha_b,
      posterior_mean = mean(x$value), posterior_sd = stats::sd(x$value),
      posterior_median = stats::median(x$value),
      lower_95 = unname(stats::quantile(x$value, .025)),
      upper_95 = unname(stats::quantile(x$value, .975)),
      stringsAsFactors = FALSE)))
}

.study05_aggregate <- function(x, value = "value",
                               keys = c("scenario", "method", "metric")) {
  groups <- split(x, interaction(x[keys], drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    y <- z[[value]][z$available %in% TRUE & is.finite(z[[value]])]
    cbind(z[1L, keys, drop = FALSE],
      data.frame(replicate_count = 5L, successful_replicates = length(y),
        mean = if (length(y)) mean(y) else NA_real_,
        sd = if (length(y) >= 2L) stats::sd(y) else NA_real_,
        median = if (length(y)) stats::median(y) else NA_real_,
        minimum = if (length(y)) min(y) else NA_real_,
        maximum = if (length(y)) max(y) else NA_real_,
        stringsAsFactors = FALSE))
  }))
}

.study05_comparison_registry <- function() data.frame(
  comparison_id = c("bed_rc_minus_r", "bed_rc_minus_r",
    "csr_rc_minus_r", "csr_rc_minus_r", "csr_rc_minus_bed_rc",
    "csr_r_minus_bed_r"),
  scenario = c("informative_annotations", "uninformative_annotations",
    "informative_annotations", "uninformative_annotations", "*", "*"),
  focal_method = c("st_bed_bayesrc", "st_bed_bayesrc",
    "st_csr_sbayesrc", "st_csr_sbayesrc", "st_csr_sbayesrc",
    "st_csr_sbayesr"),
  comparison_method = c("st_bed_bayesr", "st_bed_bayesr",
    "st_csr_sbayesr", "st_csr_sbayesr", "st_bed_bayesrc",
    "st_bed_bayesr"), stringsAsFactors = FALSE)

.study05_paired <- function(metrics) {
  registry <- .study05_comparison_registry()
  rows <- list()
  for (i in seq_len(nrow(registry))) {
    scenarios <- if (registry$scenario[i] == "*") unique(metrics$scenario) else registry$scenario[i]
    for (scenario in scenarios) {
      f <- metrics[metrics$scenario == scenario &
        metrics$method == registry$focal_method[i], ]
      c <- metrics[metrics$scenario == scenario &
        metrics$method == registry$comparison_method[i], ]
      z <- merge(f, c, by = c("scenario", "replicate", "metric"),
        suffixes = c("_focal", "_comparison"))
      if (nrow(z)) rows[[length(rows) + 1L]] <- data.frame(
        comparison_id = registry$comparison_id[i], scenario = scenario,
        replicate = z$replicate, metric = z$metric,
        focal_method = registry$focal_method[i],
        comparison_method = registry$comparison_method[i],
        difference = z$value_focal - z$value_comparison,
        complete_pair = is.finite(z$value_focal) & is.finite(z$value_comparison),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

.study05_paired_summary <- function(x) {
  groups <- split(x, interaction(x$comparison_id, x$scenario, x$metric, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    y <- z$difference[z$complete_pair]
    data.frame(comparison_id = z$comparison_id[1L], scenario = z$scenario[1L],
      metric = z$metric[1L], complete_paired_replicates = length(y),
      mean_difference = if (length(y)) mean(y) else NA_real_,
      sd_difference = if (length(y) >= 2L) stats::sd(y) else NA_real_,
      median_difference = if (length(y)) stats::median(y) else NA_real_,
      minimum_difference = if (length(y)) min(y) else NA_real_,
      maximum_difference = if (length(y)) max(y) else NA_real_,
      stringsAsFactors = FALSE)
  }))
}

.study05_interactions <- function(paired) {
  x <- paired[paired$comparison_id %in% c("bed_rc_minus_r", "csr_rc_minus_r"), ]
  info <- x[x$scenario == "informative_annotations", ]
  null <- x[x$scenario == "uninformative_annotations", ]
  z <- merge(info, null, by = c("comparison_id", "replicate", "metric"),
    suffixes = c("_informative", "_uninformative"))
  data.frame(comparison_id = z$comparison_id, replicate = z$replicate,
    metric = z$metric,
    interaction_difference = z$difference_informative - z$difference_uninformative,
    complete_pair = z$complete_pair_informative & z$complete_pair_uninformative,
    stringsAsFactors = FALSE)
}
