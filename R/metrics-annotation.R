# Stable Study 06 annotation-recovery metrics. These functions consume tidy
# extracted quantities and never inspect native fit objects.

.annotation_metric_row <- function(metadata, metric, value, status = "ok",
                                   reason = "") {
  data.frame(study = metadata$study, scenario = metadata$scenario,
    replicate = as.integer(metadata$replicate), method = metadata$method,
    metric = metric, value = as.numeric(value), status = status,
    reason = reason, stringsAsFactors = FALSE)
}

#' Summarise annotation coefficient recovery
#'
#' @param traces Tidy true alpha traces.
#' @param true_alpha True annotation-by-stick coefficient matrix.
#' @param metadata Named study/scenario/replicate/method list.
#' @return A tidy coefficient recovery data frame.
#' @export
annotation_alpha_recovery <- function(traces, true_alpha, metadata) {
  if (!is.data.frame(traces) || !nrow(traces) ||
      any(traces$status != "ok"))
    return(data.frame(study = metadata$study, scenario = metadata$scenario,
      replicate = as.integer(metadata$replicate), method = metadata$method,
      annotation = NA_character_, stick = NA_character_, truth = NA_real_,
      posterior_mean = NA_real_, posterior_sd = NA_real_, lower_95 = NA_real_,
      upper_95 = NA_real_, bias = NA_real_, squared_error = NA_real_,
      covered_95 = FALSE, interval_width_95 = NA_real_,
      status = "unavailable", reason = "True alpha traces are unavailable.",
      stringsAsFactors = FALSE))
  alpha <- traces[traces$parameter == "alpha", , drop = FALSE]
  groups <- split(alpha, interaction(alpha$annotation, alpha$stick,
    drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(x) {
    truth <- true_alpha[x$annotation[1L], x$stick[1L]]
    interval <- unname(stats::quantile(x$value, c(.025, .975)))
    estimate <- mean(x$value)
    data.frame(study = metadata$study, scenario = metadata$scenario,
      replicate = as.integer(metadata$replicate), method = metadata$method,
      annotation = x$annotation[1L], stick = x$stick[1L], truth = truth,
      posterior_mean = estimate, posterior_sd = stats::sd(x$value),
      lower_95 = interval[1L], upper_95 = interval[2L],
      bias = estimate - truth, squared_error = (estimate - truth)^2,
      covered_95 = interval[1L] <= truth && interval[2L] >= truth,
      interval_width_95 = diff(interval), status = "ok", reason = "",
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Calculate annotation-prior recovery metrics
#'
#' @param prior Result from `summarise_drawwise_annotation_prior()`.
#' @param true_marker_prior Ordered true marker component probabilities.
#' @param annotations Ordered annotation matrix.
#' @param marker_truth Ordered marker truth table.
#' @param metadata Named study/scenario/replicate/method list.
#' @return A tidy metric data frame.
#' @export
annotation_prior_recovery <- function(prior, true_marker_prior, annotations,
                                      marker_truth, metadata) {
  if (!identical(prior$status, "ok"))
    return(.annotation_metric_row(metadata,
      "annotation_prior_draw_summary", NA_real_, "unavailable", prior$reason))
  marker <- prior$marker
  truth <- true_marker_prior[match(marker$marker_id,
    rownames(true_marker_prior)), , drop = FALSE]
  estimated_columns <- grep("^posterior_mean_prior_component_", names(marker),
    value = TRUE)
  estimated <- as.matrix(marker[estimated_columns])
  enriched <- annotations[, "enriched_binary"] == 1
  true_nonnull <- 1 - truth[, 1L]
  estimated_nonnull <- marker$posterior_mean_nonnull_prior
  causal <- marker_truth$true_nonnull[match(marker$marker_id,
    marker_truth$marker_id)]
  correlation <- function(x, y)
    if (stats::sd(x) == 0 || stats::sd(y) == 0) NA_real_ else
      stats::cor(x, y)
  values <- c(
    prior_component_probability_rmse = sqrt(mean((estimated - truth)^2)),
    prior_component_probability_correlation =
      correlation(as.numeric(estimated), as.numeric(truth)),
    prior_expected_active_error = sum(estimated_nonnull) - sum(true_nonnull),
    posterior_mean_prior_enriched = mean(estimated_nonnull[enriched]),
    posterior_mean_prior_unannotated = mean(estimated_nonnull[!enriched]),
    enriched_prior_contrast =
      mean(estimated_nonnull[enriched]) - mean(estimated_nonnull[!enriched]),
    causal_vs_noncausal_prior_contrast =
      mean(estimated_nonnull[causal]) - mean(estimated_nonnull[!causal]),
    true_signal_prior_correlation = correlation(true_nonnull,
      estimated_nonnull),
    continuous_signal_prior_correlation = correlation(
      annotations[, "continuous_signal"], estimated_nonnull),
    null_annotation_prior_correlation = correlation(
      annotations[, "null_annotation"], estimated_nonnull))
  do.call(rbind, Map(function(id, value) .annotation_metric_row(metadata, id,
    value, if (is.finite(value)) "ok" else "unavailable",
    if (is.finite(value)) "" else "zero-variance correlation"),
    names(values), values))
}

#' Calculate Study 06 marker prioritisation and effect-recovery metrics
#'
#' @param marker_results Data frame with marker ID, posterior inclusion
#'   probability, and posterior mean effect.
#' @param marker_truth Data frame with marker ID, causal status, and true effect.
#' @param metadata Named study/scenario/replicate/method list.
#' @param top_k Positive integer cutoffs.
#' @return A tidy metric data frame.
#' @export
annotation_marker_recovery <- function(marker_results, marker_truth, metadata,
                                       top_k = c(10L, 50L, 100L)) {
  required <- c("marker_id", "posterior_inclusion_probability",
    "posterior_mean_effect")
  if (!is.data.frame(marker_results) ||
      !all(required %in% names(marker_results)))
    stop("marker_results lacks required Study 06 columns.", call. = FALSE)
  truth <- marker_truth[match(marker_results$marker_id,
    marker_truth$marker_id), , drop = FALSE]
  if (anyNA(truth$marker_id))
    stop("Marker recovery truth is not aligned.", call. = FALSE)
  score <- marker_results$posterior_inclusion_probability
  causal <- as.logical(truth$true_nonnull)
  effect <- marker_results$posterior_mean_effect
  rank_desc <- rank(-score, ties.method = "min")
  average_precision <- function(score, causal) {
    if (!any(causal)) return(NA_real_)
    ranked <- causal[order(-score, seq_along(score))]
    mean(cumsum(ranked)[ranked] / which(ranked))
  }
  correlation <- function(x, y)
    if (stats::sd(x) == 0 || stats::sd(y) == 0) NA_real_ else
      stats::cor(x, y)
  values <- c(
    causal_marker_pip = mean(score[causal]),
    causal_marker_rank = mean(rank_desc[causal]),
    pip_auprc = average_precision(score, causal),
    pip_auroc = if (!any(causal) || all(causal)) NA_real_ else {
      r <- rank(score, ties.method = "average")
      (sum(r[causal]) - sum(seq_len(sum(causal)))) /
        (sum(causal) * sum(!causal))
    },
    pip_brier = mean((score - causal)^2),
    effect_rmse = sqrt(mean((effect - truth$effect)^2)),
    effect_correlation = correlation(effect, truth$effect))
  for (k in as.integer(top_k)) {
    selected <- order(score, decreasing = TRUE)[seq_len(min(k,
      length(score)))]
    true_positive <- sum(causal[selected])
    values[[paste0("top_", k, "_recovery")]] <- true_positive / sum(causal)
    values[[paste0("top_", k, "_precision")]] <- true_positive /
      length(selected)
    values[[paste0("top_", k, "_recall")]] <- true_positive / sum(causal)
  }
  do.call(rbind, Map(function(id, value) .annotation_metric_row(metadata, id,
    value, if (is.finite(value)) "ok" else "unavailable",
    if (is.finite(value)) "" else "metric is undefined for this fixture"),
    names(values), values))
}

#' Calculate Study 06 variance-parameter recovery metrics
#'
#' @param estimates Tidy parameter estimates with `parameter`,
#'   `posterior_mean`, and `truth` columns.
#' @param metadata Named study/scenario/replicate/method list.
#' @return Tidy bias, absolute-error, and squared-error rows.
#' @export
annotation_parameter_recovery <- function(estimates, metadata) {
  required <- c("parameter", "posterior_mean", "truth", "status", "reason")
  if (!is.data.frame(estimates) || !all(required %in% names(estimates)))
    stop("Annotation parameter estimates lack recovery columns.",
      call. = FALSE)
  rows <- list()
  for (i in seq_len(nrow(estimates))) {
    ok <- identical(estimates$status[i], "ok") &&
      is.finite(estimates$posterior_mean[i]) && is.finite(estimates$truth[i])
    error <- estimates$posterior_mean[i] - estimates$truth[i]
    values <- c(bias = error, absolute_error = abs(error),
      squared_error = error^2)
    for (metric in names(values)) rows[[length(rows) + 1L]] <-
      .annotation_metric_row(metadata,
        paste(estimates$parameter[i], metric, sep = "_"),
        if (ok) values[[metric]] else NA_real_,
        if (ok) "ok" else "unavailable",
        if (ok) "" else estimates$reason[i])
  }
  do.call(rbind, rows)
}

#' Create matched annotation-method comparisons
#'
#' @param metrics Tidy Study 06 metric table.
#' @return Paired annotation-minus-baseline and scenario-interaction rows.
#' @export
annotation_paired_advantages <- function(metrics) {
  pairs <- data.frame(interface = c("BED", "CSR"),
    annotation_method = c("st_bed_bayesrc", "st_csr_sbayesrc"),
    baseline_method = c("st_bed_bayesr", "st_csr_sbayesr"),
    stringsAsFactors = FALSE)
  rows <- list()
  for (i in seq_len(nrow(pairs))) for (scenario in unique(metrics$scenario)) {
    a <- metrics[metrics$scenario == scenario &
      metrics$method == pairs$annotation_method[i], ]
    b <- metrics[metrics$scenario == scenario &
      metrics$method == pairs$baseline_method[i], ]
    z <- merge(a, b, by = c("study", "scenario", "replicate", "metric"),
      suffixes = c("_annotation", "_baseline"))
    if (nrow(z)) rows[[length(rows) + 1L]] <- data.frame(study = z$study,
      interface = pairs$interface[i], scenario = z$scenario,
      replicate = z$replicate, metric = z$metric,
      annotation_method = pairs$annotation_method[i],
      baseline_method = pairs$baseline_method[i],
      difference = z$value_annotation - z$value_baseline,
      status = ifelse(z$status_annotation == "ok" & z$status_baseline == "ok",
        "ok", "unavailable"), stringsAsFactors = FALSE)
  }
  paired <- if (length(rows)) do.call(rbind, rows) else NULL
  if (is.null(paired)) return(NULL)
  informative <- paired[paired$scenario == "informative_annotations", ]
  uninformative <- paired[paired$scenario == "uninformative_annotations", ]
  interaction <- merge(informative, uninformative,
    by = c("study", "interface", "replicate", "metric",
      "annotation_method", "baseline_method"),
    suffixes = c("_informative", "_uninformative"))
  interaction$scenario <- "informative_minus_uninformative"
  interaction$difference <- interaction$difference_informative -
    interaction$difference_uninformative
  interaction$status <- ifelse(interaction$status_informative == "ok" &
    interaction$status_uninformative == "ok", "ok", "unavailable")
  keep <- names(paired)
  rbind(paired, interaction[keep])
}
