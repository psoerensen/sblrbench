# Extraction distinguishes native formatted fits from benchmark result objects
# and never infers unavailable traces or probabilities.

.benchmark_native_fit <- function(fit) {
  if (inherits(fit, "sblrbench_result")) return(fit$native_fit)
  if (is.list(fit) && !is.null(fit$native_fit)) return(fit$native_fit)
  fit
}

extract_marker_effects <- function(fit) {
  value <- if (inherits(fit, "sblrbench_result")) fit$estimates$effects else
    .benchmark_native_fit(fit)$bm
  if (is.null(value)) stop("Posterior mean marker effects are unavailable.",
    call. = FALSE)
  attr(value, "benchmark_quantity") <- "posterior_mean"
  value
}

extract_marker_probabilities <- function(fit) {
  native <- .benchmark_native_fit(fit)
  posterior <- if (inherits(fit, "sblrbench_result")) fit$estimates$pip else
    native$dm
  list(posterior_inclusion = posterior,
    posterior_component = native$comp_prob,
    prior_component = native$pi,
    prior_marker_component = native$pim)
}

extract_variance_components <- function(fit) {
  native <- .benchmark_native_fit(fit)
  list(posterior_means = native[intersect(c("vbs", "vgs", "ves", "vle",
    "vld"), names(native))],
    true_traces = if (is.list(native$chains)) native$chains else NULL,
    final_states = native[intersect(c("vb", "vg", "ve"), names(native))])
}

extract_fit_metadata <- function(fit) {
  native <- .benchmark_native_fit(fit)
  data.frame(
    method = if (inherits(fit, "sblrbench_result")) fit$method_id else NA_character_,
    marker_count = if (!is.null(native$bm)) nrow(native$bm) else NA_integer_,
    trait_count = if (!is.null(native$bm)) ncol(native$bm) else NA_integer_,
    has_posterior_probabilities = !is.null(native$dm),
    has_component_probabilities = !is.null(native$comp_prob),
    has_true_chain_traces = is.list(native$chains),
    stringsAsFactors = FALSE)
}

extract_runtime <- function(fit) {
  value <- if (inherits(fit, "sblrbench_result"))
    fit$computation$elapsed_seconds else fit$elapsed_seconds
  if (is.null(value) || !length(value)) NA_real_ else as.numeric(value[[1L]])
}

extract_chain_information <- function(fit) {
  native <- .benchmark_native_fit(fit)
  diagnostics <- if (inherits(fit, "sblrbench_result"))
    fit$diagnostics$convergence else native$convergence
  list(true_traces = if (is.list(native$chains)) native$chains else NULL,
    compact_summaries = diagnostics,
    final_states = native[intersect(c("b", "d", "vb", "vg", "ve", "pi"),
      names(native))])
}

prediction_estimate_table <- function(result, scenario, replicate, method) {
  effects <- extract_marker_effects(result)
  data.frame(study = "02_prediction", scenario = scenario,
    replicate = as.integer(replicate), method = method,
    marker = rep(rownames(effects), times = ncol(effects)),
    trait = rep(colnames(effects), each = nrow(effects)),
    estimate = as.numeric(effects), quantity = "posterior_mean_effect",
    stringsAsFactors = FALSE)
}

prediction_marker_table <- function(result, scenario, replicate, method) {
  effects <- extract_marker_effects(result)
  probabilities <- extract_marker_probabilities(result)$posterior_inclusion
  out <- data.frame(study = "02_prediction", scenario = scenario,
    replicate = as.integer(replicate), method = method,
    marker = rep(rownames(effects), times = ncol(effects)),
    trait = rep(colnames(effects), each = nrow(effects)),
    posterior_mean_effect = as.numeric(effects), stringsAsFactors = FALSE)
  if (!is.null(probabilities)) {
    probabilities <- align_traits(align_markers(probabilities,
      rownames(effects)), colnames(effects))
    out$posterior_inclusion_probability <- as.numeric(probabilities)
  }
  out
}
