.study04_registry <- function() data.frame(
  estimand = c("effect_variance", "genetic_variance", "residual_variance", "heritability",
    "causal_proportion", "total_marker_effect_variance"),
  source = c("vbs", "vgs", "ves", "vgs/(vgs+ves)", "pi_trace", "vbs*pi_trace*marker_count"),
  classification = c(rep("core", 4), "unavailable", "unavailable"),
  required = c(rep(TRUE, 4), FALSE, FALSE), lower_bound = 0,
  upper_bound = c(Inf, Inf, Inf, 1, 1, Inf), stringsAsFactors = FALSE)

.study04_thresholds <- function(config) data.frame(
  diagnostic = c("rhat", "ess_bulk", "ess_tail", "relative_mcse", "chain_count", "standardized_mean_shift"),
  operator = c("<=", ">=", ">=", "<=", "==", "<="),
  threshold = unlist(config$thresholds, use.names = FALSE), stringsAsFactors = FALSE)
