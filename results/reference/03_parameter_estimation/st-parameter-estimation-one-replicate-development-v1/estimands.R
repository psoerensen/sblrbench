.study03_estimand_registry <- function() {
  all_methods <- paste(c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"), collapse = ";")
  inclusion_trace_methods <- paste(c("st_bed_bayesc", "st_csr_sbayesc", "st_csr_sbayesr"), collapse = ";")
  data.frame(
    estimand_id = c("causal_proportion", "effect_variance",
      "total_marker_effect_variance", "genetic_variance", "residual_variance", "heritability"),
    label = c("Causal proportion", "Nonzero effect variance", "Total marker-effect variance",
      "Genetic variance", "Residual variance", "Heritability"),
    description = c("Proportion of markers with nonzero effects",
      "Mean squared nonzero marker effect", "Sum of squared marker effects",
      "Sample variance of genetic values", "Sample variance of residuals",
      "Genetic variance divided by genetic plus residual variance"),
    unit_or_scale = c("proportion", "scaled-effect squared", "scaled-effect squared",
      "phenotype variance", "phenotype variance", "proportion"),
    truth_type = c("realized_quantity", "realized_quantity", "realized_quantity",
      "realized_quantity", "realized_quantity", "realized_quantity"),
    posterior_source = c("pi_trace", "vbs", "vbs*pi_trace*marker_count", "vgs", "ves", "vgs/(vgs+ves)"),
    truth_source = c("causal_count/marker_count", "mean(causal_effects^2)", "sum(marker_effects^2)",
      "var(genetic_values)", "var(residuals)", "var(G)/(var(G)+var(E))"),
    available_methods = c(inclusion_trace_methods, all_methods, inclusion_trace_methods,
      all_methods, all_methods, all_methods),
    primary = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE),
    bounded = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    lower_bound = c(0, 0, 0, 0, 0, 0), upper_bound = c(1, Inf, Inf, Inf, Inf, 1),
    relative_error_allowed = TRUE, stringsAsFactors = FALSE)
}

.study03_validate_registry <- function(x) {
  required <- c("estimand_id", "label", "description", "unit_or_scale", "truth_type",
    "posterior_source", "truth_source", "available_methods", "primary", "bounded",
    "lower_bound", "upper_bound", "relative_error_allowed")
  if (!all(required %in% names(x)) || anyDuplicated(x$estimand_id) ||
      any(!x$truth_type %in% c("generating_parameter", "realized_quantity")) ||
      any(x$lower_bound > x$upper_bound)) stop("Invalid Study 03 estimand registry.", call. = FALSE)
  x
}
