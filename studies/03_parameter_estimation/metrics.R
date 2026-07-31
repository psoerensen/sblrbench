.study03_model_specification <- function(architecture, prior_class) {
  matched <- (architecture == "sparse_homogeneous" & prior_class == "BayesC") |
    (architecture == "sparse_mixture" & prior_class == "BayesR")
  ifelse(matched, "matched", "misspecified")
}

.study03_recovery <- function(summary, truth, method, registry, tolerance) {
  x <- merge(summary, truth, by = "estimand_id", all.x = TRUE, sort = FALSE)
  x$status <- x$status.x
  x$reason <- ifelse(nzchar(x$reason.x), x$reason.x, x$reason.y)
  x$study <- "03_parameter_estimation"; x$method_label <- method$id
  x$prior_class <- method$prior_class
  x$model_specification <- .study03_model_specification(x$architecture, method$prior_class)
  x$estimand_label <- registry$label[match(x$estimand_id, registry$estimand_id)]
  x$bias <- x$posterior_mean - x$truth; x$absolute_error <- abs(x$bias); x$squared_error <- x$bias^2
  allowed <- registry$relative_error_allowed[match(x$estimand_id, registry$estimand_id)] & abs(x$truth) > tolerance
  x$relative_error <- ifelse(allowed, x$bias / x$truth, NA_real_)
  x$covered_95 <- x$lower_95 <= x$truth & x$truth <= x$upper_95
  x$interval_width_95 <- x$upper_95 - x$lower_95
  x[, c("study", "architecture", "replicate", "method", "method_label", "prior_class",
    "model_specification", "estimand_id", "estimand_label", "truth", "truth_type",
    "posterior_mean", "posterior_sd", "posterior_median", "lower_95", "upper_95",
    "bias", "absolute_error", "squared_error", "relative_error", "covered_95",
    "interval_width_95", "n_posterior_draws", "chain_count", "draws_per_chain",
    "status", "reason")]
}

.study03_complete_summary <- function(summary, registry, method) {
  if (!nrow(summary) && !"estimand_id" %in% names(summary)) summary <- data.frame(
    method = character(), estimand_id = character(), posterior_mean = numeric(),
    posterior_sd = numeric(), posterior_median = numeric(), lower_95 = numeric(),
    upper_95 = numeric(), n_posterior_draws = integer(), chain_count = integer(),
    draws_per_chain = integer(), mcse_mean = numeric(), status = character(),
    reason = character(), stringsAsFactors = FALSE)
  missing <- setdiff(registry$estimand_id, summary$estimand_id)
  if (!length(missing)) return(summary)
  absent <- data.frame(method = method, estimand_id = missing,
    posterior_mean = NA_real_, posterior_sd = NA_real_, posterior_median = NA_real_,
    lower_95 = NA_real_, upper_95 = NA_real_, n_posterior_draws = 0L,
    chain_count = 4L, draws_per_chain = 0L, mcse_mean = NA_real_,
    status = "unavailable", reason = "aligned draw-level posterior component is not retained",
    stringsAsFactors = FALSE)
  rbind(summary, absent)
}

.study03_pair_definitions <- function() data.frame(
  comparison_id = c("csr_vs_bed_bayesc", "csr_vs_bed_bayesr", "bayesr_vs_bayesc_bed", "sbayesr_vs_sbayesc_csr"),
  focal_method = c("st_csr_sbayesc", "st_csr_sbayesr", "st_bed_bayesr", "st_csr_sbayesr"),
  comparison_method = c("st_bed_bayesc", "st_bed_bayesr", "st_bed_bayesc", "st_csr_sbayesc"),
  comparison_type = c("BED-versus-CSR", "BED-versus-CSR", "prior-class", "prior-class"), stringsAsFactors = FALSE)

.study03_paired <- function(x) {
  keys <- c("architecture", "replicate", "estimand_id"); defs <- .study03_pair_definitions(); out <- list()
  if (anyDuplicated(x[c(keys, "method")])) stop("Duplicate recovery rows prevent pairing.", call. = FALSE)
  for (i in seq_len(nrow(defs))) {
    a <- x[x$method == defs$focal_method[i], ]; b <- x[x$method == defs$comparison_method[i], ]
    z <- merge(a, b, by = keys, suffixes = c("_focal", "_comparison"), all = FALSE)
    if (!nrow(z)) stop("Incomplete parameter comparison rows.", call. = FALSE)
    complete <- is.finite(z$posterior_mean_focal) & is.finite(z$posterior_mean_comparison)
    out[[i]] <- data.frame(z[keys], defs[i, ], estimate_difference = z$posterior_mean_focal - z$posterior_mean_comparison,
      absolute_error_difference = z$absolute_error_focal - z$absolute_error_comparison,
      interval_width_difference = z$interval_width_95_focal - z$interval_width_95_comparison,
      coverage_agreement = z$covered_95_focal == z$covered_95_comparison,
      complete_pair = complete,
      reason = ifelse(complete, "", "one or both method-specific estimands unavailable"),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

.study03_aggregate <- function(x) do.call(rbind, lapply(split(x,
  interaction(x$architecture, x$method, x$model_specification, x$estimand_id, drop = TRUE)), function(z) {
    ok <- z$status == "ok" & is.finite(z$posterior_mean)
    v <- z[ok, , drop = FALSE]
    data.frame(
    architecture = z$architecture[1], method = z$method[1], model_specification = z$model_specification[1],
    estimand_id = z$estimand_id[1], replicate_count = length(unique(z$replicate)),
    successful_replicates = nrow(v), mean_truth = if (nrow(v)) mean(v$truth) else NA_real_,
    mean_posterior_mean = if (nrow(v)) mean(v$posterior_mean) else NA_real_,
    mean_bias = if (nrow(v)) mean(v$bias) else NA_real_,
    rmse = if (nrow(v)) sqrt(mean(v$squared_error)) else NA_real_,
    mae = if (nrow(v)) mean(v$absolute_error) else NA_real_,
    observed_coverage_count = if (nrow(v)) sum(v$covered_95) else 0L,
    observed_coverage_proportion = if (nrow(v)) mean(v$covered_95) else NA_real_,
    coverage_note = "coarse observed coverage; each miss changes five-replicate coverage by 20 percentage points",
    empirical_coverage_95 = if (nrow(v) > 1L) mean(v$covered_95) else NA_real_,
    coverage_observation = if (nrow(v) == 1L) v$covered_95[1] else NA,
    mean_interval_width_95 = if (nrow(v)) mean(v$interval_width_95) else NA_real_,
    minimum = if (nrow(v)) min(v$posterior_mean) else NA_real_,
    maximum = if (nrow(v)) max(v$posterior_mean) else NA_real_,
    failure_count = sum(!ok), stringsAsFactors = FALSE)
  }))
