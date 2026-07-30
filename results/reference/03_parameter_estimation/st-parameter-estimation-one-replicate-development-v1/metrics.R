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
    "interval_width_95", "n_posterior_draws", "status", "reason")]
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
    if (!nrow(z) || anyNA(z$posterior_mean_focal) || anyNA(z$posterior_mean_comparison)) stop("Incomplete parameter pair.", call. = FALSE)
    out[[i]] <- data.frame(z[keys], defs[i, ], estimate_difference = z$posterior_mean_focal - z$posterior_mean_comparison,
      absolute_error_difference = z$absolute_error_focal - z$absolute_error_comparison,
      interval_width_difference = z$interval_width_95_focal - z$interval_width_95_comparison,
      coverage_agreement = z$covered_95_focal == z$covered_95_comparison, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

.study03_aggregate <- function(x) do.call(rbind, lapply(split(x,
  interaction(x$architecture, x$method, x$model_specification, x$estimand_id, drop = TRUE)), function(z) data.frame(
    architecture = z$architecture[1], method = z$method[1], model_specification = z$model_specification[1],
    estimand_id = z$estimand_id[1], replicate_count = length(unique(z$replicate)),
    successful_replicates = sum(z$status == "ok"), mean_truth = mean(z$truth),
    mean_posterior_mean = mean(z$posterior_mean), mean_bias = mean(z$bias),
    rmse = sqrt(mean(z$squared_error)), mae = mean(z$absolute_error),
    empirical_coverage_95 = if (length(unique(z$replicate)) > 1L) mean(z$covered_95) else NA_real_,
    coverage_observation = if (length(unique(z$replicate)) == 1L) z$covered_95[1] else NA,
    mean_interval_width_95 = mean(z$interval_width_95), minimum = min(z$posterior_mean),
    maximum = max(z$posterior_mean), failure_count = sum(z$status != "ok"), stringsAsFactors = FALSE)))
