.study02_validate_promotion_tables <- function(config, metrics, paired,
                                                computation, status,
                                                simulations, manifest) {
  expected_methods <- c("st_bed_bayesc", "st_bed_bayesr",
    "st_csr_sbayesc", "st_csr_sbayesr")
  architectures <- c("sparse_homogeneous", "sparse_mixture")
  if (!identical(config$methods, expected_methods) ||
      any(grepl("^mt_", unique(c(metrics$method, computation$method, status$method))))) {
    stop("Promotion refuses active MT rows or an invalid active method set.", call. = FALSE)
  }
  replicate_count <- as.integer(manifest$replicate_count)
  expected <- expand.grid(architecture = architectures,
    replicate = seq_len(replicate_count), method = expected_methods,
    stringsAsFactors = FALSE)
  key <- function(x) paste(x$architecture, x$replicate, x$method, sep = "|")
  if (nrow(status) != nrow(expected) || !setequal(key(status), key(expected)) ||
      any(status$status != "ok")) stop("Promotion requires the complete successful active fit grid.", call. = FALSE)
  if (nrow(computation) != nrow(expected) || !setequal(key(computation), key(expected)) ||
      any(computation$status != "ok") || any(!is.finite(computation$runtime))) stop("Computational summaries are incomplete or failed.", call. = FALSE)
  required_metrics <- c("prediction_correlation", "prediction_mse",
    "prediction_nmse", "phenotype_prediction_correlation",
    "prediction_calibration_intercept", "prediction_calibration_slope",
    "effect_rmse")
  observed <- unique(metrics[, c("architecture", "replicate", "method", "metric")])
  expected_metric_n <- nrow(expected) * length(required_metrics)
  if (nrow(observed) != expected_metric_n ||
      !setequal(unique(metrics$metric), required_metrics) ||
      any(metrics$status != "ok") || any(!is.finite(metrics$value))) stop("Required prediction metrics are incomplete or invalid.", call. = FALSE)
  if (nrow(simulations) != length(architectures) * replicate_count ||
      !setequal(simulations$architecture, architectures) ||
      any(simulations$causal_count != config$simulation$n_causal) ||
      any(!simulations$oracle_ok) ||
      any(abs(simulations$realized_h2 - config$simulation$h2) > 1e-10)) stop("Simulation validation is incomplete.", call. = FALSE)
  required_comparisons <- c("bayesr_vs_bayesc_bed", "sbayesr_vs_sbayesc_csr",
    "csr_vs_bed_bayesc", "csr_vs_bed_bayesr")
  if (!setequal(unique(paired$comparison_id), required_comparisons) ||
      any(!paired$complete_pair) || any(!is.finite(paired$advantage))) stop("Paired comparisons are incomplete.", call. = FALSE)
  if (!identical(manifest$task, "single_trait_prediction") ||
      isTRUE(manifest$multitrait$enabled) ||
      !identical(unname(manifest$active_methods), expected_methods)) stop("Manifest does not describe the active single-trait benchmark.", call. = FALSE)
  invisible(TRUE)
}
