# Prediction-specific metric selection and summaries. Generic metric kernels
# remain in R/metrics.R.

prediction_metric_table <- function(simulation, result, metrics) {
  out <- evaluate_metrics(simulation, result, metrics)
  names(out)[names(out) == "method_id"] <- "method"
  names(out)[names(out) == "scenario"] <- "scenario"
  out
}

prediction_effect_recovery <- function(simulation, result) {
  estimated <- align_traits(align_markers(extract_marker_effects(result),
    simulation$data$marker_ids), simulation$data$trait_names)
  truth <- simulation$truth$effects
  do.call(rbind, lapply(simulation$data$trait_names, function(trait) {
    data.frame(study = simulation$scenario$study,
      scenario = simulation$scenario$architecture,
      replicate = simulation$scenario$replicate, method = result$method_id,
      trait = trait,
      metric = c("effect_rmse", "effect_correlation"),
      value = c(sqrt(mean((estimated[, trait] - truth[, trait])^2)),
        stats::cor(estimated[, trait], truth[, trait])), status = "ok",
      reason = "", stringsAsFactors = FALSE)
  }))
}

prediction_genetic_value_recovery <- function(simulation, result) {
  predicted <- align_traits(align_samples(result$predictions$genetic_value,
    simulation$data$sample_ids), simulation$data$trait_names)
  truth <- simulation$truth$genetic_values
  do.call(rbind, lapply(simulation$data$trait_names, function(trait) {
    data.frame(study = simulation$scenario$study,
      scenario = simulation$scenario$architecture,
      replicate = simulation$scenario$replicate, method = result$method_id,
      trait = trait,
      metric = c("genetic_value_correlation", "genetic_value_rmse"),
      value = c(stats::cor(predicted[, trait], truth[, trait]),
        sqrt(mean((predicted[, trait] - truth[, trait])^2))), status = "ok",
      reason = "", stringsAsFactors = FALSE)
  }))
}

prediction_comparisons <- function() data.frame(
  comparison_id = c("bayesr_vs_bayesc_bed", "sbayesr_vs_sbayesc_csr",
    "csr_vs_bed_bayesc", "csr_vs_bed_bayesr"),
  focal_method = c("st_bed_bayesr", "st_csr_sbayesr",
    "st_csr_sbayesc", "st_csr_sbayesr"),
  comparison_method = c("st_bed_bayesc", "st_csr_sbayesc",
    "st_bed_bayesc", "st_bed_bayesr"), stringsAsFactors = FALSE)

prediction_paired_metrics <- function(metrics) {
  z <- metrics
  if ("scenario" %in% names(z) && !"architecture" %in% names(z))
    names(z)[names(z) == "scenario"] <- "architecture"
  paired_method_advantages(z, prediction_comparisons())
}

prediction_method_labels <- function() c(
  st_bed_bayesc = "ST-BED BayesC", st_bed_bayesr = "ST-BED BayesR",
  st_csr_sbayesc = "ST-CSR SBayesC", st_csr_sbayesr = "ST-CSR SBayesR")

prediction_comparison_labels <- function() c(
  bayesr_vs_bayesc_bed = "BED BayesR versus BayesC",
  sbayesr_vs_sbayesc_csr = "CSR SBayesR versus SBayesC",
  csr_vs_bed_bayesc = "CSR versus BED BayesC class",
  csr_vs_bed_bayesr = "CSR versus BED BayesR class")

#' Summarize prediction metrics across replicates
#'
#' @param metrics Tidy prediction metric rows.
#' @return A scenario, method, and metric summary data frame.
#' @export
prediction_metric_summary <- function(metrics) {
  architecture <- if ("architecture" %in% names(metrics))
    metrics$architecture else metrics$scenario
  metrics$architecture <- architecture
  key <- interaction(metrics$architecture, metrics$method, metrics$metric,
    drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(metrics, key), function(x) {
    values <- x$value[x$status == "ok"]
    data.frame(architecture = x$architecture[[1L]], method = x$method[[1L]],
      method_label = unname(prediction_method_labels()[x$method[[1L]]]),
      metric = x$metric[[1L]],
      value = if (length(values) == 1L) values else if (length(values))
        mean(values) else NA_real_,
      replicate_count = length(unique(x$replicate)),
      successful_replicates = length(values),
      mean = if (length(values)) mean(values) else NA_real_,
      sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median = if (length(values)) stats::median(values) else NA_real_,
      minimum = if (length(values)) min(values) else NA_real_,
      maximum = if (length(values)) max(values) else NA_real_)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

prediction_paired_summary <- function(paired) {
  key <- interaction(paired$architecture, paired$comparison_id,
    paired$paired_metric, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(paired, key), function(x) {
    values <- x$advantage[x$complete_pair]
    data.frame(architecture = x$architecture[[1L]],
      comparison_id = x$comparison_id[[1L]],
      comparison_label = unname(prediction_comparison_labels()[
        x$comparison_id[[1L]]]), focal_method = x$focal_method[[1L]],
      comparison_method = x$comparison_method[[1L]],
      paired_metric = x$paired_metric[[1L]],
      advantage = if (length(values) == 1L) values else mean(values),
      replicate_count = length(unique(x$replicate)),
      successful_replicates = length(values),
      mean_advantage = if (length(values)) mean(values) else NA_real_,
      sd_advantage = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median_advantage = if (length(values)) stats::median(values) else NA_real_,
      minimum_advantage = if (length(values)) min(values) else NA_real_,
      maximum_advantage = if (length(values)) max(values) else NA_real_,
      complete_pairs = length(values))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

prediction_runtime_summary <- function(runtime) {
  benchmark_runtime_summary(runtime)
}
