# Shared, table-only presentation helpers for benchmark reports.

sblrbench_method_levels <- function() c(
  "st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"
)

sblrbench_method_labels <- function() c(
  st_bed_bayesc = "ST-BED BayesC", st_bed_bayesr = "ST-BED BayesR",
  st_csr_sbayesc = "ST-CSR SBayesC", st_csr_sbayesr = "ST-CSR SBayesR"
)

sblrbench_method_colours <- function() c(
  st_bed_bayesc = "#0072B2", st_bed_bayesr = "#D55E00",
  st_csr_sbayesc = "#009E73", st_csr_sbayesr = "#CC79A7"
)

sblrbench_method_shapes <- function() c(
  st_bed_bayesc = 16, st_bed_bayesr = 16,
  st_csr_sbayesc = 17, st_csr_sbayesr = 17
)

sblrbench_method_linetypes <- function() c(
  st_bed_bayesc = "solid", st_bed_bayesr = "solid",
  st_csr_sbayesc = "dashed", st_csr_sbayesr = "dashed"
)

sblrbench_architecture_levels <- function() c("sparse_homogeneous", "sparse_mixture")
sblrbench_architecture_labels <- function() c(
  sparse_homogeneous = "Sparse homogeneous",
  sparse_mixture = "Sparse variance mixture"
)

sblrbench_method_factor <- function(x) factor(x, levels = sblrbench_method_levels(),
  labels = unname(sblrbench_method_labels()[sblrbench_method_levels()]))
sblrbench_architecture_factor <- function(x) factor(x,
  levels = sblrbench_architecture_levels(),
  labels = unname(sblrbench_architecture_labels()[sblrbench_architecture_levels()]))

sblrbench_model_specification_labels <- function() c(
  matched = "Matched prior class", misspecified = "Misspecified prior class"
)

format_sblrbench_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits, big.mark = ","))
}
format_sblrbench_proportion <- function(x, digits = 3L) format_sblrbench_number(x, digits)
format_sblrbench_runtime <- function(x) {
  ifelse(is.na(x), "—", ifelse(x < 60, sprintf("%.1f s", x), sprintf("%.1f min", x / 60)))
}
format_sblrbench_interval <- function(mean, lower, upper, digits = 3L) paste0(
  format_sblrbench_number(mean, digits), " [",
  format_sblrbench_number(lower, digits), ", ",
  format_sblrbench_number(upper, digits), "]"
)
format_sblrbench_status <- function(x) {
  labels <- c(ok = "Passed", pass = "Passed", complete = "Complete",
    failed = "Failed", fail = "Failed", unavailable = "Unavailable",
    indeterminate = "Indeterminate")
  out <- unname(labels[as.character(x)]); out[is.na(out)] <- as.character(x)[is.na(out)]; out
}

theme_sblrbench <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E5E7EB", linewidth = .3),
      plot.title.position = "plot", plot.title = ggplot2::element_text(face = "bold"),
      plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 0),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", colour = NA),
      panel.spacing = grid::unit(.8, "lines"), legend.position = "bottom",
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

sblrbench_method_scales <- function() list(
  ggplot2::scale_colour_manual(values = sblrbench_method_colours(),
    breaks = sblrbench_method_levels(), labels = sblrbench_method_labels()),
  ggplot2::scale_shape_manual(values = sblrbench_method_shapes(),
    breaks = sblrbench_method_levels(), labels = sblrbench_method_labels()),
  ggplot2::scale_linetype_manual(values = sblrbench_method_linetypes(),
    breaks = sblrbench_method_levels(), labels = sblrbench_method_labels())
)

prepare_sblrbench_replicates <- function(data, group_cols, value_col = "value",
                                         replicate_col = "replicate") {
  needed <- c(group_cols, value_col, replicate_col)
  if (!all(needed %in% names(data))) stop("Replicate data are missing required columns.", call. = FALSE)
  observations <- data[needed]
  key <- interaction(observations[group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(observations, key)
  summary <- do.call(rbind, lapply(pieces, function(z) {
    values <- z[[value_col]]
    row <- z[1L, group_cols, drop = FALSE]
    row$n_replicates <- length(unique(z[[replicate_col]]))
    row$mean <- mean(values); row$sd <- if (row$n_replicates > 1L) stats::sd(values) else NA_real_
    row$minimum <- min(values); row$maximum <- max(values); row
  }))
  rownames(summary) <- NULL
  list(observations = observations, summary = summary)
}

.sblrbench_script_text <- function(path) {
  if (length(path) != 1L || !file.exists(path)) stop("Capsule script does not exist: ", path, call. = FALSE)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

display_capsule_script <- function(path) {
  cat("```r\n", .sblrbench_script_text(path), "\n```\n", sep = "")
  invisible(path)
}

display_capsule_script_collapsed <- function(path, summary = "Show developer contract smoke test") {
  cat("<details><summary>", summary, "</summary>\n\n```r\n",
    .sblrbench_script_text(path), "\n```\n\n</details>\n", sep = "")
  invisible(path)
}

.benchmark_collapse_values <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  paste(format(x, trim = TRUE, scientific = FALSE), collapse = "; ")
}

.benchmark_scenario_name <- function(data) {
  if ("scenario" %in% names(data)) "scenario" else if (
    "architecture" %in% names(data)) "architecture" else
      stop("Data require a scenario or architecture column.", call. = FALSE)
}

.benchmark_require_columns <- function(data, columns, label) {
  if (!is.data.frame(data)) stop(label, " must be a data frame.", call. = FALSE)
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(label, " is missing columns: ",
    paste(missing, collapse = ", "), ".", call. = FALSE)
  invisible(data)
}

.benchmark_ok_rows <- function(data) {
  if (!"status" %in% names(data)) rep(TRUE, nrow(data)) else
    !is.na(data$status) & data$status == "ok"
}

#' Summarize benchmark data design
#'
#' @param spec A validated benchmark specification.
#' @param profile A supported profile name.
#' @return A one-row data frame describing data preparation without loading it.
#' @export
benchmark_data_summary <- function(spec, profile = "benchmark") {
  validate_benchmark_spec(spec)
  resolved <- resolve_benchmark_profile(spec, profile)
  split <- spec$split
  sample_count <- spec$validation$benchmark_sample_count
  training_count <- if (is.null(split)) sample_count else
    spec$validation$benchmark_training_count
  test_count <- if (is.null(split)) NA_integer_ else
    spec$validation$benchmark_test_count
  data.frame(study = spec$study, profile = profile,
    data_source = paste0(spec$data$source, " at ",
      spec$data$example_data$commit),
    sample_count = sample_count, training_count = training_count,
    test_count = test_count,
    split = if (is.null(split)) "all samples used for estimation" else
      paste0(format(100 * split$train_fraction, trim = TRUE),
        "/", format(100 * (1 - split$train_fraction), trim = TRUE),
        " train/test; seed ", split$seed),
    marker_count = spec$validation$benchmark_marker_count,
    marker_selection = paste0("chromosome ", spec$data$chromosome,
      "; QC: MAF >= ", spec$markers$qc$excludeMAF,
      ", missingness <= ", spec$markers$qc$excludeMISS,
      ", canonical QC order retained"),
    alignment_policy = "strict sample, marker, and trait identifiers; requested order preserved",
    preparation_rule = if (is.null(split))
      "full-sample genotype scaling, sparse LD, and summary statistics" else
      "training-only allele frequencies, scaling, sparse LD, and summary statistics",
    replicate_count = resolved$replicate_count, stringsAsFactors = FALSE)
}

#' Tabulate benchmark scenarios
#'
#' @inheritParams benchmark_data_summary
#' @return A tidy data frame with one row per scenario.
#' @export
benchmark_scenario_table <- function(spec, profile = "benchmark") {
  validate_benchmark_spec(spec)
  replicate_count <- resolve_benchmark_profile(spec, profile)$replicate_count
  do.call(rbind, lapply(names(spec$scenarios), function(id) {
    scenario <- spec$scenarios[[id]]
    data.frame(scenario = id,
      target_heritability = spec$controls$simulation$h2,
      causal_markers = spec$controls$simulation$n_causal,
      effect_distribution = scenario$effect_distribution,
      mixture_probabilities = .benchmark_collapse_values(
        scenario$mixture_prob),
      mixture_multipliers = .benchmark_collapse_values(
        scenario$mixture_var),
      replicate_count = replicate_count, stringsAsFactors = FALSE)
  }))
}

#' Tabulate benchmark coordinates and seeds
#'
#' @inheritParams benchmark_data_summary
#' @return A tidy coordinate table with chain seeds collapsed for display.
#' @export
benchmark_coordinate_table <- function(spec, profile = "benchmark") {
  coordinates <- benchmark_seeds(spec, profile)
  data.frame(scenario = coordinates$scenario,
    replicate = coordinates$replicate, method = coordinates$method,
    simulation_seed = coordinates$simulation_seed,
    method_seed = coordinates$fit_seed,
    chain_seeds = vapply(coordinates$chain_seeds,
      .benchmark_collapse_values, character(1)), stringsAsFactors = FALSE)
}

#' Tabulate benchmark methods and resolved controls
#'
#' @inheritParams benchmark_data_summary
#' @return A tidy data frame with one row per method.
#' @export
benchmark_method_table <- function(spec, profile = "benchmark") {
  coordinates <- benchmark_seeds(spec, profile)
  do.call(rbind, lapply(names(spec$methods), function(id) {
    method <- spec$methods[[id]]
    coordinate <- coordinates[coordinates$method == id, , drop = FALSE][1L, ]
    controls <- benchmark_method_controls(spec, id, profile,
      coordinate$fit_seed, coordinate$chain_seeds[[1L]])
    data.frame(method = id, label = method$label,
      interface = method$interface, model = method$native_method,
      representation = method$representation,
      prior_class = method$prior_class, nburn = controls$nburn,
      nit = controls$nit, nthin = controls$nthin,
      nchains = controls$nchains, ncores = controls$ncores,
      inclusion_prior = if (identical(spec$task,"finemapping"))
        "backend default (recorded in fit metadata)" else if (!is.null(controls$pi_init))
        paste0("pi_init=", controls$pi_init) else
        paste0("pi=", .benchmark_collapse_values(controls$pi)),
      mixture_multipliers = .benchmark_collapse_values(controls$mixture_var),
      update_flags = "backend defaults, recorded in fit metadata",
      scheduler_policy = if (identical(method$representation, "BED"))
        "validated sblr BED scheduler defaults" else "not applicable",
      ld_policy = if (identical(method$representation, "CSR"))
        paste0("max variants ", spec$data$sparse_ld$max_distance_variants,
          "; r2 >= ", spec$data$sparse_ld$r2_threshold,
          "; block size ", spec$data$sparse_ld$block_size) else
        "individual-level BED genotype operator",
      stringsAsFactors = FALSE)
  }))
}

#' Tabulate parameter-estimation estimands
#'
#' @param spec A parameter-estimation benchmark specification.
#' @return The authoritative estimand definitions in display-ready form.
#' @export
benchmark_estimand_table <- function(spec) {
  validate_benchmark_spec(spec)
  if (is.null(spec$estimands))
    stop("The benchmark specification does not define estimands.",
      call. = FALSE)
  data.frame(estimand = spec$estimands$estimand_id,
    label = spec$estimands$label,
    formula_or_fit_field = spec$estimands$posterior_source,
    required_fit_fields = vapply(spec$estimands$posterior_source,
      function(x) paste(unique(regmatches(x,
        gregexpr("vbs|vgs|ves|pi_trace", x))[[1L]]), collapse = "; "),
      character(1)),
    truth_definition = spec$estimands$truth_source,
    available_methods = spec$estimands$available_methods,
    primary = spec$estimands$primary, stringsAsFactors = FALSE)
}

#' Inventory benchmark output paths
#'
#' @param results A result index returned by `run_benchmark()`.
#' @return A data frame of named output paths and current existence status.
#' @export
benchmark_output_inventory <- function(results) {
  if (!is.list(results) || !is.list(results$paths))
    stop("results must contain a named paths list.", call. = FALSE)
  paths <- results$paths
  keep <- vapply(paths, function(x) is.character(x) && length(x) == 1L,
    logical(1))
  paths <- paths[keep]
  data.frame(output = names(paths), path = unname(unlist(paths)),
    exists = file.exists(unname(unlist(paths))), stringsAsFactors = FALSE)
}

#' Summarize benchmark runtime
#'
#' @param runtime Tidy runtime rows from `run_benchmark()`.
#' @return A data frame summarized by scenario and method.
#' @export
benchmark_runtime_summary <- function(runtime) {
  .benchmark_require_columns(runtime, c("study", "scenario", "method",
    "elapsed_seconds", "status"), "runtime")
  key <- interaction(runtime$scenario, runtime$method, drop = TRUE,
    lex.order = TRUE)
  out <- do.call(rbind, lapply(split(runtime, key), function(x) {
    values <- x$elapsed_seconds[x$status == "ok" &
      is.finite(x$elapsed_seconds)]
    data.frame(study = x$study[[1L]], scenario = x$scenario[[1L]],
      method = x$method[[1L]], successful_fits = length(values),
      mean_seconds = if (length(values)) mean(values) else NA_real_,
      median_seconds = if (length(values)) stats::median(values) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

.benchmark_plot_data <- function(data) {
  scenario <- .benchmark_scenario_name(data)
  data$scenario_display <- sblrbench_architecture_factor(data[[scenario]])
  data$method_display <- sblrbench_method_factor(data$method)
  data
}

#' Plot prediction performance metrics
#'
#' @param metrics Tidy prediction metric rows.
#' @param metric_ids Metrics to include.
#' @return A `ggplot2` plot.
#' @export
plot_prediction_metrics <- function(metrics, metric_ids = c(
    "prediction_correlation", "prediction_nmse")) {
  .benchmark_require_columns(metrics, c("method", "metric", "value"),
    "metrics")
  data <- metrics[metrics$metric %in% metric_ids &
    .benchmark_ok_rows(metrics), , drop = FALSE]
  if (!nrow(data)) stop("No requested prediction metrics are available.",
    call. = FALSE)
  data <- .benchmark_plot_data(data)
  ggplot2::ggplot(data, ggplot2::aes(method_display, value,
      colour = method)) +
    ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .08),
      alpha = .7, show.legend = FALSE) +
    ggplot2::facet_grid(metric ~ scenario_display, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Metric value",
      title = "Prediction performance") + theme_sblrbench()
}

#' Plot prediction calibration
#'
#' @inheritParams plot_prediction_metrics
#' @return A `ggplot2` plot.
#' @export
plot_prediction_calibration <- function(metrics) {
  plot_prediction_metrics(metrics, c("prediction_calibration_intercept",
    "prediction_calibration_slope")) +
    ggplot2::labs(title = "Prediction calibration")
}

#' Plot prediction effect recovery
#'
#' @inheritParams plot_prediction_metrics
#' @return A `ggplot2` plot.
#' @export
plot_effect_recovery <- function(metrics) {
  plot_prediction_metrics(metrics, c("effect_rmse",
    "genetic_value_correlation", "genetic_value_rmse")) +
    ggplot2::labs(title = "Effect and genetic-value recovery")
}

#' Plot parameter recovery against truth
#'
#' @param estimates Tidy parameter-estimate rows.
#' @param estimand_id One estimand identifier.
#' @return A `ggplot2` plot.
#' @export
plot_parameter_recovery <- function(estimates, estimand_id) {
  .benchmark_require_columns(estimates, c("method", "estimand_id", "truth",
    "posterior_mean"), "estimates")
  data <- estimates[estimates$estimand_id == estimand_id &
    .benchmark_ok_rows(estimates), , drop = FALSE]
  if (!nrow(data)) stop("No estimates are available for `", estimand_id,
    "`.", call. = FALSE)
  data <- .benchmark_plot_data(data)
  ggplot2::ggplot(data, ggplot2::aes(truth, posterior_mean,
      colour = method, shape = method)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2,
      colour = "grey45") + ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(~scenario_display, scales = "free") +
    sblrbench_method_scales() +
    ggplot2::labs(x = "Realized truth", y = "Posterior mean",
      title = paste("Recovery:", estimand_id), colour = NULL, shape = NULL) +
    theme_sblrbench()
}

#' Plot parameter bias
#'
#' @param estimates Tidy parameter-estimate rows.
#' @return A `ggplot2` plot.
#' @export
plot_parameter_bias <- function(estimates) {
  .benchmark_require_columns(estimates,
    c("method", "estimand_id", "bias"), "estimates")
  data <- estimates[.benchmark_ok_rows(estimates), , drop = FALSE]
  data <- .benchmark_plot_data(data)
  ggplot2::ggplot(data, ggplot2::aes(method_display, bias,
      colour = method)) + ggplot2::geom_hline(yintercept = 0,
      linetype = 2, colour = "grey45") +
    ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .08),
      alpha = .65, show.legend = FALSE) +
    ggplot2::facet_grid(estimand_id ~ scenario_display, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Posterior mean minus truth",
      title = "Parameter bias") + theme_sblrbench()
}

#' Plot active or component probabilities
#'
#' @param estimates Tidy parameter-estimate rows.
#' @return A `ggplot2` plot.
#' @export
plot_component_probabilities <- function(estimates) {
  .benchmark_require_columns(estimates,
    c("method", "estimand_id", "truth", "posterior_mean"), "estimates")
  ids <- grepl("causal_proportion|active_probability|component.*probability",
    estimates$estimand_id)
  data <- estimates[ids & .benchmark_ok_rows(estimates), , drop = FALSE]
  if (!nrow(data)) stop("No active/component probability estimates are available.",
    call. = FALSE)
  data <- .benchmark_plot_data(data)
  ggplot2::ggplot(data, ggplot2::aes(method_display, posterior_mean,
      colour = method)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .08),
      alpha = .75, show.legend = FALSE) +
    ggplot2::facet_grid(estimand_id ~ scenario_display, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Posterior mean probability",
      title = "Active-marker probability recovery") + theme_sblrbench()
}

#' Plot benchmark runtime
#'
#' @param runtime Tidy runtime rows.
#' @return A `ggplot2` plot.
#' @export
plot_benchmark_runtime <- function(runtime) {
  .benchmark_require_columns(runtime,
    c("method", "elapsed_seconds"), "runtime")
  data <- runtime[.benchmark_ok_rows(runtime), , drop = FALSE]
  if (!nrow(data)) stop("No successful runtime rows are available.",
    call. = FALSE)
  data <- .benchmark_plot_data(data)
  ggplot2::ggplot(data, ggplot2::aes(method_display, elapsed_seconds,
      colour = method)) +
    ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = .08),
      alpha = .65, show.legend = FALSE) +
    ggplot2::facet_wrap(~scenario_display, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Elapsed seconds", title = "Runtime") +
    theme_sblrbench()
}

#' Plot causal-marker posterior inclusion probabilities
#' @param marker_results Tidy fine-mapping marker results.
#' @return A `ggplot2` object.
#' @export
plot_causal_marker_pip <- function(marker_results) {
  .benchmark_require_columns(marker_results,c("method","causal",
    "posterior_inclusion_probability"),"marker_results")
  data <- marker_results[marker_results$causal,,drop=FALSE]
  data$method_display <- sblrbench_method_factor(data$method)
  ggplot2::ggplot(data,ggplot2::aes(method_display,
      posterior_inclusion_probability,colour=method)) +
    ggplot2::geom_boxplot(outlier.shape=NA,show.legend=FALSE) +
    ggplot2::geom_point(position=ggplot2::position_jitter(width=.08),
      alpha=.6,show.legend=FALSE) +
    ggplot2::labs(x=NULL,y="Causal-marker PIP",title="Causal-marker recovery") +
    theme_sblrbench()
}

#' Plot causal-marker ranks
#' @inheritParams plot_causal_marker_pip
#' @return A `ggplot2` object.
#' @export
plot_causal_marker_rank <- function(marker_results) {
  .benchmark_require_columns(marker_results,c("method","causal","causal_rank"),
    "marker_results")
  data <- marker_results[marker_results$causal,,drop=FALSE]
  data$method_display <- sblrbench_method_factor(data$method)
  ggplot2::ggplot(data,ggplot2::aes(method_display,causal_rank,colour=method)) +
    ggplot2::geom_boxplot(outlier.shape=NA,show.legend=FALSE) +
    ggplot2::scale_y_log10() + ggplot2::labs(x=NULL,y="Causal-marker rank",
      title="Causal-marker ranking (lower is better)") + theme_sblrbench()
}

#' Plot credible-set sizes
#' @param credible_sets Tidy credible-set results.
#' @return A `ggplot2` object.
#' @export
plot_credible_set_size <- function(credible_sets) {
  .benchmark_require_columns(credible_sets,c("method","set_size"),
    "credible_sets")
  credible_sets$method_display <- sblrbench_method_factor(credible_sets$method)
  ggplot2::ggplot(credible_sets,ggplot2::aes(method_display,set_size,
      colour=method)) + ggplot2::geom_boxplot(show.legend=FALSE) +
    ggplot2::labs(x=NULL,y="Markers",title="Credible-set size") +
    theme_sblrbench()
}

#' Plot credible-set coverage
#' @inheritParams plot_credible_set_size
#' @return A `ggplot2` object.
#' @export
plot_credible_set_coverage <- function(credible_sets) {
  .benchmark_require_columns(credible_sets,c("method","exact_covered",
    "ld_proxy_covered"),"credible_sets")
  long <- rbind(transform(credible_sets,coverage="Exact",covered=exact_covered),
    transform(credible_sets,coverage="LD proxy",covered=ld_proxy_covered))
  aggregate <- stats::aggregate(covered ~ method + coverage,long,mean)
  aggregate$method_display <- sblrbench_method_factor(aggregate$method)
  ggplot2::ggplot(aggregate,ggplot2::aes(method_display,covered,fill=coverage)) +
    ggplot2::geom_col(position="dodge") +
    ggplot2::labs(x=NULL,y="Coverage fraction",fill=NULL,
      title="Credible-set coverage") + theme_sblrbench()
}

#' Plot posterior-inclusion-probability calibration
#' @inheritParams plot_causal_marker_pip
#' @return A `ggplot2` object.
#' @export
plot_pip_calibration <- function(marker_results) {
  .benchmark_require_columns(marker_results,c("posterior_inclusion_probability",
    "causal","method"),"marker_results")
  marker_results$bin <- cut(marker_results$posterior_inclusion_probability,
    breaks=seq(0,1,by=.1),include.lowest=TRUE)
  data <- stats::aggregate(cbind(pip=posterior_inclusion_probability,
    observed=as.numeric(causal)) ~ method + bin,marker_results,mean)
  ggplot2::ggplot(data,ggplot2::aes(pip,observed,colour=method)) +
    ggplot2::geom_abline(slope=1,intercept=0,linetype=2) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    sblrbench_method_scales()[[1L]] + ggplot2::labs(x="Mean PIP",
      y="Observed causal fraction",title="PIP calibration") + theme_sblrbench()
}

.benchmark_convergence_plot_data <- function(data, value) {
  .benchmark_require_columns(data, c("method", "quantity",
    "retained_draw_candidate", value), "convergence data")
  data$method_display <- sblrbench_method_factor(data$method)
  data
}

#' Plot Study 04 R-hat trajectories
#' @param convergence Tidy convergence diagnostics.
#' @param threshold Prespecified maximum R-hat.
#' @return A `ggplot2` object.
#' @export
plot_convergence_rhat <- function(convergence, threshold = 1.01) {
  data <- .benchmark_convergence_plot_data(convergence, "rhat")
  ggplot2::ggplot(data, ggplot2::aes(retained_draw_candidate, rhat,
      colour = quantity, group = quantity)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_grid(burnin_candidate ~ method_display) +
    ggplot2::labs(x = "Retained draws per chain", y = "R-hat",
      colour = "Quantity") + theme_sblrbench()
}

#' Plot Study 04 effective sample sizes
#' @param convergence Tidy convergence diagnostics.
#' @param threshold Prespecified minimum ESS.
#' @return A `ggplot2` object.
#' @export
plot_convergence_ess <- function(convergence, threshold = 400) {
  data <- .benchmark_convergence_plot_data(convergence, "ess_bulk")
  long <- rbind(transform(data, diagnostic = "Bulk ESS", value = ess_bulk),
    transform(data, diagnostic = "Tail ESS", value = ess_tail))
  ggplot2::ggplot(long, ggplot2::aes(retained_draw_candidate, value,
      colour = quantity, group = quantity)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_grid(diagnostic + burnin_candidate ~ method_display,
      scales = "free_y") +
    ggplot2::labs(x = "Retained draws per chain",
      y = "Effective sample size", colour = "Quantity") +
    theme_sblrbench()
}

#' Plot Study 04 relative MCSE trajectories
#' @param convergence Tidy convergence diagnostics.
#' @param threshold Prespecified maximum relative MCSE.
#' @return A `ggplot2` object.
#' @export
plot_convergence_mcse <- function(convergence, threshold = 0.05) {
  data <- .benchmark_convergence_plot_data(convergence, "relative_mcse")
  ggplot2::ggplot(data, ggplot2::aes(retained_draw_candidate, relative_mcse,
      colour = quantity, group = quantity)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_grid(burnin_candidate ~ method_display) +
    ggplot2::labs(x = "Retained draws per chain", y = "Relative MCSE",
      colour = "Quantity") + theme_sblrbench()
}

#' Plot Study 04 burn-in stability
#' @param stability Tidy standardized-mean-shift rows.
#' @param threshold Prespecified maximum standardized shift.
#' @return A `ggplot2` object.
#' @export
plot_convergence_stability <- function(stability, threshold = 0.10) {
  .benchmark_require_columns(stability, c("method", "quantity",
    "burnin_candidate", "standardized_mean_shift"), "stability")
  stability$method_display <- sblrbench_method_factor(stability$method)
  ggplot2::ggplot(stability, ggplot2::aes(burnin_candidate,
      abs(standardized_mean_shift), colour = quantity, group = quantity)) +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~method_display) +
    ggplot2::labs(x = "Burn-in candidate",
      y = "Absolute standardized mean shift", colour = "Quantity") +
    theme_sblrbench()
}
