for (f in c("annotation_design.R", "simulation.R", "methods.R",
            "chain_extraction.R", "diagnostics.R", "metrics.R", "pilot.R"))
  source(file.path("studies", "06_annotation_models", f), local = TRUE)

.study06_common_targets <- function() list(
  targets::tar_target(study06_config_file,
    file.path("studies", "06_annotation_models", "config.R"), format = "file"),
  targets::tar_target(study06_config,
    source(study06_config_file, local = TRUE)$value),
  targets::tar_target(study06_paths, .study06_paths()),
  targets::tar_target(study06_example_files,
    sblrbench:::benchmark_example_files(study06_paths$data_dir,
      list(example_data=study06_config$example_data))),
  targets::tar_target(study06_base_glist,
    sblrbench:::benchmark_load_glist(study06_paths, study06_example_files)),
  targets::tar_target(study06_filtered_markers,
    sblrbench:::benchmark_filter_markers(study06_base_glist,
      study06_config$chr,study06_config$qc,study06_config$sparse_ld)),
  targets::tar_target(study06_sample_ids,
    sblrbench:::benchmark_selected_ids(study06_base_glist,
      study06_config$sample_limit)),
  targets::tar_target(study06_split, sblrbench::make_prediction_split(
    study06_sample_ids, study06_config$split$train_fraction,
    study06_config$split$seed)),
  targets::tar_target(study06_scaled_genotypes,
    .study06_load_scaled_genotypes(study06_base_glist, study06_config$chr,
      study06_sample_ids, study06_filtered_markers$marker_ids, study06_split)),
  targets::tar_target(study06_working_glist,
    sblrbench:::benchmark_set_training_af(
    sblrbench:::benchmark_set_glist_marker_order(study06_base_glist, study06_config$chr,
      study06_filtered_markers$marker_ids), study06_config$chr,
    study06_filtered_markers$marker_ids,
    study06_scaled_genotypes$allele_frequency)),
  targets::tar_target(study06_ld_bundle, {
    started <- proc.time()[["elapsed"]]
    x <- sblrbench:::benchmark_make_training_ld(study06_working_glist,
      study06_split,
      study06_filtered_markers$marker_ids, list(
        chromosome = study06_config$chr, sparse_ld = study06_config$sparse_ld),
      study06_paths$genotype_output_dir)
    list(Glist = x, elapsed_seconds = proc.time()[["elapsed"]] - started)
  }),
  targets::tar_target(study06_annotation,
    .study06_annotation_design(study06_filtered_markers$marker_ids,
      study06_config)),
  targets::tar_target(study06_alpha,
    .study06_true_alpha(study06_annotation, study06_config)),
  targets::tar_target(study06_annotation_summary,
    .study06_annotation_summary(study06_annotation, study06_alpha,
      study06_config)),
  targets::tar_target(study06_method_specs,
    .study06_method_specs(study06_config), iteration = "list"),
  targets::tar_target(study06_method_spec, study06_method_specs,
    pattern = map(study06_method_specs), iteration = "list")
)

.study06_convergence_targets <- function() c(.study06_common_targets(), list(
  targets::tar_target(study06_convergence_simulation,
    .study06_simulate("informative_annotations", 1L,
      study06_scaled_genotypes$all, study06_split$train_rows,
      study06_annotation, study06_alpha, study06_config)),
  targets::tar_target(study06_convergence_stats,
    .study06_summary_stats(study06_convergence_simulation,
      study06_ld_bundle$Glist, study06_split, study06_config)),
  targets::tar_target(study06_convergence_methods,
    Filter(function(x) isTRUE(x$annotation), study06_method_specs),
    iteration = "list"),
  targets::tar_target(study06_convergence_method,
    study06_convergence_methods,
    pattern = map(study06_convergence_methods), iteration = "list"),
  targets::tar_target(study06_convergence_run,
    .study06_fit(study06_convergence_method, study06_convergence_simulation,
      study06_convergence_stats$stats, study06_ld_bundle$Glist,
      study06_split, study06_annotation, study06_alpha, study06_config,
      phase = "convergence"),
    pattern = map(study06_convergence_method), iteration = "list"),
  targets::tar_target(study06_convergence_checkpoint,
    .study06_fit_checkpoint(study06_convergence_run, "convergence"),
    pattern = map(study06_convergence_run), format = "file"),
  targets::tar_target(study06_convergence_draws_branch,
    {
      study06_convergence_checkpoint
      .study06_extract_chain_draws(study06_convergence_run, study06_annotation)
    },
    pattern = map(study06_convergence_run, study06_convergence_checkpoint),
    iteration = "list"),
  targets::tar_target(study06_convergence_selection_branch,
    .study06_select_recommendation(study06_convergence_draws_branch,
      study06_config),
    pattern = map(study06_convergence_draws_branch), iteration = "list"),
  targets::tar_target(study06_convergence_draws,
    do.call(rbind, study06_convergence_draws_branch)),
  targets::tar_target(study06_convergence_diagnostics,
    do.call(rbind, lapply(study06_convergence_selection_branch, `[[`,
      "diagnostics"))),
  targets::tar_target(study06_convergence_candidates,
    do.call(rbind, lapply(study06_convergence_selection_branch, `[[`,
      "candidates"))),
  targets::tar_target(study06_convergence_recommendations,
    do.call(rbind, lapply(study06_convergence_selection_branch, `[[`,
      "recommendation"))),
  targets::tar_target(study06_convergence_status, do.call(rbind,
    lapply(study06_convergence_run, function(x)
      data.frame(method = x$method$id, scenario = x$scenario,
        replicate = x$replicate, status = x$status, reason = x$reason,
        runtime = x$runtime, chain_count = if (x$status == "ok") 4L else 0L,
        stringsAsFactors = FALSE)))),
  targets::tar_target(study06_convergence_seed_registry, {
    x <- .study06_effective_seed_registry(study06_convergence_run,
      study06_config)
    x[x$scenario == "informative_annotations" & x$replicate == 1L &
      x$method %in% c("st_bed_bayesrc", "st_csr_sbayesrc"), ]
  }),
  targets::tar_target(study06_convergence_output_dir,
    Sys.getenv("SBLR_BENCH_OUTPUT_DIR", file.path(study06_paths$local_dir,
      "convergence_outputs"))),
  targets::tar_target(study06_convergence_draws_file,
    .study06_write_csv(study06_convergence_draws,
      file.path(study06_convergence_output_dir, "scalar_chain_draws.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_diagnostics_file,
    .study06_write_csv(study06_convergence_diagnostics,
      file.path(study06_convergence_output_dir, "convergence_diagnostics.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_candidates_file,
    .study06_write_csv(study06_convergence_candidates,
      file.path(study06_convergence_output_dir, "candidate_settings.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_recommendations_file,
    .study06_write_csv(study06_convergence_recommendations,
      file.path(study06_convergence_output_dir, "method_recommendations.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_status_file,
    .study06_write_csv(study06_convergence_status,
      file.path(study06_convergence_output_dir, "fit_status.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_computational_file,
    .study06_write_csv(study06_convergence_status,
      file.path(study06_convergence_output_dir, "computational_summary.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_seed_file,
    .study06_write_csv(study06_convergence_seed_registry,
      file.path(study06_convergence_output_dir, "seed_registry.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_annotation_file,
    .study06_write_csv(study06_annotation_summary,
      file.path(study06_convergence_output_dir, "annotation_design_summary.csv")),
    format = "file"),
  targets::tar_target(study06_convergence_alpha_file, {
    x <- do.call(rbind, lapply(study06_config$scenarios, function(s)
      data.frame(scenario = s, annotation_name = rownames(study06_alpha[[s]]),
        study06_alpha[[s]], check.names = FALSE)))
    .study06_write_csv(x, file.path(study06_convergence_output_dir,
      "true_alpha.csv"))
  }, format = "file"),
  targets::tar_target(study06_convergence_manifest_file,
    .study06_write_json(list(
      study = "06_annotation_models", task = "annotation_model_convergence_selection",
      benchmark_scope = "two_method_maximum_history_development",
      benchmark_status = if (nrow(study06_convergence_status) == 2L &&
        all(study06_convergence_status$status == "ok") &&
        nrow(study06_convergence_recommendations) == 2L) "complete" else "incomplete",
      method_count = 2L, expected_fit_count = 2L, nchains = 4L,
      maximum_history = study06_config$convergence_pilot,
      mixture_var = study06_config$mixture_var,
      methods = c("st_bed_bayesrc", "st_csr_sbayesrc"),
      installed_sblr_version = as.character(packageVersion("sblr")),
      installed_sblr_commit = packageDescription("sblr")$RemoteSha,
      package_versions = as.list(vapply(c("sblr", "sblrbench", "targets",
        "posterior"), function(x) as.character(packageVersion(x)), "")),
      repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      qgdata = study06_config$example_data,
      seed_registry_file = "seed_registry.csv",
      source_status = "cached_checksum_validated_offline",
      completion_counts = list(successful_fits =
        sum(study06_convergence_status$status == "ok"), failed_fits =
        sum(study06_convergence_status$status != "ok")),
      failures = unname(study06_convergence_status$reason[
        study06_convergence_status$status != "ok"]),
      provenance = list(installed_package_only = TRUE,
        sibling_source_access = "read_only_audit"),
      recommendations = study06_convergence_recommendations,
      validation_status = "complete_grid_validated"),
      file.path(study06_convergence_output_dir, "benchmark_manifest.json")),
    format = "file")
))

.study06_benchmark_targets <- function() c(.study06_common_targets(), list(
  targets::tar_target(study06_specs, .study06_specs(study06_config),
    iteration = "list"),
  targets::tar_target(study06_spec, study06_specs,
    pattern = map(study06_specs), iteration = "list"),
  targets::tar_target(study06_simulation,
    .study06_simulate(study06_spec$scenario, study06_spec$replicate,
      study06_scaled_genotypes$all, study06_split$train_rows,
      study06_annotation, study06_alpha, study06_config),
    pattern = map(study06_spec), iteration = "list"),
  targets::tar_target(study06_bundle, {
    stats <- .study06_summary_stats(study06_simulation,
      study06_ld_bundle$Glist, study06_split, study06_config)
    list(simulation = study06_simulation, stats = stats)
  },
    pattern = map(study06_simulation), iteration = "list"),
  targets::tar_target(study06_run,
    .study06_fit(study06_method_spec, study06_bundle$simulation,
      study06_bundle$stats$stats, study06_ld_bundle$Glist, study06_split,
      study06_annotation, study06_alpha, study06_config, phase = "benchmark"),
    pattern = cross(study06_bundle, study06_method_spec),
    iteration = "list"),
  targets::tar_target(study06_run_checkpoint,
    .study06_fit_checkpoint(study06_run, "benchmark"),
    pattern = map(study06_run), format = "file"),
  targets::tar_target(study06_outputs_branch, {
    study06_run_checkpoint
    metric_started <- proc.time()[["elapsed"]]
    draws <- .study06_extract_chain_draws(study06_run, study06_annotation)
    parameter <- .study06_scalar_estimates(study06_run, draws,
      study06_bundle$simulation)
    annotation <- .study06_annotation_estimates(study06_run, draws,
      study06_bundle$simulation)
    sigma <- .study06_sigma_alpha(study06_run, draws)
    prediction_started <- proc.time()[["elapsed"]]
    prediction <- .study06_prediction_metrics(study06_run,
      study06_bundle$simulation, study06_scaled_genotypes$test,
      study06_split$test_ids)
    prediction_time <- proc.time()[["elapsed"]] - prediction_started
    marker <- .study06_marker_metrics(study06_run, study06_bundle$simulation,
      study06_annotation, study06_config)
    inventory <- .study06_probability_inventory(study06_run,
      study06_annotation, length(study06_config$mixture_var))
    diagnostics <- .study06_diagnostics(draws, 0L,
      max(draws$iteration), study06_config)
    list(draws = draws, prediction = prediction, parameter = parameter,
      annotation = annotation, sigma = sigma, marker = marker$marker,
      probability = marker$probability, enrichment = marker$enrichment,
      probability_inventory = inventory,
      diagnostics = diagnostics,
      computational = data.frame(scenario = study06_run$scenario,
        replicate = study06_run$replicate, method = study06_run$method$id,
        status = study06_run$status, error_message = study06_run$reason,
        preprocessing_time = study06_ld_bundle$elapsed_seconds,
        summary_statistic_time = study06_bundle$stats$elapsed_seconds,
        ld_construction_time = study06_ld_bundle$elapsed_seconds,
        model_fitting_time = study06_run$runtime,
        prediction_time = prediction_time,
        metric_generation_time = proc.time()[["elapsed"]] - metric_started,
        total_time = study06_run$runtime +
          study06_bundle$stats$elapsed_seconds + prediction_time,
        output_size_bytes_estimate = as.numeric(object.size(study06_run$fit)),
        memory_bytes_estimate = if (is.null(
          study06_run$fit$memory_estimate$estimated_total_bytes)) NA_real_ else
          as.numeric(study06_run$fit$memory_estimate$estimated_total_bytes),
        warnings = paste(study06_run$result$diagnostics$warnings, collapse = " | "),
        stringsAsFactors = FALSE))
  }, pattern = map(study06_run, study06_bundle, study06_run_checkpoint),
  iteration = "list"),
  targets::tar_target(study06_prediction_metrics,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "prediction"))),
  targets::tar_target(study06_parameter_estimates,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "parameter"))),
  targets::tar_target(study06_annotation_estimates,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "annotation"))),
  targets::tar_target(study06_sigma_alpha,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "sigma"))),
  targets::tar_target(study06_marker_metrics,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "marker"))),
  targets::tar_target(study06_probability_metrics,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "probability"))),
  targets::tar_target(study06_probability_inventory,
    do.call(rbind, lapply(study06_outputs_branch, `[[`,
      "probability_inventory"))),
  targets::tar_target(study06_enrichment_metrics,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "enrichment"))),
  targets::tar_target(study06_convergence_diagnostics,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "diagnostics"))),
  targets::tar_target(study06_computational,
    do.call(rbind, lapply(study06_outputs_branch, `[[`, "computational"))),
  targets::tar_target(study06_fit_status, study06_computational[
    c("scenario", "replicate", "method", "status", "error_message")]),
  targets::tar_target(study06_simulation_summary,
    do.call(rbind, lapply(study06_simulation, .study06_simulation_summary))),
  targets::tar_target(study06_simulation_marker_truth,
    do.call(rbind, lapply(study06_simulation, function(x)
      data.frame(scenario = x$scenario, replicate = x$replicate,
        x$marker_truth, check.names = FALSE, stringsAsFactors = FALSE)))),
  targets::tar_target(study06_parameter_recovery_summary, {
    z <- study06_parameter_estimates[study06_parameter_estimates$available, ]
    groups <- split(z, interaction(z$scenario, z$method, z$estimand, drop = TRUE))
    do.call(rbind, lapply(groups, function(x) data.frame(
      scenario = x$scenario[1L], method = x$method[1L],
      estimand = x$estimand[1L], replicate_count = 5L,
      successful_replicates = nrow(x), mean_bias = mean(x$bias),
      rmse = sqrt(mean(x$squared_error)), mae = mean(abs(x$bias)),
      observed_coverage_count = sum(x$interval_coverage),
      observed_coverage_proportion = mean(x$interval_coverage),
      coverage_note = "coarse five-replicate observed coverage; each miss changes 20 percentage points",
      mean_interval_width_95 = mean(x$upper_95 - x$lower_95),
      stringsAsFactors = FALSE)))
  }),
  targets::tar_target(study06_annotation_recovery_summary, {
    x <- study06_annotation_estimates
    x$available <- TRUE
    .study06_aggregate(transform(x, metric = paste(annotation_name,
      stick_name, sep = ":")), value = "bias",
      keys = c("scenario", "method", "metric"))
  }),
  targets::tar_target(study06_prediction_summary,
    .study06_aggregate(study06_prediction_metrics)),
  targets::tar_target(study06_probability_summary,
    .study06_aggregate(study06_probability_metrics,
      keys = c("scenario", "method", "probability_object", "metric"))),
  targets::tar_target(study06_enrichment_summary,
    .study06_aggregate(study06_enrichment_metrics,
      keys = c("scenario", "method", "probability_object", "metric"))),
  targets::tar_target(study06_comparison_metrics, {
    p <- transform(study06_prediction_metrics,
      metric = sub("^phenotype_prediction_", "prediction_", metric))
    h <- study06_parameter_estimates[
      study06_parameter_estimates$estimand == "heritability", ]
    h <- data.frame(scenario = h$scenario, replicate = h$replicate,
      method = h$method, metric = "heritability_absolute_error",
      value = abs(h$bias), available = h$available, reason = h$reason)
    m <- study06_marker_metrics[study06_marker_metrics$metric == "marker_effect_rmse", ]
    a <- study06_probability_metrics[
      study06_probability_metrics$metric == "posterior_nonnull_auprc", ]
    a <- a[c("scenario", "replicate", "method", "metric", "value",
      "available", "reason")]
    rbind(p, h, m, a)
  }),
  targets::tar_target(study06_paired_differences,
    .study06_paired(study06_comparison_metrics)),
  targets::tar_target(study06_paired_summary,
    .study06_paired_summary(study06_paired_differences)),
  targets::tar_target(study06_interactions,
    .study06_interactions(study06_paired_differences)),
  targets::tar_target(study06_interaction_summary, {
    x <- study06_interactions
    groups <- split(x, interaction(x$comparison_id, x$metric, drop = TRUE))
    do.call(rbind, lapply(groups, function(z) {
      y <- z$interaction_difference[z$complete_pair]
      data.frame(comparison_id = z$comparison_id[1L],
        metric = z$metric[1L], complete_paired_replicates = length(y),
        mean_difference = mean(y), sd_difference = stats::sd(y),
        median_difference = stats::median(y), minimum_difference = min(y),
        maximum_difference = max(y), stringsAsFactors = FALSE)
    }))
  }),
  targets::tar_target(study06_convergence_validation_summary,
    .study06_validation_summary(study06_convergence_diagnostics)),
  targets::tar_target(study06_seed_registry,
    .study06_effective_seed_registry(study06_run, study06_config)),
  targets::tar_target(study06_method_availability, {
    x <- expand.grid(
      method = study06_config$methods,
      output = c("marker_effects", "posterior_component_probabilities",
        "annotation_implied_prior_probabilities", "alpha_draws",
        "sigmaSqAlpha_draws", "global_component_proportion_draws"),
      stringsAsFactors = FALSE)
    x$available <- x$output %in% c("marker_effects",
      "posterior_component_probabilities") |
      (grepl("bayesrc", x$method) & x$output %in% c(
        "annotation_implied_prior_probabilities", "alpha_draws",
        "sigmaSqAlpha_draws"))
    x$reason <- ifelse(x$available, "",
      "not returned by current installed model")
    x
  }),
  targets::tar_target(study06_output_dir,
    Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
      file.path(study06_paths$local_dir, "benchmark_outputs"))),
  targets::tar_target(study06_output_files, {
    tables <- list(
      annotation_design_summary.csv = study06_annotation_summary,
      annotation_column_metadata.csv = data.frame(
        column = colnames(study06_annotation),
        role = c("explicit_intercept", "seeded_binary_signal",
          "seeded_standardized_signal", "seeded_standardized_null"),
        preprocessing = c("constant_one", "exact_ten_percent",
          "sample_mean_zero_sample_sd_one",
          "sample_mean_zero_sample_sd_one"),
        column_order = seq_len(ncol(study06_annotation))),
      sample_split.csv = data.frame(
        sample_id = study06_sample_ids,
        split = ifelse(study06_sample_ids %in% study06_split$train_ids,
          "training", "test"),
        split_seed = study06_config$split$seed),
      simulation_summary.csv = study06_simulation_summary,
      simulation_marker_truth.csv = study06_simulation_marker_truth,
      fit_status.csv = study06_fit_status,
      prediction_metrics.csv = study06_prediction_metrics,
      prediction_summary.csv = study06_prediction_summary,
      parameter_estimates.csv = study06_parameter_estimates,
      parameter_recovery_summary.csv = study06_parameter_recovery_summary,
      annotation_coefficient_estimates.csv = study06_annotation_estimates,
      annotation_recovery_summary.csv = study06_annotation_recovery_summary,
      annotation_variance_summary.csv = study06_sigma_alpha,
      marker_effect_metrics.csv = study06_marker_metrics,
      probability_recovery_metrics.csv = study06_probability_metrics,
      probability_recovery_summary.csv = study06_probability_summary,
      probability_output_inventory.csv = study06_probability_inventory,
      enrichment_metrics.csv = study06_enrichment_metrics,
      enrichment_summary.csv = study06_enrichment_summary,
      paired_replicate_differences.csv = study06_paired_differences,
      paired_comparison_summary.csv = study06_paired_summary,
      annotation_value_interactions.csv = study06_interactions,
      annotation_value_interaction_summary.csv = study06_interaction_summary,
      convergence_diagnostics.csv = study06_convergence_diagnostics,
      convergence_validation_summary.csv = study06_convergence_validation_summary,
      computational_summary.csv = study06_computational,
      seed_registry.csv = study06_seed_registry,
      method_output_availability.csv = study06_method_availability)
    paths <- character()
    for (name in names(tables)) {
      paths <- c(paths, .study06_write_csv(tables[[name]],
        file.path(study06_output_dir, name)))
    }
    alpha <- do.call(rbind, lapply(study06_config$scenarios, function(s)
      data.frame(scenario = s, annotation_name = rownames(study06_alpha[[s]]),
        study06_alpha[[s]], check.names = FALSE)))
    paths <- c(paths, .study06_write_csv(alpha,
      file.path(study06_output_dir, "true_alpha.csv")))
    annotation_design <- data.frame(marker_id = rownames(study06_annotation),
      study06_annotation, check.names = FALSE)
    paths <- c(paths, .study06_write_csv(annotation_design,
      file.path(study06_output_dir, "annotation_design.csv")))
    paths
  }, format = "file"),
  targets::tar_target(study06_manifest_file,
    .study06_write_json(list(study = "06_annotation_models",
      task = "single_trait_annotation_informed_models",
      benchmark_scope = "five_replicate_development",
      benchmark_status = if (nrow(study06_fit_status) == 40L &&
        all(study06_fit_status$status == "ok")) "complete" else "incomplete",
      scenarios = study06_config$scenarios, replicate_count = 5L,
      methods = study06_config$methods, expected_fit_count = 40L,
      successful_fit_count = sum(study06_fit_status$status == "ok"),
      failed_fit_count = sum(study06_fit_status$status != "ok"),
      chains_per_fit = 4L, expected_chain_count = 160L,
      mixture_var = study06_config$mixture_var,
      sample_count = length(study06_sample_ids),
      training_sample_count = length(study06_split$train_ids),
      test_sample_count = length(study06_split$test_ids),
      marker_count = length(study06_filtered_markers$marker_ids),
      annotation_columns = colnames(study06_annotation),
      annotation_dimensions = dim(study06_annotation),
      target_h2 = study06_config$simulation$h2,
      baseline_recommendations = .study06_baseline_recommendations(study06_config),
      annotation_recommendations = .study06_annotation_recommendations(study06_config),
      installed_sblr_version = as.character(packageVersion("sblr")),
      installed_sblr_commit = packageDescription("sblr")$RemoteSha,
      package_versions = as.list(vapply(c("sblr", "sblrbench", "targets",
        "posterior"), function(x) as.character(packageVersion(x)), "")),
      repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      qgdata = study06_config$example_data,
      seed_registry_file = "seed_registry.csv",
      source_status = "cached_checksum_validated_offline",
      completion_counts = list(successful_fits =
        sum(study06_fit_status$status == "ok"), failed_fits =
        sum(study06_fit_status$status != "ok")),
      failures = unname(study06_fit_status$error_message[
        study06_fit_status$status != "ok"]),
      provenance = list(installed_package_only = TRUE,
        sibling_source_access = "read_only_audit"),
      validation_status = "complete_grid_validated",
      coverage_interpretation = "coarse five-replicate observed coverage",
      limitations = c("five paired simulations are descriptive evidence",
        "not a definitive method ranking",
        "not universal convergence validation",
        "current installed sblr behaviour")),
      file.path(study06_output_dir, "benchmark_manifest.json")),
    format = "file")
))

phase <- Sys.getenv("SBLR_BENCH_STUDY06_PHASE", "benchmark")
if (identical(phase, "convergence")) {
  .study06_convergence_targets()
} else if (identical(phase, "benchmark")) {
  .study06_benchmark_targets()
} else {
  stop("SBLR_BENCH_STUDY06_PHASE must be convergence or benchmark.")
}
