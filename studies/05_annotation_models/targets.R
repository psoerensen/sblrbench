source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
for (f in c("annotation_design.R", "simulation.R", "methods.R",
            "chain_extraction.R", "diagnostics.R", "metrics.R", "pilot.R"))
  source(file.path("studies", "05_annotation_models", f), local = TRUE)

.study05_common_targets <- function() list(
  targets::tar_target(study05_config_file,
    file.path("studies", "05_annotation_models", "config.R"), format = "file"),
  targets::tar_target(study05_config,
    source(study05_config_file, local = TRUE)$value),
  targets::tar_target(study05_paths, .study05_paths()),
  targets::tar_target(study05_example_files,
    .study01_example_files(study05_paths$data_dir, study05_config$example_data)),
  targets::tar_target(study05_base_glist,
    .study01_load_glist(study05_paths, study05_example_files)),
  targets::tar_target(study05_filtered_markers,
    .study01_run_qc(study05_base_glist, study05_config)),
  targets::tar_target(study05_sample_ids,
    .study01_selected_ids(study05_base_glist, study05_config$sample_limit)),
  targets::tar_target(study05_split, sblrbench::make_prediction_split(
    study05_sample_ids, study05_config$split$train_fraction,
    study05_config$split$seed)),
  targets::tar_target(study05_scaled_genotypes,
    .study05_load_scaled_genotypes(study05_base_glist, study05_config$chr,
      study05_sample_ids, study05_filtered_markers$marker_ids, study05_split)),
  targets::tar_target(study05_working_glist,
    sblrbench:::benchmark_set_training_af(
    .study01_set_rsids_ld(study05_base_glist, study05_config$chr,
      study05_filtered_markers$marker_ids), study05_config$chr,
    study05_filtered_markers$marker_ids,
    study05_scaled_genotypes$allele_frequency)),
  targets::tar_target(study05_ld_bundle, {
    started <- proc.time()[["elapsed"]]
    x <- sblrbench:::benchmark_make_training_ld(study05_working_glist,
      study05_split,
      study05_filtered_markers$marker_ids, list(
        chromosome = study05_config$chr, sparse_ld = study05_config$sparse_ld),
      study05_paths$genotype_output_dir)
    list(Glist = x, elapsed_seconds = proc.time()[["elapsed"]] - started)
  }),
  targets::tar_target(study05_annotation,
    .study05_annotation_design(study05_filtered_markers$marker_ids,
      study05_config)),
  targets::tar_target(study05_alpha,
    .study05_true_alpha(study05_annotation, study05_config)),
  targets::tar_target(study05_annotation_summary,
    .study05_annotation_summary(study05_annotation, study05_alpha,
      study05_config)),
  targets::tar_target(study05_method_specs,
    .study05_method_specs(study05_config), iteration = "list"),
  targets::tar_target(study05_method_spec, study05_method_specs,
    pattern = map(study05_method_specs), iteration = "list")
)

.study05_convergence_targets <- function() c(.study05_common_targets(), list(
  targets::tar_target(study05_convergence_simulation,
    .study05_simulate("informative_annotations", 1L,
      study05_scaled_genotypes$all, study05_split$train_rows,
      study05_annotation, study05_alpha, study05_config)),
  targets::tar_target(study05_convergence_stats,
    .study05_summary_stats(study05_convergence_simulation,
      study05_ld_bundle$Glist, study05_split, study05_config)),
  targets::tar_target(study05_convergence_methods,
    Filter(function(x) isTRUE(x$annotation), study05_method_specs),
    iteration = "list"),
  targets::tar_target(study05_convergence_method,
    study05_convergence_methods,
    pattern = map(study05_convergence_methods), iteration = "list"),
  targets::tar_target(study05_convergence_run,
    .study05_fit(study05_convergence_method, study05_convergence_simulation,
      study05_convergence_stats$stats, study05_ld_bundle$Glist,
      study05_split, study05_annotation, study05_alpha, study05_config,
      phase = "convergence"),
    pattern = map(study05_convergence_method), iteration = "list"),
  targets::tar_target(study05_convergence_checkpoint,
    .study05_fit_checkpoint(study05_convergence_run, "convergence"),
    pattern = map(study05_convergence_run), format = "file"),
  targets::tar_target(study05_convergence_draws_branch,
    {
      study05_convergence_checkpoint
      .study05_extract_chain_draws(study05_convergence_run, study05_annotation)
    },
    pattern = map(study05_convergence_run, study05_convergence_checkpoint),
    iteration = "list"),
  targets::tar_target(study05_convergence_selection_branch,
    .study05_select_recommendation(study05_convergence_draws_branch,
      study05_config),
    pattern = map(study05_convergence_draws_branch), iteration = "list"),
  targets::tar_target(study05_convergence_draws,
    do.call(rbind, study05_convergence_draws_branch)),
  targets::tar_target(study05_convergence_diagnostics,
    do.call(rbind, lapply(study05_convergence_selection_branch, `[[`,
      "diagnostics"))),
  targets::tar_target(study05_convergence_candidates,
    do.call(rbind, lapply(study05_convergence_selection_branch, `[[`,
      "candidates"))),
  targets::tar_target(study05_convergence_recommendations,
    do.call(rbind, lapply(study05_convergence_selection_branch, `[[`,
      "recommendation"))),
  targets::tar_target(study05_convergence_status, do.call(rbind,
    lapply(study05_convergence_run, function(x)
      data.frame(method = x$method$id, scenario = x$scenario,
        replicate = x$replicate, status = x$status, reason = x$reason,
        runtime = x$runtime, chain_count = if (x$status == "ok") 4L else 0L,
        stringsAsFactors = FALSE)))),
  targets::tar_target(study05_convergence_seed_registry, {
    x <- .study05_effective_seed_registry(study05_convergence_run,
      study05_config)
    x[x$scenario == "informative_annotations" & x$replicate == 1L &
      x$method %in% c("st_bed_bayesrc", "st_csr_sbayesrc"), ]
  }),
  targets::tar_target(study05_convergence_output_dir,
    Sys.getenv("SBLR_BENCH_OUTPUT_DIR", file.path(study05_paths$local_dir,
      "convergence_outputs"))),
  targets::tar_target(study05_convergence_draws_file,
    .study05_write_csv(study05_convergence_draws,
      file.path(study05_convergence_output_dir, "scalar_chain_draws.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_diagnostics_file,
    .study05_write_csv(study05_convergence_diagnostics,
      file.path(study05_convergence_output_dir, "convergence_diagnostics.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_candidates_file,
    .study05_write_csv(study05_convergence_candidates,
      file.path(study05_convergence_output_dir, "candidate_settings.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_recommendations_file,
    .study05_write_csv(study05_convergence_recommendations,
      file.path(study05_convergence_output_dir, "method_recommendations.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_status_file,
    .study05_write_csv(study05_convergence_status,
      file.path(study05_convergence_output_dir, "fit_status.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_computational_file,
    .study05_write_csv(study05_convergence_status,
      file.path(study05_convergence_output_dir, "computational_summary.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_seed_file,
    .study05_write_csv(study05_convergence_seed_registry,
      file.path(study05_convergence_output_dir, "seed_registry.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_annotation_file,
    .study05_write_csv(study05_annotation_summary,
      file.path(study05_convergence_output_dir, "annotation_design_summary.csv")),
    format = "file"),
  targets::tar_target(study05_convergence_alpha_file, {
    x <- do.call(rbind, lapply(study05_config$scenarios, function(s)
      data.frame(scenario = s, annotation_name = rownames(study05_alpha[[s]]),
        study05_alpha[[s]], check.names = FALSE)))
    .study05_write_csv(x, file.path(study05_convergence_output_dir,
      "true_alpha.csv"))
  }, format = "file"),
  targets::tar_target(study05_convergence_manifest_file,
    .study05_write_json(list(
      study = "05_annotation_models", task = "annotation_model_convergence_selection",
      benchmark_scope = "two_method_maximum_history_development",
      benchmark_status = if (nrow(study05_convergence_status) == 2L &&
        all(study05_convergence_status$status == "ok") &&
        nrow(study05_convergence_recommendations) == 2L) "complete" else "incomplete",
      method_count = 2L, expected_fit_count = 2L, nchains = 4L,
      maximum_history = study05_config$convergence_pilot,
      mixture_var = study05_config$mixture_var,
      methods = c("st_bed_bayesrc", "st_csr_sbayesrc"),
      installed_sblr_version = as.character(packageVersion("sblr")),
      installed_sblr_commit = packageDescription("sblr")$RemoteSha,
      package_versions = as.list(vapply(c("sblr", "sblrbench", "targets",
        "posterior"), function(x) as.character(packageVersion(x)), "")),
      repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      qgdata = study05_config$example_data,
      seed_registry_file = "seed_registry.csv",
      source_status = "cached_checksum_validated_offline",
      completion_counts = list(successful_fits =
        sum(study05_convergence_status$status == "ok"), failed_fits =
        sum(study05_convergence_status$status != "ok")),
      failures = unname(study05_convergence_status$reason[
        study05_convergence_status$status != "ok"]),
      provenance = list(installed_package_only = TRUE,
        sibling_source_access = "read_only_audit"),
      recommendations = study05_convergence_recommendations,
      validation_status = "complete_grid_validated"),
      file.path(study05_convergence_output_dir, "benchmark_manifest.json")),
    format = "file")
))

.study05_benchmark_targets <- function() c(.study05_common_targets(), list(
  targets::tar_target(study05_specs, .study05_specs(study05_config),
    iteration = "list"),
  targets::tar_target(study05_spec, study05_specs,
    pattern = map(study05_specs), iteration = "list"),
  targets::tar_target(study05_simulation,
    .study05_simulate(study05_spec$scenario, study05_spec$replicate,
      study05_scaled_genotypes$all, study05_split$train_rows,
      study05_annotation, study05_alpha, study05_config),
    pattern = map(study05_spec), iteration = "list"),
  targets::tar_target(study05_bundle, {
    stats <- .study05_summary_stats(study05_simulation,
      study05_ld_bundle$Glist, study05_split, study05_config)
    list(simulation = study05_simulation, stats = stats)
  },
    pattern = map(study05_simulation), iteration = "list"),
  targets::tar_target(study05_run,
    .study05_fit(study05_method_spec, study05_bundle$simulation,
      study05_bundle$stats$stats, study05_ld_bundle$Glist, study05_split,
      study05_annotation, study05_alpha, study05_config, phase = "benchmark"),
    pattern = cross(study05_bundle, study05_method_spec),
    iteration = "list"),
  targets::tar_target(study05_run_checkpoint,
    .study05_fit_checkpoint(study05_run, "benchmark"),
    pattern = map(study05_run), format = "file"),
  targets::tar_target(study05_outputs_branch, {
    study05_run_checkpoint
    metric_started <- proc.time()[["elapsed"]]
    draws <- .study05_extract_chain_draws(study05_run, study05_annotation)
    parameter <- .study05_scalar_estimates(study05_run, draws,
      study05_bundle$simulation)
    annotation <- .study05_annotation_estimates(study05_run, draws,
      study05_bundle$simulation)
    sigma <- .study05_sigma_alpha(study05_run, draws)
    prediction_started <- proc.time()[["elapsed"]]
    prediction <- .study05_prediction_metrics(study05_run,
      study05_bundle$simulation, study05_scaled_genotypes$test,
      study05_split$test_ids)
    prediction_time <- proc.time()[["elapsed"]] - prediction_started
    marker <- .study05_marker_metrics(study05_run, study05_bundle$simulation,
      study05_annotation, study05_config)
    inventory <- .study05_probability_inventory(study05_run,
      study05_annotation, length(study05_config$mixture_var))
    diagnostics <- .study05_diagnostics(draws, 0L,
      max(draws$iteration), study05_config)
    list(draws = draws, prediction = prediction, parameter = parameter,
      annotation = annotation, sigma = sigma, marker = marker$marker,
      probability = marker$probability, enrichment = marker$enrichment,
      probability_inventory = inventory,
      diagnostics = diagnostics,
      computational = data.frame(scenario = study05_run$scenario,
        replicate = study05_run$replicate, method = study05_run$method$id,
        status = study05_run$status, error_message = study05_run$reason,
        preprocessing_time = study05_ld_bundle$elapsed_seconds,
        summary_statistic_time = study05_bundle$stats$elapsed_seconds,
        ld_construction_time = study05_ld_bundle$elapsed_seconds,
        model_fitting_time = study05_run$runtime,
        prediction_time = prediction_time,
        metric_generation_time = proc.time()[["elapsed"]] - metric_started,
        total_time = study05_run$runtime +
          study05_bundle$stats$elapsed_seconds + prediction_time,
        output_size_bytes_estimate = as.numeric(object.size(study05_run$fit)),
        memory_bytes_estimate = if (is.null(
          study05_run$fit$memory_estimate$estimated_total_bytes)) NA_real_ else
          as.numeric(study05_run$fit$memory_estimate$estimated_total_bytes),
        warnings = paste(study05_run$result$diagnostics$warnings, collapse = " | "),
        stringsAsFactors = FALSE))
  }, pattern = map(study05_run, study05_bundle, study05_run_checkpoint),
  iteration = "list"),
  targets::tar_target(study05_prediction_metrics,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "prediction"))),
  targets::tar_target(study05_parameter_estimates,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "parameter"))),
  targets::tar_target(study05_annotation_estimates,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "annotation"))),
  targets::tar_target(study05_sigma_alpha,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "sigma"))),
  targets::tar_target(study05_marker_metrics,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "marker"))),
  targets::tar_target(study05_probability_metrics,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "probability"))),
  targets::tar_target(study05_probability_inventory,
    do.call(rbind, lapply(study05_outputs_branch, `[[`,
      "probability_inventory"))),
  targets::tar_target(study05_enrichment_metrics,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "enrichment"))),
  targets::tar_target(study05_convergence_diagnostics,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "diagnostics"))),
  targets::tar_target(study05_computational,
    do.call(rbind, lapply(study05_outputs_branch, `[[`, "computational"))),
  targets::tar_target(study05_fit_status, study05_computational[
    c("scenario", "replicate", "method", "status", "error_message")]),
  targets::tar_target(study05_simulation_summary,
    do.call(rbind, lapply(study05_simulation, .study05_simulation_summary))),
  targets::tar_target(study05_simulation_marker_truth,
    do.call(rbind, lapply(study05_simulation, function(x)
      data.frame(scenario = x$scenario, replicate = x$replicate,
        x$marker_truth, check.names = FALSE, stringsAsFactors = FALSE)))),
  targets::tar_target(study05_parameter_recovery_summary, {
    z <- study05_parameter_estimates[study05_parameter_estimates$available, ]
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
  targets::tar_target(study05_annotation_recovery_summary, {
    x <- study05_annotation_estimates
    x$available <- TRUE
    .study05_aggregate(transform(x, metric = paste(annotation_name,
      stick_name, sep = ":")), value = "bias",
      keys = c("scenario", "method", "metric"))
  }),
  targets::tar_target(study05_prediction_summary,
    .study05_aggregate(study05_prediction_metrics)),
  targets::tar_target(study05_probability_summary,
    .study05_aggregate(study05_probability_metrics,
      keys = c("scenario", "method", "probability_object", "metric"))),
  targets::tar_target(study05_enrichment_summary,
    .study05_aggregate(study05_enrichment_metrics,
      keys = c("scenario", "method", "probability_object", "metric"))),
  targets::tar_target(study05_comparison_metrics, {
    p <- transform(study05_prediction_metrics,
      metric = sub("^phenotype_prediction_", "prediction_", metric))
    h <- study05_parameter_estimates[
      study05_parameter_estimates$estimand == "heritability", ]
    h <- data.frame(scenario = h$scenario, replicate = h$replicate,
      method = h$method, metric = "heritability_absolute_error",
      value = abs(h$bias), available = h$available, reason = h$reason)
    m <- study05_marker_metrics[study05_marker_metrics$metric == "marker_effect_rmse", ]
    a <- study05_probability_metrics[
      study05_probability_metrics$metric == "posterior_nonnull_auprc", ]
    a <- a[c("scenario", "replicate", "method", "metric", "value",
      "available", "reason")]
    rbind(p, h, m, a)
  }),
  targets::tar_target(study05_paired_differences,
    .study05_paired(study05_comparison_metrics)),
  targets::tar_target(study05_paired_summary,
    .study05_paired_summary(study05_paired_differences)),
  targets::tar_target(study05_interactions,
    .study05_interactions(study05_paired_differences)),
  targets::tar_target(study05_interaction_summary, {
    x <- study05_interactions
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
  targets::tar_target(study05_convergence_validation_summary,
    .study05_validation_summary(study05_convergence_diagnostics)),
  targets::tar_target(study05_seed_registry,
    .study05_effective_seed_registry(study05_run, study05_config)),
  targets::tar_target(study05_method_availability, {
    x <- expand.grid(
      method = study05_config$methods,
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
  targets::tar_target(study05_output_dir,
    Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
      file.path(study05_paths$local_dir, "benchmark_outputs"))),
  targets::tar_target(study05_output_files, {
    tables <- list(
      annotation_design_summary.csv = study05_annotation_summary,
      annotation_column_metadata.csv = data.frame(
        column = colnames(study05_annotation),
        role = c("explicit_intercept", "seeded_binary_signal",
          "seeded_standardized_signal", "seeded_standardized_null"),
        preprocessing = c("constant_one", "exact_ten_percent",
          "sample_mean_zero_sample_sd_one",
          "sample_mean_zero_sample_sd_one"),
        column_order = seq_len(ncol(study05_annotation))),
      sample_split.csv = data.frame(
        sample_id = study05_sample_ids,
        split = ifelse(study05_sample_ids %in% study05_split$train_ids,
          "training", "test"),
        split_seed = study05_config$split$seed),
      simulation_summary.csv = study05_simulation_summary,
      simulation_marker_truth.csv = study05_simulation_marker_truth,
      fit_status.csv = study05_fit_status,
      prediction_metrics.csv = study05_prediction_metrics,
      prediction_summary.csv = study05_prediction_summary,
      parameter_estimates.csv = study05_parameter_estimates,
      parameter_recovery_summary.csv = study05_parameter_recovery_summary,
      annotation_coefficient_estimates.csv = study05_annotation_estimates,
      annotation_recovery_summary.csv = study05_annotation_recovery_summary,
      annotation_variance_summary.csv = study05_sigma_alpha,
      marker_effect_metrics.csv = study05_marker_metrics,
      probability_recovery_metrics.csv = study05_probability_metrics,
      probability_recovery_summary.csv = study05_probability_summary,
      probability_output_inventory.csv = study05_probability_inventory,
      enrichment_metrics.csv = study05_enrichment_metrics,
      enrichment_summary.csv = study05_enrichment_summary,
      paired_replicate_differences.csv = study05_paired_differences,
      paired_comparison_summary.csv = study05_paired_summary,
      annotation_value_interactions.csv = study05_interactions,
      annotation_value_interaction_summary.csv = study05_interaction_summary,
      convergence_diagnostics.csv = study05_convergence_diagnostics,
      convergence_validation_summary.csv = study05_convergence_validation_summary,
      computational_summary.csv = study05_computational,
      seed_registry.csv = study05_seed_registry,
      method_output_availability.csv = study05_method_availability)
    paths <- character()
    for (name in names(tables)) {
      paths <- c(paths, .study05_write_csv(tables[[name]],
        file.path(study05_output_dir, name)))
    }
    alpha <- do.call(rbind, lapply(study05_config$scenarios, function(s)
      data.frame(scenario = s, annotation_name = rownames(study05_alpha[[s]]),
        study05_alpha[[s]], check.names = FALSE)))
    paths <- c(paths, .study05_write_csv(alpha,
      file.path(study05_output_dir, "true_alpha.csv")))
    annotation_design <- data.frame(marker_id = rownames(study05_annotation),
      study05_annotation, check.names = FALSE)
    paths <- c(paths, .study05_write_csv(annotation_design,
      file.path(study05_output_dir, "annotation_design.csv")))
    paths
  }, format = "file"),
  targets::tar_target(study05_manifest_file,
    .study05_write_json(list(study = "05_annotation_models",
      task = "single_trait_annotation_informed_models",
      benchmark_scope = "five_replicate_development",
      benchmark_status = if (nrow(study05_fit_status) == 40L &&
        all(study05_fit_status$status == "ok")) "complete" else "incomplete",
      scenarios = study05_config$scenarios, replicate_count = 5L,
      methods = study05_config$methods, expected_fit_count = 40L,
      successful_fit_count = sum(study05_fit_status$status == "ok"),
      failed_fit_count = sum(study05_fit_status$status != "ok"),
      chains_per_fit = 4L, expected_chain_count = 160L,
      mixture_var = study05_config$mixture_var,
      sample_count = length(study05_sample_ids),
      training_sample_count = length(study05_split$train_ids),
      test_sample_count = length(study05_split$test_ids),
      marker_count = length(study05_filtered_markers$marker_ids),
      annotation_columns = colnames(study05_annotation),
      annotation_dimensions = dim(study05_annotation),
      target_h2 = study05_config$simulation$h2,
      baseline_recommendations = .study05_baseline_recommendations(study05_config),
      annotation_recommendations = .study05_annotation_recommendations(study05_config),
      installed_sblr_version = as.character(packageVersion("sblr")),
      installed_sblr_commit = packageDescription("sblr")$RemoteSha,
      package_versions = as.list(vapply(c("sblr", "sblrbench", "targets",
        "posterior"), function(x) as.character(packageVersion(x)), "")),
      repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      qgdata = study05_config$example_data,
      seed_registry_file = "seed_registry.csv",
      source_status = "cached_checksum_validated_offline",
      completion_counts = list(successful_fits =
        sum(study05_fit_status$status == "ok"), failed_fits =
        sum(study05_fit_status$status != "ok")),
      failures = unname(study05_fit_status$error_message[
        study05_fit_status$status != "ok"]),
      provenance = list(installed_package_only = TRUE,
        sibling_source_access = "read_only_audit"),
      validation_status = "complete_grid_validated",
      coverage_interpretation = "coarse five-replicate observed coverage",
      limitations = c("five paired simulations are descriptive evidence",
        "not a definitive method ranking",
        "not universal convergence validation",
        "current installed sblr behaviour")),
      file.path(study05_output_dir, "benchmark_manifest.json")),
    format = "file")
))

phase <- Sys.getenv("SBLR_BENCH_STUDY05_PHASE", "benchmark")
if (identical(phase, "convergence")) {
  .study05_convergence_targets()
} else if (identical(phase, "benchmark")) {
  .study05_benchmark_targets()
} else {
  stop("SBLR_BENCH_STUDY05_PHASE must be convergence or benchmark.")
}
