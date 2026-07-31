source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
for (f in c("estimands.R", "simulation.R", "methods.R", "metrics.R", "pilot.R"))
  source(file.path("studies", "03_parameter_estimation", f), local = TRUE)
source(file.path("studies", "five_replicate_helpers.R"), local = TRUE)

list(
  targets::tar_target(parameter_config_file, file.path("studies", "03_parameter_estimation", "config.R"), format = "file"),
  targets::tar_target(parameter_config, source(parameter_config_file, local = TRUE)$value),
  targets::tar_target(parameter_registry, .study03_validate_registry(.study03_estimand_registry())),
  targets::tar_target(parameter_paths, .study03_paths()),
  targets::tar_target(parameter_example_files, .study01_example_files(parameter_paths$data_dir, parameter_config$example_data)),
  targets::tar_target(parameter_base_glist, .study01_load_glist(parameter_paths, parameter_example_files)),
  targets::tar_target(parameter_markers, .study01_run_qc(parameter_base_glist, parameter_config)),
  targets::tar_target(parameter_sample_ids, .study01_selected_ids(parameter_base_glist, parameter_config$sample_limit)),
  targets::tar_target(parameter_working_glist, .study01_set_rsids_ld(parameter_base_glist,
    parameter_config$chr, parameter_markers$marker_ids)),
  targets::tar_target(parameter_genotypes, .study01_extract_genotypes(parameter_working_glist,
    parameter_config$chr, parameter_sample_ids, parameter_markers$marker_ids)),
  targets::tar_target(parameter_sparse_ld_glist, .study01_make_sparse_ld(parameter_working_glist,
    parameter_markers, parameter_config, parameter_paths$output_dir)),
  targets::tar_target(parameter_specs, .study03_replicate_specs(parameter_config), iteration = "list"),
  targets::tar_target(parameter_spec, parameter_specs, pattern = map(parameter_specs), iteration = "list"),
  targets::tar_target(parameter_simulation_bundle, {
    simulation <- .study03_simulate(parameter_spec, parameter_genotypes, parameter_config)
    oracle <- sblrbench::check_oracle_genetic_values(simulation, tolerance = parameter_config$oracle_tolerance)
    if (!oracle$ok) stop("Study 03 simulation oracle failed.", call. = FALSE)
    truth <- .study03_truth(simulation, parameter_config)
    stats <- .study03_summary_stats(simulation, parameter_sparse_ld_glist, parameter_config)
    list(simulation = simulation, oracle = oracle, truth = truth, stats = stats)
  }, pattern = map(parameter_spec), iteration = "list"),
  targets::tar_target(parameter_methods, .study03_method_specs(parameter_config), iteration = "list"),
  targets::tar_target(parameter_method, parameter_methods, pattern = map(parameter_methods), iteration = "list"),
  targets::tar_target(parameter_seed_registry, .study03_seed_registry(parameter_specs, parameter_methods, parameter_config)),
  targets::tar_target(parameter_method_run, {
    fit <- .study03_fit(parameter_method, parameter_simulation_bundle$simulation,
      parameter_simulation_bundle$stats, parameter_sparse_ld_glist, parameter_config)
    summary <- data.frame()
    if (fit$status == "ok") tryCatch({
      draws <- if (identical(parameter_config$profile, "five_replicate_development"))
        .study03_extract_multichain_draws(fit$native_fit, fit$method$id,
          parameter_registry, length(parameter_simulation_bundle$simulation$data$marker_ids)) else
        .study03_extract_draws(fit$native_fit, fit$method$id, parameter_registry,
          length(parameter_simulation_bundle$simulation$data$marker_ids))
      summary <- .study03_summarise_draws(draws)
    }, error = function(e) {
      fit$status <<- "failed"; fit$reason <<- conditionMessage(e)
    })
    summary <- .study03_complete_summary(summary, parameter_registry, fit$method$id)
    recovery <- .study03_recovery(summary, parameter_simulation_bundle$truth,
      fit$method, parameter_registry, parameter_config$relative_error_tolerance)
    computational <- data.frame(architecture = parameter_simulation_bundle$simulation$scenario$architecture,
      replicate = parameter_simulation_bundle$simulation$scenario$replicate, method = fit$method$id,
      runtime = if (fit$status == "ok") fit$result$computation$elapsed_seconds else NA_real_,
      status = fit$status, error_class = if (fit$status == "ok") "" else "fit_error",
      error_message = fit$reason, nit = fit$controls$nit, nburn = fit$controls$nburn,
      nthin = fit$controls$nthin, nchains = fit$controls$nchains,
      ncores = fit$controls$ncores,
      chain_count = if (fit$status == "ok") max(recovery$chain_count, na.rm = TRUE) else 0L,
      draws_per_chain = if (fit$status == "ok") max(recovery$draws_per_chain, na.rm = TRUE) else 0L,
      posterior_draw_count = if (fit$status == "ok") max(recovery$n_posterior_draws, na.rm = TRUE) else 0L,
      fit_seed = fit$seed, stringsAsFactors = FALSE)
    list(fit = fit, recovery = recovery, computational = computational)
  }, pattern = cross(parameter_simulation_bundle, parameter_method), iteration = "list"),
  targets::tar_target(parameter_recovery_branch, parameter_method_run$recovery,
    pattern = map(parameter_method_run), iteration = "list"),
  targets::tar_target(parameter_computational_branch, parameter_method_run$computational,
    pattern = map(parameter_method_run), iteration = "list"),
  targets::tar_target(parameter_truth_branch, parameter_simulation_bundle$truth,
    pattern = map(parameter_simulation_bundle), iteration = "list"),
  targets::tar_target(parameter_estimates, do.call(rbind, parameter_recovery_branch)),
  targets::tar_target(parameter_computational_summary,
    do.call(rbind, parameter_computational_branch)),
  targets::tar_target(parameter_recovery_summary, .study03_aggregate(parameter_estimates)),
  targets::tar_target(parameter_paired_differences, .study03_paired(parameter_estimates)),
  targets::tar_target(parameter_paired_summary, .study03_paired_summary(parameter_paired_differences)),
  targets::tar_target(parameter_simulation_truth, do.call(rbind, parameter_truth_branch)),
  targets::tar_target(parameter_replicate_status, parameter_computational_summary[, c("architecture", "replicate", "method", "status", "error_class", "error_message")]),
  targets::tar_target(parameter_output_dir,
    Sys.getenv("SBLR_BENCH_OUTPUT_DIR", file.path("results", "local", "03_parameter_estimation"))),
  targets::tar_target(parameter_registry_file, .study03_write_csv(parameter_registry, file.path(parameter_output_dir, "estimand_registry.csv")), format = "file"),
  targets::tar_target(parameter_truth_file, .study03_write_csv(parameter_simulation_truth, file.path(parameter_output_dir, "simulation_truth.csv")), format = "file"),
  targets::tar_target(parameter_estimates_file, .study03_write_csv(parameter_estimates, file.path(parameter_output_dir, "parameter_estimates.csv")), format = "file"),
  targets::tar_target(parameter_summary_file, .study03_write_csv(parameter_recovery_summary, file.path(parameter_output_dir, "parameter_recovery_summary.csv")), format = "file"),
  targets::tar_target(parameter_paired_file, .study03_write_csv(parameter_paired_differences, file.path(parameter_output_dir, "paired_parameter_differences.csv")), format = "file"),
  targets::tar_target(parameter_paired_summary_file, .study03_write_csv(parameter_paired_summary, file.path(parameter_output_dir, "paired_comparison_summary.csv")), format = "file"),
  targets::tar_target(parameter_computational_file, .study03_write_csv(parameter_computational_summary, file.path(parameter_output_dir, "computational_summary.csv")), format = "file"),
  targets::tar_target(parameter_status_file, .study03_write_csv(parameter_replicate_status, file.path(parameter_output_dir, "replicate_status.csv")), format = "file"),
  targets::tar_target(parameter_seed_file, .study03_write_csv(parameter_seed_registry, file.path(parameter_output_dir, "seed_registry.csv")), format = "file"),
  targets::tar_target(parameter_manifest_file, {
    path <- file.path(parameter_output_dir, "benchmark_manifest.json")
    nrep <- parameter_config$profiles[[parameter_config$profile]]$replicate_count
    expected_fits <- 2L * nrep * length(parameter_config$methods)
    sblr_provenance <- .five_replicate_sblr_provenance()
    jsonlite::write_json(list(study_id = parameter_config$study, task = parameter_config$task,
      benchmark_scope = if (nrep == 5L) "five_replicate_development" else "one_replicate_development",
      benchmark_status = if (nrow(parameter_replicate_status) == expected_fits &&
        all(parameter_replicate_status$status == "ok")) "complete" else "incomplete",
      development_settings = TRUE, created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      sblr_version = sblr_provenance$version, sblr_commit = sblr_provenance$commit,
      sblr_source_status = sblr_provenance$source_status,
      qgdata_commit = parameter_config$example_data$commit,
      active_methods = parameter_config$methods,
      architectures = names(parameter_config$simulation$architectures), replicate_count = nrep,
      expected_fit_count = expected_fits, successful_fit_count = sum(parameter_replicate_status$status == "ok"),
      failed_fit_count = sum(parameter_replicate_status$status != "ok"),
      analysis_sample_count = length(parameter_sample_ids), canonical_marker_count = length(parameter_markers$marker_ids),
      causal_marker_count = parameter_config$simulation$n_causal, target_h2 = parameter_config$simulation$h2,
      mcmc = if (nrep == 5L) .five_replicate_recommendations() else parameter_config$profiles$development,
      seeds = parameter_seed_registry, thread_settings = .five_replicate_thread_settings(),
      coverage_interpretation = "coarse observed coverage; each miss changes five-replicate coverage by 20 percentage points",
      source_status = "working-tree five-replicate source; source inventory retained in capsule",
      completion_counts = list(expected = expected_fits,
        successful = sum(parameter_replicate_status$status == "ok"),
        failed = sum(parameter_replicate_status$status != "ok")),
      failures = parameter_replicate_status[parameter_replicate_status$status != "ok", ],
      provenance = parameter_config$example_data,
      validation_status = if (nrow(parameter_replicate_status) == expected_fits &&
        all(parameter_replicate_status$status == "ok")) "complete_grid_validated" else "incomplete_not_promotable",
      estimands = parameter_registry,
      primary_estimands = parameter_registry$estimand_id[parameter_registry$primary],
      truth_definitions = unique(parameter_simulation_truth[c("estimand_id", "truth_type", "truth_definition")]),
      data_provenance = parameter_config$example_data), path, pretty = TRUE, auto_unbox = TRUE)
    path
  }, format = "file")
)
