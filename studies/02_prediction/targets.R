source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
source(file.path("studies", "02_prediction", "pilot.R"), local = TRUE)

list(
  targets::tar_target(prediction_config_file, file.path("studies", "02_prediction", "config.R"), format = "file"),
  targets::tar_target(prediction_config, source(prediction_config_file, local = TRUE)$value),
  targets::tar_target(prediction_paths, .study02_paths()),
  targets::tar_target(prediction_example_files, .study01_example_files(prediction_paths$data_dir, prediction_config$example_data)),
  targets::tar_target(prediction_base_glist, .study01_load_glist(prediction_paths, prediction_example_files)),
  targets::tar_target(prediction_filtered_markers, .study01_run_qc(prediction_base_glist, prediction_config)),
  targets::tar_target(prediction_sample_ids, .study01_selected_ids(prediction_base_glist, prediction_config$sample_limit)),
  targets::tar_target(prediction_split, sblrbench::make_prediction_split(prediction_sample_ids,
    prediction_config$split$train_fraction, prediction_config$split$seed)),
  targets::tar_target(prediction_raw_genotypes, .study02_extract_raw(prediction_base_glist,
    prediction_config$chr, prediction_sample_ids, prediction_filtered_markers$marker_ids)),
  targets::tar_target(prediction_scaled_genotypes, sblrbench::training_scaled_genotypes(
    prediction_raw_genotypes, prediction_split$train_rows)),
  targets::tar_target(prediction_working_glist, .study02_set_training_af(
    .study01_set_rsids_ld(prediction_base_glist, prediction_config$chr,
      prediction_filtered_markers$marker_ids), prediction_config$chr,
    prediction_filtered_markers$marker_ids, prediction_scaled_genotypes$allele_frequency)),
  targets::tar_target(prediction_sparse_ld_glist, .study02_make_ld(prediction_working_glist,
    prediction_split, prediction_filtered_markers$marker_ids,
    prediction_config, prediction_paths$output_dir)),
  targets::tar_target(prediction_replicate_override, Sys.getenv("SBLR_BENCH_REPLICATES", ""), cue = targets::tar_cue(mode = "always")),
  targets::tar_target(prediction_architecture_override, Sys.getenv("SBLR_BENCH_ARCHITECTURE", ""), cue = targets::tar_cue(mode = "always")),
  targets::tar_target(prediction_replicate_specs, .study02_replicate_specs(prediction_config,
    prediction_replicate_override, prediction_architecture_override), iteration = "list"),
  targets::tar_target(prediction_replicate_spec, prediction_replicate_specs,
    pattern = map(prediction_replicate_specs), iteration = "list"),
  targets::tar_target(prediction_simulation_bundle, {
    simulation <- .study02_simulate(prediction_replicate_spec,
      prediction_scaled_genotypes$all, prediction_config)
    simulation$data$train_ids <- prediction_split$train_ids
    simulation$data$test_ids <- prediction_split$test_ids
    sblrbench::validate_sblrbench_simulation(simulation)
    oracle <- sblrbench::check_oracle_genetic_values(simulation,
      tolerance = prediction_config$oracle_tolerance)
    test_simulation <- sblrbench::subset_sblrbench_simulation_samples(
      simulation, prediction_split$test_ids)
    stats <- .study02_summary_stats(simulation, prediction_sparse_ld_glist,
      prediction_split, prediction_config)
    list(simulation = simulation, test_simulation = test_simulation,
      oracle = oracle, stats = stats, spec = prediction_replicate_spec)
  }, pattern = map(prediction_replicate_spec), iteration = "list"),
  targets::tar_target(prediction_method_specs,
    .study02_method_specs(prediction_config), iteration = "list"),
  targets::tar_target(prediction_method_spec, prediction_method_specs,
    pattern = map(prediction_method_specs), iteration = "list"),
  targets::tar_target(prediction_method_run, {
    fit <- .study02_fit(prediction_method_spec,
      prediction_simulation_bundle$simulation, prediction_simulation_bundle$stats,
      prediction_sparse_ld_glist, prediction_split, prediction_config)
    fit <- .study02_predict(fit, prediction_simulation_bundle$simulation,
      prediction_simulation_bundle$test_simulation, prediction_scaled_genotypes$test)
    comp <- data.frame(
      architecture = prediction_simulation_bundle$simulation$scenario$architecture,
      replicate = prediction_simulation_bundle$simulation$scenario$replicate,
      method = fit$method$id,
      runtime = if (fit$status == "ok") fit$result$computation$elapsed_seconds else NA_real_,
      status = fit$status,
      warnings = if (fit$status == "ok") paste(fit$result$diagnostics$warnings, collapse = " | ") else "",
      reason = fit$reason,
      simulation_seed = prediction_simulation_bundle$simulation$provenance$seed,
      split_seed = prediction_split$split_seed, mcmc_seed = fit$mcmc_seed,
      sblr_version = as.character(packageVersion("sblr")),
      qgg_version = as.character(packageVersion("qgg")),
      sblrbench_commit = sblrbench::sblrbench_git_commit("."),
      method_controls = jsonlite::toJSON(fit$controls, auto_unbox = TRUE),
      stringsAsFactors = FALSE)
    list(fit = fit,
      metrics = .study02_metrics(fit, prediction_simulation_bundle$test_simulation),
      computational = comp)
  }, pattern = cross(prediction_simulation_bundle, prediction_method_spec), iteration = "list"),
  targets::tar_target(prediction_computational_summary,
    do.call(rbind, lapply(prediction_method_run, `[[`, "computational"))),
  targets::tar_target(prediction_metrics, {
    z <- do.call(rbind, lapply(prediction_method_run, `[[`, "metrics"))
    names(z)[names(z) == "method_id"] <- "method"
    names(z)[names(z) == "scenario"] <- "architecture"
    z
  }),
  targets::tar_target(prediction_paired_method_differences,
    .study02_paired(prediction_metrics)),
  targets::tar_target(prediction_replicate_status, {
    expected <- expand.grid(architecture = names(prediction_config$simulation$architectures),
      replicate = seq_len(length(prediction_replicate_specs) /
        length(prediction_config$simulation$architectures)),
      method = prediction_config$methods, stringsAsFactors = FALSE)
    z <- merge(expected, prediction_computational_summary[, c("architecture", "replicate",
      "method", "status", "reason")], by = c("architecture", "replicate", "method"),
      all.x = TRUE, sort = TRUE)
    z$status[is.na(z$status)] <- "missing"
    z$reason[is.na(z$reason)] <- "required active method branch is absent"
    z$expected_fit_count <- nrow(expected)
    z
  }),
  targets::tar_target(prediction_simulation_summary,
    do.call(rbind, lapply(prediction_simulation_bundle, function(x) data.frame(
      architecture = x$simulation$scenario$architecture,
      replicate = x$simulation$scenario$replicate,
      simulation_seed = x$simulation$provenance$seed,
      causal_count = length(x$simulation$truth$causal$all),
      target_h2 = x$simulation$truth$parameters$h2_target,
      realized_h2 = x$simulation$truth$parameters$h2_observed,
      effect_distribution = x$simulation$extras$effect_distribution,
      component_counts = paste(names(table(x$simulation$extras$effect_components$component)),
        as.integer(table(x$simulation$extras$effect_components$component)), sep = ":", collapse = ";"),
      oracle_ok = x$oracle$ok, oracle_max_abs_error = x$oracle$max_abs_error,
      stringsAsFactors = FALSE)))),
  targets::tar_target(prediction_output_dir,
    file.path("results", "local", "02_prediction")),
  targets::tar_target(prediction_metrics_file,
    .study02_write_csv(prediction_metrics, file.path(prediction_output_dir,
      "prediction_metrics.csv")), format = "file"),
  targets::tar_target(prediction_paired_file,
    .study02_write_csv(prediction_paired_method_differences,
      file.path(prediction_output_dir, "paired_method_differences.csv")), format = "file"),
  targets::tar_target(prediction_computational_file,
    .study02_write_csv(prediction_computational_summary,
      file.path(prediction_output_dir, "computational_summary.csv")), format = "file"),
  targets::tar_target(prediction_status_file,
    .study02_write_csv(prediction_replicate_status,
      file.path(prediction_output_dir, "replicate_status.csv")), format = "file"),
  targets::tar_target(prediction_simulation_file,
    .study02_write_csv(prediction_simulation_summary,
      file.path(prediction_output_dir, "simulation_summary.csv")), format = "file"),
  targets::tar_target(prediction_manifest_file, {
    path <- file.path(prediction_output_dir, "prediction_manifest.json")
    jsonlite::write_json(list(study = prediction_config$study,
      task = prediction_config$task, development_settings = TRUE,
      active_architectures = names(prediction_config$simulation$architectures),
      active_methods = prediction_config$methods,
      multitrait = prediction_config$multitrait,
      replicate_count = length(prediction_replicate_specs) /
        length(prediction_config$simulation$architectures),
      expected_fit_count = length(prediction_replicate_specs) * length(prediction_config$methods),
      split = prediction_config$split,
      policies = list(genotype_scaling = "training individuals only",
        sparse_ld = "training individuals only",
        summary_statistics = "training individuals and phenotypes only",
        test_phenotypes = "evaluation only"),
      simulation = prediction_config$simulation, mcmc = prediction_config$mcmc,
      example_data = prediction_config$example_data,
      versions = list(sblr = as.character(packageVersion("sblr")),
        qgg = as.character(packageVersion("qgg")))), path,
      pretty = TRUE, auto_unbox = TRUE)
    path
  }, format = "file")
)
