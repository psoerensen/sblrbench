source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
source(file.path("studies", "02_prediction", "pilot.R"), local = TRUE)
source(file.path("studies", "02_prediction", "promotion.R"), local = TRUE)
source(file.path("studies", "five_replicate_helpers.R"), local = TRUE)

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
  targets::tar_target(prediction_benchmark_summary,
    .study02_benchmark_summary(prediction_metrics)),
  targets::tar_target(prediction_paired_comparison_summary,
    .study02_paired_summary(prediction_paired_method_differences)),
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
  targets::tar_target(prediction_seed_registry, do.call(rbind,
    lapply(prediction_method_run, function(x) do.call(rbind,
      lapply(seq_along(x$fit$chain_seeds), function(chain) data.frame(
        architecture = x$computational$architecture,
        replicate = x$computational$replicate, method = x$computational$method,
        data_selection_seed = prediction_config$split$seed,
        architecture_seed = prediction_config$simulation$base_seed +
          match(x$computational$architecture,
            names(prediction_config$simulation$architectures)) * 1000L,
        simulation_seed = x$computational$simulation_seed,
        fit_seed = x$fit$mcmc_seed, chain = chain,
        chain_seed = x$fit$chain_seeds[[chain]], stringsAsFactors = FALSE)))))),
  targets::tar_target(prediction_output_dir,
    Sys.getenv("SBLR_BENCH_OUTPUT_DIR", file.path("results", "local", "02_prediction"))),
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
  targets::tar_target(prediction_benchmark_summary_file,
    .study02_write_csv(prediction_benchmark_summary,
      file.path(prediction_output_dir, "benchmark_summary.csv")), format = "file"),
  targets::tar_target(prediction_paired_summary_file,
    .study02_write_csv(prediction_paired_comparison_summary,
      file.path(prediction_output_dir, "paired_comparison_summary.csv")), format = "file"),
  targets::tar_target(prediction_seed_file,
    .study02_write_csv(prediction_seed_registry,
      file.path(prediction_output_dir, "seed_registry.csv")), format = "file"),
  targets::tar_target(prediction_manifest_file, {
    path <- file.path(prediction_output_dir, "prediction_manifest.json")
    source_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    source_clean <- !length(system2("git", c("status", "--porcelain"), stdout = TRUE))
    sblr_provenance <- .five_replicate_sblr_provenance()
    nrep <- length(unique(prediction_simulation_summary$replicate))
    expected_fits <- 2L * nrep * length(prediction_config$methods)
    complete <- nrow(prediction_replicate_status) == expected_fits &&
      all(prediction_replicate_status$status == "ok")
    method_seeds <- unique(prediction_computational_summary[, c("architecture", "replicate", "method", "mcmc_seed")])
    jsonlite::write_json(list(study = prediction_config$study, task = prediction_config$task,
      benchmark_title = if (nrep == 5L) "Single-trait prediction five-replicate development benchmark" else "Single-trait prediction one-replicate development benchmark",
      benchmark_scope = if (nrep == 5L) "five_replicate_development" else "one_replicate_development",
      benchmark_status = if (complete) "complete" else "incomplete",
      development_settings = TRUE, replicate_count = nrep,
      expected_fit_count = expected_fits,
      successful_fit_count = sum(prediction_replicate_status$status == "ok"),
      failed_fit_count = sum(prediction_replicate_status$status != "ok"),
      architectures = names(prediction_config$simulation$architectures),
      active_methods = prediction_config$methods,
      train_fraction = prediction_config$split$train_fraction,
      training_sample_count = length(prediction_split$train_ids),
      test_sample_count = length(prediction_split$test_ids),
      canonical_marker_count = length(prediction_filtered_markers$marker_ids),
      causal_count = prediction_config$simulation$n_causal,
      target_h2 = prediction_config$simulation$h2,
      realized_h2 = prediction_simulation_summary[c("architecture", "replicate", "realized_h2")],
      simulation_seeds = prediction_simulation_summary[c("architecture", "replicate", "simulation_seed")],
      split_seed = prediction_split$split_seed,
      method_seeds = method_seeds, chain_seeds = prediction_seed_registry,
      mcmc_controls = if (nrep == 5L) .five_replicate_recommendations() else prediction_config$mcmc,
      bayesr_mixture_multipliers = prediction_config$priors$bayesr_mixture_var,
      bayesc_initial_pi = prediction_config$priors$bayesc_inclusion_probability,
      policies = list(genotype_scaling = "training individuals only",
        sparse_ld = "training individuals only",
        summary_statistics = "training individuals and phenotypes only",
        test_phenotypes = "evaluation only"),
      oracle_results = prediction_simulation_summary[, c("architecture", "oracle_ok", "oracle_max_abs_error")],
      qgdata = prediction_config$example_data,
      versions = list(sblr = as.character(packageVersion("sblr")),
        sblr_source_commit = sblr_provenance$commit,
        sblr_source_status = sblr_provenance$source_status,
        qgg = as.character(packageVersion("qgg")),
        sblrbench = as.character(packageVersion("sblrbench"))),
      sblrbench_source_commit = source_commit, source_tree_clean = source_clean,
      source_status = "working-tree five-replicate source; source inventory retained in capsule",
      completion_counts = list(expected = expected_fits,
        successful = sum(prediction_replicate_status$status == "ok"),
        failed = sum(prediction_replicate_status$status != "ok")),
      failures = prediction_replicate_status[prediction_replicate_status$status != "ok", ],
      provenance = prediction_config$example_data,
      validation_status = if (complete) "complete_grid_validated" else "incomplete_not_promotable",
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      thread_settings = .five_replicate_thread_settings(),
      limitations = if (nrep == 5L) c("five-replicate development evidence is not universal validation",
        "aggregate uncertainty is descriptive across five simulations") else
        c("one replicate per architecture", "short single-chain development MCMC",
          "no replicate-to-replicate uncertainty", "no convergence or general method-ranking claims")), path,
      pretty = TRUE, auto_unbox = TRUE)
    path
  }, format = "file")
)
