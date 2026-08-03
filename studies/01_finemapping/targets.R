source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
source(file.path("studies", "01_finemapping", "pilot.R"), local = TRUE)
source(file.path("studies", "five_replicate_helpers.R"), local = TRUE)

list(
  targets::tar_target(config_file, file.path("studies", "01_finemapping", "config.R"), format = "file"),
  targets::tar_target(config, source(config_file, local = TRUE)$value),
  targets::tar_target(data_paths, .study01_paths()),
  targets::tar_target(example_files, if (identical(data_paths$source, "qgg_example")) .study01_example_files(data_paths$data_dir, config$example_data) else character()),
  targets::tar_target(Glist, .study01_load_glist(data_paths, example_files)),
  targets::tar_target(filtered_markers, .study01_run_qc(Glist, config)),
  targets::tar_target(working_glist, .study01_set_rsids_ld(Glist, config$chr, filtered_markers$marker_ids)),
  targets::tar_target(sparse_ld_glist, .study01_make_sparse_ld(working_glist, filtered_markers, config, data_paths$output_dir)),
  targets::tar_target(selected_ids, .study01_selected_ids(sparse_ld_glist, config$sample_limit)),
  targets::tar_target(standardized_genotypes, .study01_extract_genotypes(sparse_ld_glist, config$chr, selected_ids, filtered_markers$marker_ids)),
  targets::tar_target(replicate_override, Sys.getenv("SBLR_BENCH_REPLICATES", ""), cue = targets::tar_cue(mode = "always")),
  targets::tar_target(replicate_specs, .study01_replicate_specs(config, replicate_override), iteration = "list"),
  targets::tar_target(replicate_spec, replicate_specs, pattern = map(replicate_specs), iteration = "list"),
  targets::tar_target(replicate_bundle, {
    selection <- select_separated_causal_markers(sparse_ld_glist, config$chr, filtered_markers$marker_ids,
      config$simulation$n_causal, config$causal_selection$min_distance_bp, replicate_spec$causal_seed,
      config$causal_selection$min_maf, config$causal_selection$max_maf)
    simulation <- .study01_simulate(replicate_spec, selection, standardized_genotypes, config)
    oracle <- sblrbench::check_oracle_genetic_values(simulation, tolerance = config$oracle_tolerance, stop_on_failure = TRUE)
    summary_stats <- .study01_summary_stats(simulation, sparse_ld_glist, config)
    list(spec = replicate_spec, selection = selection, simulation = simulation, oracle = oracle, stats = summary_stats)
  }, pattern = map(replicate_spec), iteration = "list"),
  targets::tar_target(method_specs, .study01_method_specs(config), iteration = "list"),
  targets::tar_target(method_spec, method_specs, pattern = map(method_specs), iteration = "list"),
  targets::tar_target(method_run, {
    fit <- .study01_fit(method_spec, replicate_bundle$simulation, replicate_bundle$stats, sparse_ld_glist, config)
    metrics <- .study01_marker_metrics(fit, replicate_bundle$simulation)
    cs <- .study01_credible_sets(fit, sparse_ld_glist, config)
    cs_metrics <- .study01_evaluate_cs(cs, fit, replicate_bundle$simulation, standardized_genotypes, sparse_ld_glist, config)
    list(fit = fit, marker_metrics = metrics, credible_sets = cs, credible_set_metrics = cs_metrics,
      computational = .study01_compact_fit(fit, replicate_bundle$simulation, config))
  }, pattern = cross(replicate_bundle, method_spec), iteration = "list"),
  targets::tar_target(computational_summary_branch, method_run$computational,
    pattern = map(method_run), iteration = "list"),
  targets::tar_target(marker_metrics_branch, {
    z <- method_run$marker_metrics
    names(z)[names(z) == "method_id"] <- "method"
    merge(z, method_run$computational,
      by = c("replicate", "method", "trait"), all.x = TRUE, sort = TRUE)
  }, pattern = map(method_run), iteration = "list"),
  targets::tar_target(credible_set_metrics_branch, {
    z <- method_run$credible_set_metrics
    names(z)[names(z) == "method_id"] <- "method"
    merge(z, method_run$computational,
      by = c("replicate", "method", "trait"), all.x = TRUE, sort = TRUE)
  }, pattern = map(method_run), iteration = "list"),
  targets::tar_target(credible_set_summary_branch, {
    z <- .study01_cs_summary(method_run$credible_set_metrics)
    names(z)[names(z) == "method_id"] <- "method"
    merge(z, method_run$computational,
      by = c("replicate", "method", "trait"), all.x = TRUE, sort = TRUE)
  }, pattern = map(method_run), iteration = "list"),
  targets::tar_target(computational_summary, do.call(rbind, computational_summary_branch)),
  targets::tar_target(marker_metrics, do.call(rbind, marker_metrics_branch)),
  targets::tar_target(credible_set_metrics, do.call(rbind, credible_set_metrics_branch)),
  targets::tar_target(credible_set_summary, do.call(rbind, credible_set_summary_branch)),
  targets::tar_target(replicate_status, {
    expected <- expand.grid(replicate = seq_along(replicate_specs), method = config$methods, stringsAsFactors = FALSE)
    z <- merge(expected, computational_summary[, c("replicate", "method", "status", "reason")], by = c("replicate", "method"), all.x = TRUE, sort = TRUE)
    z$status[is.na(z$status)] <- "missing"; z$reason[is.na(z$reason)] <- "required method branch is absent"; z
  }),
  targets::tar_target(simulation_summary, do.call(rbind, lapply(replicate_bundle, function(x) {
    pars <- x$simulation$truth$parameters
    data.frame(replicate = x$simulation$scenario$replicate,
      causal_selection_seed = x$spec$causal_seed,
      simulation_seed = x$spec$simulation_seed,
      causal_count = length(x$simulation$truth$causal$all),
      causal_markers = paste(x$simulation$truth$causal$all, collapse = ";"),
      target_h2 = config$simulation$h2,
      realized_h2 = pars$h2_observed %||% pars$h2_realized %||% NA_real_,
      oracle_ok = isTRUE(x$oracle$ok),
      oracle_max_abs_error = x$oracle$max_abs_error,
      stringsAsFactors = FALSE)
  }))),
  targets::tar_target(seed_registry, do.call(rbind, lapply(method_run, function(x)
    data.frame(replicate = x$computational$replicate,
      method = x$computational$method,
      causal_selection_seed = x$computational$causal_selection_seed,
      simulation_seed = x$computational$simulation_seed,
      fit_seed = x$computational$mcmc_seed,
      chain = 1L, chain_seed = x$computational$mcmc_seed,
      stringsAsFactors = FALSE)))),
  targets::tar_target(benchmark_summary,
    .study01_metric_summary(marker_metrics, credible_set_summary)),
  targets::tar_target(output_dir, Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
    file.path("results", "local", "01_finemapping", "separated"))),
  targets::tar_target(marker_metrics_file, .study01_write_csv(marker_metrics, file.path(output_dir, "marker_metrics.csv")), format = "file"),
  targets::tar_target(credible_set_metrics_file, .study01_write_csv(credible_set_metrics, file.path(output_dir, "credible_set_metrics.csv")), format = "file"),
  targets::tar_target(credible_set_summary_file, .study01_write_csv(credible_set_summary, file.path(output_dir, "credible_set_summary.csv")), format = "file"),
  targets::tar_target(computational_summary_file, .study01_write_csv(computational_summary, file.path(output_dir, "computational_summary.csv")), format = "file"),
  targets::tar_target(replicate_status_file, .study01_write_csv(replicate_status, file.path(output_dir, "replicate_status.csv")), format = "file"),
  targets::tar_target(simulation_summary_file, .study01_write_csv(simulation_summary, file.path(output_dir, "simulation_summary.csv")), format = "file"),
  targets::tar_target(seed_registry_file, .study01_write_csv(seed_registry, file.path(output_dir, "seed_registry.csv")), format = "file"),
  targets::tar_target(benchmark_summary_file, .study01_write_csv(benchmark_summary, file.path(output_dir, "benchmark_summary.csv")), format = "file"),
  targets::tar_target(pilot_manifest_file, {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(output_dir, "pilot_manifest.json")
    sblr_provenance <- .five_replicate_sblr_provenance()
    expected_fits <- length(replicate_specs) * length(config$methods)
    complete <- nrow(replicate_status) == expected_fits && all(replicate_status$status == "ok")
    jsonlite::write_json(list(study = config$study, architecture = config$architecture,
      development_settings = TRUE, replicate_count = length(replicate_specs), methods = config$methods,
      simulation = config$simulation, causal_selection = config$causal_selection, mcmc = config$mcmc,
      credible_sets = config$credible_sets, package_versions = list(sblr = as.character(utils::packageVersion("sblr")),
      sblrbench = as.character(utils::packageVersion("sblrbench")), qgg = as.character(utils::packageVersion("qgg"))),
      sblr_source_commit = sblr_provenance$commit,
      sblr_source_status = sblr_provenance$source_status,
      sblrbench_source_commit = sblrbench::sblrbench_git_commit("."),
      canonical_marker_count = length(filtered_markers$marker_ids),
      sample_count = length(selected_ids), expected_fit_count = expected_fits,
      successful_fit_count = sum(replicate_status$status == "ok"),
      failed_fit_count = sum(replicate_status$status != "ok"),
      validation_status = if (complete) "complete_grid_validated" else "incomplete_not_promotable",
      qgdata = config$example_data,
      execution_policy = "established Study 01 one-chain fine-mapping policy retained explicitly"),
      path, pretty = TRUE, auto_unbox = TRUE, null = "null")
    path
  }, format = "file")
)
