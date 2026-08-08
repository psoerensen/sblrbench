targets::tar_option_set(packages = c("sblrbench", "sblr", "qgg",
  "posterior", "jsonlite", "Matrix"))

for (f in c("interface_audit.R", "state_contract.R", "simulation.R",
  "alignment.R", "operators.R", "methods.R", "chain_extraction.R",
  "diagnostics.R", "runtime_scaling.R", "metrics.R", "pilot.R",
  "promotion.R")) source(file.path("studies", "08_mt_validation", f),
    local = TRUE)

list(
  targets::tar_target(study08_config_file,
    file.path("studies", "08_mt_validation", "config.R"), format = "file"),
  targets::tar_target(study08_config,
    source(study08_config_file, local = TRUE)$value),
  targets::tar_target(study08_paths, .study08_paths(study08_config)),
  targets::tar_target(study08_resources,
    .study08_base_resources(study08_config)),
  targets::tar_target(study08_raw_genotypes,
    .study08_raw_genotypes(study08_resources, study08_config)),
  targets::tar_target(study08_interface_audit,
    .study08_interface_audit(study08_config)),
  targets::tar_target(study08_output_semantics,
    .study08_output_semantics()),
  targets::tar_target(study08_state_map, {
    .study08_validate_state_contract(study08_config)
    .study08_state_table(study08_config$trait_names)
  }),

  # 08A deterministic contract validation.
  targets::tar_target(study08_contract_scaled,
    .study08_scaled_data(study08_raw_genotypes, 2000L, 1:1000)),
  targets::tar_target(study08_contract_evidence,
    .study08_contract_evidence(study08_contract_scaled$train,
      study08_config)),
  targets::tar_target(study08_contract_simulations,
    lapply(study08_config$contract_architectures, function(a) {
      x <- .study08_simulate(study08_contract_scaled$train, a, 1L,
        study08_config)
      .study08_validate_simulation(x, study08_contract_scaled$train,
        study08_config)
      x
    }), iteration = "list"),
  targets::tar_target(study08_tiny_scaled,
    .study08_scaled_data(study08_raw_genotypes,
      study08_config$tiny$marker_count,
      seq_len(study08_config$tiny$sample_count))),
  targets::tar_target(study08_tiny_simulation, {
    x <- .study08_simulate(study08_tiny_scaled$train,
      "partially_shared", 1L, study08_config)
    .study08_validate_simulation(x, study08_tiny_scaled$train,
      study08_config)
    x
  }),
  targets::tar_target(study08_tiny_working_glist,
    .study08_working_glist(study08_resources$base_glist,
      study08_tiny_simulation$marker_ids,
      study08_tiny_scaled$allele_frequency, study08_config)),
  targets::tar_target(study08_tiny_ld_glist,
    .study08_make_ld(study08_tiny_working_glist,
      seq_len(study08_config$tiny$sample_count),
      study08_tiny_simulation$marker_ids, study08_config,
      study08_paths$ld_dir, "tiny")),
  targets::tar_target(study08_tiny_stats,
    .study08_make_stats(study08_tiny_simulation,
      study08_tiny_ld_glist, seq_len(study08_config$tiny$sample_count),
      study08_config)),
  targets::tar_target(study08_tiny_full_csr,
    sblr::sparseLD_read_CSR(study08_tiny_ld_glist$sparseLD$prefix,
      one_based = FALSE)),
  targets::tar_target(study08_tiny_operator,
    .study08_operator_bundle(study08_tiny_ld_glist, study08_tiny_stats,
      study08_tiny_full_csr, study08_config,
      file.path(study08_paths$operator_dir, "tiny"),
      study08_tiny_simulation$effects)),
  targets::tar_target(study08_tiny_runs, {
    controls <- list(nit = study08_config$tiny$nit,
      nburn = study08_config$tiny$nburn, nthin = 1L,
      nchains = study08_config$tiny$nchains,
      ncores = study08_config$tiny$ncores)
    lapply(study08_config$contract_implementations, function(id)
      .study08_cached_fit(id, study08_tiny_simulation,
        study08_tiny_stats, study08_tiny_ld_glist,
        study08_tiny_operator$runtime_glist,
        seq_len(study08_config$tiny$sample_count),
        study08_tiny_operator$blocks, study08_config, controls, "tiny"))
  }, iteration = "list"),
  targets::tar_target(study08_tiny_fit_status,
    do.call(rbind, lapply(study08_tiny_runs, function(x) data.frame(
      implementation = x$implementation$id, status = x$status,
      marker_count = nrow(x$fit$bm), trait_count = ncol(x$fit$bm),
      chain_count = length(x$fit$chains),
      covariance_history_available = any(
        x$fit$convergence_traces$quantities$group == "cov_g"),
      probability_history_available = any(
        x$fit$convergence_traces$quantities$group %in%
          c("pattern_pi", "pi_active")), warnings = paste(x$warnings,
            collapse = " | "), stringsAsFactors = FALSE)))),

  # 08B nested runtime scaling.
  targets::tar_target(study08_physical_memory,
    .study08_detect_memory()),
  targets::tar_target(study08_runtime_runs, {
    .study08_assert_phase_allowed("runtime")
    rows <- list()
    for (m in study08_config$marker_candidates) {
      scaled <- .study08_scaled_data(study08_raw_genotypes, m,
        seq_len(study08_config$sample_count))
      sim <- .study08_simulate(scaled$all, "partially_shared", 1L,
        study08_config)
      .study08_validate_simulation(sim, scaled$all, study08_config)
      glist <- .study08_working_glist(study08_resources$base_glist,
        sim$marker_ids, scaled$allele_frequency, study08_config)
      ld <- .study08_make_ld(glist, seq_len(study08_config$sample_count),
        sim$marker_ids, study08_config, study08_paths$ld_dir,
        paste0("runtime_", m))
      stats <- .study08_make_stats(sim, ld,
        seq_len(study08_config$sample_count), study08_config)
      csr <- sblr::sparseLD_read_CSR(ld$sparseLD$prefix, one_based = FALSE)
      operator <- .study08_operator_bundle(ld, stats, csr, study08_config,
        file.path(study08_paths$operator_dir, paste0("runtime_", m)),
        sim$effects)
      controls <- study08_config$timing[c("nit", "nburn", "nthin",
        "nchains", "ncores")]
      for (id in study08_config$runtime_implementations) {
        run <- .study08_cached_fit(id, sim, stats, ld,
          operator$runtime_glist, seq_len(study08_config$sample_count),
          operator$blocks, study08_config, controls, "runtime")
        rows[[length(rows) + 1L]] <- list(run = run,
          operator_summary = .study08_operator_summary(operator, m),
          marker_count = m)
      }
    }
    rows
  }, iteration = "list"),
  targets::tar_target(study08_runtime_table,
    do.call(rbind, lapply(study08_runtime_runs, function(x)
      .study08_runtime_row(x$run, x$marker_count,
        study08_config$sample_count, x$run$controls,
        study08_physical_memory)))),
  targets::tar_target(study08_runtime_projection,
    .study08_runtime_projection(study08_runtime_table, study08_config)),
  targets::tar_target(study08_selected_marker_count,
    .study08_select_marker_count(study08_runtime_projection,
      study08_config)),
  targets::tar_target(study08_runtime_operator_summaries,
    do.call(rbind, lapply(study08_runtime_runs, function(x)
      x$operator_summary))),
  targets::tar_target(study08_contract_output_files, {
    tables <- list(interface_audit.csv = study08_interface_audit,
      output_semantics.csv = study08_output_semantics,
      joint_state_map.csv = study08_state_map,
      permutation_contracts.csv = study08_contract_evidence,
      deterministic_simulation_summary.csv = do.call(rbind,
        lapply(study08_contract_simulations, .study08_simulation_summary)),
      tiny_fit_status.csv = study08_tiny_fit_status,
      operator_equivalence.csv = .study08_operator_summary(
        study08_tiny_operator, study08_config$tiny$marker_count),
      runtime_scaling.csv = study08_runtime_table,
      runtime_projections.csv = study08_runtime_projection,
      selected_marker_count.csv = study08_selected_marker_count,
      computational_limits.csv = as.data.frame(study08_config$runtime_limits,
        stringsAsFactors = FALSE))
    vapply(names(tables), function(name) .study08_write_csv(tables[[name]],
      file.path(study08_paths$contract_output, name)), "")
  }, format = "file"),

  # Shared selected-marker data for convergence and the main grid.
  targets::tar_target(study08_selected_scaled,
    .study08_scaled_data(study08_raw_genotypes,
      study08_selected_marker_count$marker_count,
      study08_resources$split$train_rows)),
  targets::tar_target(study08_selected_working_glist,
    .study08_working_glist(study08_resources$base_glist,
      colnames(study08_selected_scaled$all),
      study08_selected_scaled$allele_frequency, study08_config)),
  targets::tar_target(study08_selected_ld_glist,
    .study08_make_ld(study08_selected_working_glist,
      study08_resources$split$train_rows,
      colnames(study08_selected_scaled$all), study08_config,
      study08_paths$ld_dir, "selected_training")),
  targets::tar_target(study08_convergence_simulation, {
    x <- .study08_simulate(study08_selected_scaled$all,
      "partially_shared", 1L, study08_config)
    .study08_validate_simulation(x, study08_selected_scaled$all,
      study08_config)
    x
  }),
  targets::tar_target(study08_convergence_stats,
    .study08_make_stats(study08_convergence_simulation,
      study08_selected_ld_glist, study08_resources$split$train_rows,
      study08_config)),
  targets::tar_target(study08_selected_csr,
    sblr::sparseLD_read_CSR(study08_selected_ld_glist$sparseLD$prefix,
      one_based = FALSE)),
  targets::tar_target(study08_selected_operator,
    .study08_operator_bundle(study08_selected_ld_glist,
      study08_convergence_stats, study08_selected_csr, study08_config,
      file.path(study08_paths$operator_dir, "selected_training"),
      study08_convergence_simulation$effects)),

  # 08C maximum-history convergence.
  targets::tar_target(study08_convergence_runs, {
    .study08_assert_phase_allowed("convergence")
    controls <- list(nit = study08_config$convergence$maximum_nit,
      nburn = study08_config$convergence$maximum_nburn,
      nthin = 1L, nchains = 4L, ncores = 4L)
    lapply(study08_config$runtime_implementations, function(id)
      .study08_cached_fit(id, study08_convergence_simulation,
        study08_convergence_stats, study08_selected_ld_glist,
        study08_selected_operator$runtime_glist,
        study08_resources$split$train_rows,
        study08_selected_operator$blocks, study08_config, controls,
        "convergence"))
  }, iteration = "list"),
  targets::tar_target(study08_convergence_draws,
    lapply(study08_convergence_runs, function(x)
      .study08_extract_draws(x, "partially_shared", 1L)),
    iteration = "list"),
  targets::tar_target(study08_convergence_selection,
    lapply(study08_convergence_draws, function(x)
      .study08_select_recommendation(x, study08_config)),
    iteration = "list"),
  targets::tar_target(study08_convergence_recommendations,
    do.call(rbind, lapply(study08_convergence_selection, `[[`,
      "recommendation"))),
  targets::tar_target(study08_convergence_diagnostics,
    do.call(rbind, lapply(study08_convergence_selection, `[[`,
      "diagnostics"))),
  targets::tar_target(study08_convergence_candidates,
    do.call(rbind, lapply(study08_convergence_selection, `[[`,
      "candidates"))),
  targets::tar_target(study08_convergence_output_files, {
    tables <- list(method_recommendations.csv =
        study08_convergence_recommendations,
      convergence_diagnostics.csv = study08_convergence_diagnostics,
      candidate_settings.csv = study08_convergence_candidates,
      fit_status.csv = do.call(rbind, lapply(study08_convergence_runs,
        function(x) data.frame(implementation = x$implementation$id,
          status = x$status, elapsed_seconds = x$runtime,
          warnings = paste(x$warnings, collapse = " | "),
          error_message = x$error, stringsAsFactors = FALSE))),
      selected_marker_count.csv = study08_selected_marker_count)
    vapply(names(tables), function(name) .study08_write_csv(tables[[name]],
      file.path(study08_paths$convergence_output, name)), "")
  }, format = "file"),

  # 08D 30-fit five-replicate development benchmark.
  targets::tar_target(study08_main_runs, {
    .study08_assert_phase_allowed("benchmark")
    runs <- list()
    for (architecture in study08_config$main_architectures)
      for (replicate in seq_len(study08_config$replicate_count)) {
        sim <- .study08_simulate(study08_selected_scaled$all,
          architecture, replicate, study08_config)
        .study08_validate_simulation(sim, study08_selected_scaled$all,
          study08_config)
        stats <- .study08_make_stats(sim, study08_selected_ld_glist,
          study08_resources$split$train_rows, study08_config)
        for (id in study08_config$implementations) {
          rec <- study08_convergence_recommendations[
            study08_convergence_recommendations$implementation == id, ]
          if (nrow(rec) != 1L || rec$recommendation_status != "available")
            stop("Missing Study 08 operational recommendation: ", id)
          controls <- list(nit = rec$nit, nburn = rec$nburn,
            nthin = rec$nthin, nchains = rec$nchains, ncores = rec$ncores)
          run <- .study08_cached_fit(id, sim, stats,
            study08_selected_ld_glist,
            study08_selected_operator$runtime_glist,
            study08_resources$split$train_rows,
            study08_selected_operator$blocks, study08_config, controls,
            "benchmark")
          runs[[length(runs) + 1L]] <- list(run = run,
            simulation = sim, stats = stats)
        }
      }
    runs
  }, iteration = "list"),
  targets::tar_target(study08_prediction_metrics,
    do.call(rbind, lapply(study08_main_runs, function(x)
      .study08_prediction_metrics(x$run, x$simulation,
        study08_selected_scaled$test, study08_resources$split$test_rows)))),
  targets::tar_target(study08_marker_metrics,
    do.call(rbind, lapply(study08_main_runs, function(x)
      .study08_marker_metrics(x$run, x$simulation,
        study08_selected_scaled$all)))),
  targets::tar_target(study08_internal_identities,
    do.call(rbind, lapply(study08_main_runs, function(x)
      .study08_internal_consistency(x$run, x$simulation)))),
  targets::tar_target(study08_main_draws,
    lapply(study08_main_runs, function(x)
      .study08_extract_draws(x$run, x$simulation$architecture,
        x$simulation$replicate)), iteration = "list"),
  targets::tar_target(study08_parameter_estimates,
    do.call(rbind, Map(function(x, d)
      .study08_parameter_estimates(d, x$simulation),
      study08_main_runs, study08_main_draws))),
  targets::tar_target(study08_paired_differences,
    .study08_paired_differences(rbind(study08_prediction_metrics,
      study08_marker_metrics))),
  targets::tar_target(study08_paired_summary,
    .study08_paired_summary(study08_paired_differences)),
  targets::tar_target(study08_main_fit_status,
    do.call(rbind, lapply(study08_main_runs, function(x) data.frame(
      architecture = x$simulation$architecture,
      replicate = x$simulation$replicate,
      implementation = x$run$implementation$id,
      marker_count = length(x$simulation$marker_ids),
      chain_count = length(x$run$fit$chains), status = x$run$status,
      elapsed_seconds = x$run$runtime,
      warnings = paste(x$run$warnings, collapse = " | "),
      error_message = x$run$error, stringsAsFactors = FALSE)))),
  targets::tar_target(study08_main_output_files, {
    if (any(!study08_internal_identities$passed))
      stop("A completed Study 08 fit failed an internal identity.")
    tables <- list(
      fit_status.csv = study08_main_fit_status,
      simulation_truth.csv = do.call(rbind, lapply(
        study08_main_runs[!duplicated(vapply(study08_main_runs,
          function(x) paste(x$simulation$architecture,
            x$simulation$replicate), ""))],
        function(x) .study08_simulation_summary(x$simulation))),
      prediction_metrics.csv = study08_prediction_metrics,
      trait_marker_metrics.csv = study08_marker_metrics,
      parameter_estimates.csv = study08_parameter_estimates,
      internal_consistency.csv = study08_internal_identities,
      paired_replicate_differences.csv = study08_paired_differences,
      paired_comparison_summary.csv = study08_paired_summary,
      computational_summary.csv = study08_main_fit_status[, c(
        "architecture", "replicate", "implementation", "elapsed_seconds")],
      selected_marker_count.csv = study08_selected_marker_count)
    vapply(names(tables), function(name) .study08_write_csv(tables[[name]],
      file.path(study08_paths$benchmark_output, name)), "")
  }, format = "file")
)
