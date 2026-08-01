targets::tar_option_set(packages = c("sblrbench", "sblr", "qgg",
  "posterior", "jsonlite", "Matrix"))

for (f in c("interface_audit.R", "state_contract.R", "simulation.R",
  "alignment.R", "operators.R", "methods.R", "chain_extraction.R",
  "diagnostics.R", "runtime_scaling.R", "metrics.R", "pilot.R",
  "promotion.R")) source(file.path("studies", "07_mtblr_validation", f),
    local = TRUE)

list(
  targets::tar_target(study07_config_file,
    file.path("studies", "07_mtblr_validation", "config.R"), format = "file"),
  targets::tar_target(study07_config,
    source(study07_config_file, local = TRUE)$value),
  targets::tar_target(study07_paths, .study07_paths(study07_config)),
  targets::tar_target(study07_resources,
    .study07_base_resources(study07_config)),
  targets::tar_target(study07_raw_genotypes,
    .study07_raw_genotypes(study07_resources, study07_config)),
  targets::tar_target(study07_interface_audit,
    .study07_interface_audit(study07_config)),
  targets::tar_target(study07_output_semantics,
    .study07_output_semantics()),
  targets::tar_target(study07_state_map, {
    .study07_validate_state_contract(study07_config)
    .study07_state_table(study07_config$trait_names)
  }),

  # 07A deterministic contract validation.
  targets::tar_target(study07_contract_scaled,
    .study07_scaled_data(study07_raw_genotypes, 2000L, 1:1000)),
  targets::tar_target(study07_contract_evidence,
    .study07_contract_evidence(study07_contract_scaled$train,
      study07_config)),
  targets::tar_target(study07_contract_simulations,
    lapply(study07_config$contract_architectures, function(a) {
      x <- .study07_simulate(study07_contract_scaled$train, a, 1L,
        study07_config)
      .study07_validate_simulation(x, study07_contract_scaled$train,
        study07_config)
      x
    }), iteration = "list"),
  targets::tar_target(study07_tiny_scaled,
    .study07_scaled_data(study07_raw_genotypes,
      study07_config$tiny$marker_count,
      seq_len(study07_config$tiny$sample_count))),
  targets::tar_target(study07_tiny_simulation, {
    x <- .study07_simulate(study07_tiny_scaled$train,
      "partially_shared", 1L, study07_config)
    .study07_validate_simulation(x, study07_tiny_scaled$train,
      study07_config)
    x
  }),
  targets::tar_target(study07_tiny_working_glist,
    .study07_working_glist(study07_resources$base_glist,
      study07_tiny_simulation$marker_ids,
      study07_tiny_scaled$allele_frequency, study07_config)),
  targets::tar_target(study07_tiny_ld_glist,
    .study07_make_ld(study07_tiny_working_glist,
      seq_len(study07_config$tiny$sample_count),
      study07_tiny_simulation$marker_ids, study07_config,
      study07_paths$ld_dir, "tiny")),
  targets::tar_target(study07_tiny_stats,
    .study07_make_stats(study07_tiny_simulation,
      study07_tiny_ld_glist, seq_len(study07_config$tiny$sample_count),
      study07_config)),
  targets::tar_target(study07_tiny_full_csr,
    sblr::sparseLD_read_CSR(study07_tiny_ld_glist$sparseLD$prefix,
      one_based = FALSE)),
  targets::tar_target(study07_tiny_operator,
    .study07_operator_bundle(study07_tiny_ld_glist, study07_tiny_stats,
      study07_tiny_full_csr, study07_config,
      file.path(study07_paths$operator_dir, "tiny"),
      study07_tiny_simulation$effects)),
  targets::tar_target(study07_tiny_runs, {
    controls <- list(nit = study07_config$tiny$nit,
      nburn = study07_config$tiny$nburn, nthin = 1L,
      nchains = study07_config$tiny$nchains,
      ncores = study07_config$tiny$ncores)
    lapply(study07_config$contract_implementations, function(id)
      .study07_cached_fit(id, study07_tiny_simulation,
        study07_tiny_stats, study07_tiny_ld_glist,
        study07_tiny_operator$runtime_glist,
        seq_len(study07_config$tiny$sample_count),
        study07_tiny_operator$blocks, study07_config, controls, "tiny"))
  }, iteration = "list"),
  targets::tar_target(study07_tiny_fit_status,
    do.call(rbind, lapply(study07_tiny_runs, function(x) data.frame(
      implementation = x$implementation$id, status = x$status,
      marker_count = nrow(x$fit$bm), trait_count = ncol(x$fit$bm),
      chain_count = length(x$fit$chains),
      covariance_history_available = any(
        x$fit$convergence_traces$quantities$group == "cov_g"),
      probability_history_available = any(
        x$fit$convergence_traces$quantities$group %in%
          c("pattern_pi", "pi_active")), warnings = paste(x$warnings,
            collapse = " | "), stringsAsFactors = FALSE)))),

  # 07B nested runtime scaling.
  targets::tar_target(study07_physical_memory,
    .study07_detect_memory()),
  targets::tar_target(study07_runtime_runs, {
    .study07_assert_phase_allowed("runtime")
    rows <- list()
    for (m in study07_config$marker_candidates) {
      scaled <- .study07_scaled_data(study07_raw_genotypes, m,
        seq_len(study07_config$sample_count))
      sim <- .study07_simulate(scaled$all, "partially_shared", 1L,
        study07_config)
      .study07_validate_simulation(sim, scaled$all, study07_config)
      glist <- .study07_working_glist(study07_resources$base_glist,
        sim$marker_ids, scaled$allele_frequency, study07_config)
      ld <- .study07_make_ld(glist, seq_len(study07_config$sample_count),
        sim$marker_ids, study07_config, study07_paths$ld_dir,
        paste0("runtime_", m))
      stats <- .study07_make_stats(sim, ld,
        seq_len(study07_config$sample_count), study07_config)
      csr <- sblr::sparseLD_read_CSR(ld$sparseLD$prefix, one_based = FALSE)
      operator <- .study07_operator_bundle(ld, stats, csr, study07_config,
        file.path(study07_paths$operator_dir, paste0("runtime_", m)),
        sim$effects)
      controls <- study07_config$timing[c("nit", "nburn", "nthin",
        "nchains", "ncores")]
      for (id in study07_config$runtime_implementations) {
        run <- .study07_cached_fit(id, sim, stats, ld,
          operator$runtime_glist, seq_len(study07_config$sample_count),
          operator$blocks, study07_config, controls, "runtime")
        rows[[length(rows) + 1L]] <- list(run = run,
          operator_summary = .study07_operator_summary(operator, m),
          marker_count = m)
      }
    }
    rows
  }, iteration = "list"),
  targets::tar_target(study07_runtime_table,
    do.call(rbind, lapply(study07_runtime_runs, function(x)
      .study07_runtime_row(x$run, x$marker_count,
        study07_config$sample_count, x$run$controls,
        study07_physical_memory)))),
  targets::tar_target(study07_runtime_projection,
    .study07_runtime_projection(study07_runtime_table, study07_config)),
  targets::tar_target(study07_selected_marker_count,
    .study07_select_marker_count(study07_runtime_projection,
      study07_config)),
  targets::tar_target(study07_runtime_operator_summaries,
    do.call(rbind, lapply(study07_runtime_runs, function(x)
      x$operator_summary))),
  targets::tar_target(study07_contract_output_files, {
    tables <- list(interface_audit.csv = study07_interface_audit,
      output_semantics.csv = study07_output_semantics,
      joint_state_map.csv = study07_state_map,
      permutation_contracts.csv = study07_contract_evidence,
      deterministic_simulation_summary.csv = do.call(rbind,
        lapply(study07_contract_simulations, .study07_simulation_summary)),
      tiny_fit_status.csv = study07_tiny_fit_status,
      operator_equivalence.csv = .study07_operator_summary(
        study07_tiny_operator, study07_config$tiny$marker_count),
      runtime_scaling.csv = study07_runtime_table,
      runtime_projections.csv = study07_runtime_projection,
      selected_marker_count.csv = study07_selected_marker_count,
      computational_limits.csv = as.data.frame(study07_config$runtime_limits,
        stringsAsFactors = FALSE))
    vapply(names(tables), function(name) .study07_write_csv(tables[[name]],
      file.path(study07_paths$contract_output, name)), "")
  }, format = "file"),

  # Shared selected-marker data for convergence and the main grid.
  targets::tar_target(study07_selected_scaled,
    .study07_scaled_data(study07_raw_genotypes,
      study07_selected_marker_count$marker_count,
      study07_resources$split$train_rows)),
  targets::tar_target(study07_selected_working_glist,
    .study07_working_glist(study07_resources$base_glist,
      colnames(study07_selected_scaled$all),
      study07_selected_scaled$allele_frequency, study07_config)),
  targets::tar_target(study07_selected_ld_glist,
    .study07_make_ld(study07_selected_working_glist,
      study07_resources$split$train_rows,
      colnames(study07_selected_scaled$all), study07_config,
      study07_paths$ld_dir, "selected_training")),
  targets::tar_target(study07_convergence_simulation, {
    x <- .study07_simulate(study07_selected_scaled$all,
      "partially_shared", 1L, study07_config)
    .study07_validate_simulation(x, study07_selected_scaled$all,
      study07_config)
    x
  }),
  targets::tar_target(study07_convergence_stats,
    .study07_make_stats(study07_convergence_simulation,
      study07_selected_ld_glist, study07_resources$split$train_rows,
      study07_config)),
  targets::tar_target(study07_selected_csr,
    sblr::sparseLD_read_CSR(study07_selected_ld_glist$sparseLD$prefix,
      one_based = FALSE)),
  targets::tar_target(study07_selected_operator,
    .study07_operator_bundle(study07_selected_ld_glist,
      study07_convergence_stats, study07_selected_csr, study07_config,
      file.path(study07_paths$operator_dir, "selected_training"),
      study07_convergence_simulation$effects)),

  # 07C maximum-history convergence.
  targets::tar_target(study07_convergence_runs, {
    .study07_assert_phase_allowed("convergence")
    controls <- list(nit = study07_config$convergence$maximum_nit,
      nburn = study07_config$convergence$maximum_nburn,
      nthin = 1L, nchains = 4L, ncores = 4L)
    lapply(study07_config$runtime_implementations, function(id)
      .study07_cached_fit(id, study07_convergence_simulation,
        study07_convergence_stats, study07_selected_ld_glist,
        study07_selected_operator$runtime_glist,
        study07_resources$split$train_rows,
        study07_selected_operator$blocks, study07_config, controls,
        "convergence"))
  }, iteration = "list"),
  targets::tar_target(study07_convergence_draws,
    lapply(study07_convergence_runs, function(x)
      .study07_extract_draws(x, "partially_shared", 1L)),
    iteration = "list"),
  targets::tar_target(study07_convergence_selection,
    lapply(study07_convergence_draws, function(x)
      .study07_select_recommendation(x, study07_config)),
    iteration = "list"),
  targets::tar_target(study07_convergence_recommendations,
    do.call(rbind, lapply(study07_convergence_selection, `[[`,
      "recommendation"))),
  targets::tar_target(study07_convergence_diagnostics,
    do.call(rbind, lapply(study07_convergence_selection, `[[`,
      "diagnostics"))),
  targets::tar_target(study07_convergence_candidates,
    do.call(rbind, lapply(study07_convergence_selection, `[[`,
      "candidates"))),
  targets::tar_target(study07_convergence_output_files, {
    tables <- list(method_recommendations.csv =
        study07_convergence_recommendations,
      convergence_diagnostics.csv = study07_convergence_diagnostics,
      candidate_settings.csv = study07_convergence_candidates,
      fit_status.csv = do.call(rbind, lapply(study07_convergence_runs,
        function(x) data.frame(implementation = x$implementation$id,
          status = x$status, elapsed_seconds = x$runtime,
          warnings = paste(x$warnings, collapse = " | "),
          error_message = x$error, stringsAsFactors = FALSE))),
      selected_marker_count.csv = study07_selected_marker_count)
    vapply(names(tables), function(name) .study07_write_csv(tables[[name]],
      file.path(study07_paths$convergence_output, name)), "")
  }, format = "file"),

  # 07D 30-fit five-replicate development benchmark.
  targets::tar_target(study07_main_runs, {
    .study07_assert_phase_allowed("benchmark")
    runs <- list()
    for (architecture in study07_config$main_architectures)
      for (replicate in seq_len(study07_config$replicate_count)) {
        sim <- .study07_simulate(study07_selected_scaled$all,
          architecture, replicate, study07_config)
        .study07_validate_simulation(sim, study07_selected_scaled$all,
          study07_config)
        stats <- .study07_make_stats(sim, study07_selected_ld_glist,
          study07_resources$split$train_rows, study07_config)
        for (id in study07_config$implementations) {
          rec <- study07_convergence_recommendations[
            study07_convergence_recommendations$implementation == id, ]
          if (nrow(rec) != 1L || rec$recommendation_status != "available")
            stop("Missing Study 07 operational recommendation: ", id)
          controls <- list(nit = rec$nit, nburn = rec$nburn,
            nthin = rec$nthin, nchains = rec$nchains, ncores = rec$ncores)
          run <- .study07_cached_fit(id, sim, stats,
            study07_selected_ld_glist,
            study07_selected_operator$runtime_glist,
            study07_resources$split$train_rows,
            study07_selected_operator$blocks, study07_config, controls,
            "benchmark")
          runs[[length(runs) + 1L]] <- list(run = run,
            simulation = sim, stats = stats)
        }
      }
    runs
  }, iteration = "list"),
  targets::tar_target(study07_prediction_metrics,
    do.call(rbind, lapply(study07_main_runs, function(x)
      .study07_prediction_metrics(x$run, x$simulation,
        study07_selected_scaled$test, study07_resources$split$test_rows)))),
  targets::tar_target(study07_marker_metrics,
    do.call(rbind, lapply(study07_main_runs, function(x)
      .study07_marker_metrics(x$run, x$simulation,
        study07_selected_scaled$all)))),
  targets::tar_target(study07_internal_identities,
    do.call(rbind, lapply(study07_main_runs, function(x)
      .study07_internal_consistency(x$run, x$simulation)))),
  targets::tar_target(study07_main_draws,
    lapply(study07_main_runs, function(x)
      .study07_extract_draws(x$run, x$simulation$architecture,
        x$simulation$replicate)), iteration = "list"),
  targets::tar_target(study07_parameter_estimates,
    do.call(rbind, Map(function(x, d)
      .study07_parameter_estimates(d, x$simulation),
      study07_main_runs, study07_main_draws))),
  targets::tar_target(study07_paired_differences,
    .study07_paired_differences(rbind(study07_prediction_metrics,
      study07_marker_metrics))),
  targets::tar_target(study07_paired_summary,
    .study07_paired_summary(study07_paired_differences)),
  targets::tar_target(study07_main_fit_status,
    do.call(rbind, lapply(study07_main_runs, function(x) data.frame(
      architecture = x$simulation$architecture,
      replicate = x$simulation$replicate,
      implementation = x$run$implementation$id,
      marker_count = length(x$simulation$marker_ids),
      chain_count = length(x$run$fit$chains), status = x$run$status,
      elapsed_seconds = x$run$runtime,
      warnings = paste(x$run$warnings, collapse = " | "),
      error_message = x$run$error, stringsAsFactors = FALSE)))),
  targets::tar_target(study07_main_output_files, {
    if (any(!study07_internal_identities$passed))
      stop("A completed Study 07 fit failed an internal identity.")
    tables <- list(
      fit_status.csv = study07_main_fit_status,
      simulation_truth.csv = do.call(rbind, lapply(
        study07_main_runs[!duplicated(vapply(study07_main_runs,
          function(x) paste(x$simulation$architecture,
            x$simulation$replicate), ""))],
        function(x) .study07_simulation_summary(x$simulation))),
      prediction_metrics.csv = study07_prediction_metrics,
      trait_marker_metrics.csv = study07_marker_metrics,
      parameter_estimates.csv = study07_parameter_estimates,
      internal_consistency.csv = study07_internal_identities,
      paired_replicate_differences.csv = study07_paired_differences,
      paired_comparison_summary.csv = study07_paired_summary,
      computational_summary.csv = study07_main_fit_status[, c(
        "architecture", "replicate", "implementation", "elapsed_seconds")],
      selected_marker_count.csv = study07_selected_marker_count)
    vapply(names(tables), function(name) .study07_write_csv(tables[[name]],
      file.path(study07_paths$benchmark_output, name)), "")
  }, format = "file")
)
