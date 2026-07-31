targets::tar_option_set(packages = c("sblrbench", "sblr", "qgg",
  "posterior", "jsonlite", "Matrix"))

source(file.path("studies", "01_finemapping",
  "setup_example_data.R"), local = TRUE)
source(file.path("studies", "02_prediction", "pilot.R"), local = TRUE)
for (f in c("blocks.R", "operators.R", "operator_validation.R",
            "simulation.R", "methods.R", "chain_extraction.R",
            "diagnostics.R", "metrics.R", "pilot.R"))
  source(file.path("studies", "06_ld_operator", f), local = TRUE)

list(
  targets::tar_target(study06_config_file,
    file.path("studies", "06_ld_operator", "config.R"),
    format = "file"),
  targets::tar_target(study06_config,
    source(study06_config_file, local = TRUE)$value),
  targets::tar_target(study06_paths, .study06_paths(study06_config)),
  targets::tar_target(study06_example_files,
    .study01_example_files(study06_paths$data_dir,
      study06_config$example_data)),
  targets::tar_target(study06_base_glist,
    .study01_load_glist(study06_paths, study06_example_files)),
  targets::tar_target(study06_filtered_markers,
    .study01_run_qc(study06_base_glist, study06_config)),
  targets::tar_target(study06_sample_ids,
    .study01_selected_ids(study06_base_glist,
      study06_config$sample_limit)),
  targets::tar_target(study06_split,
    sblrbench::make_prediction_split(study06_sample_ids,
      study06_config$split$train_fraction,
      study06_config$split$seed)),
  targets::tar_target(study06_scaled_genotypes,
    .study06_load_scaled_genotypes(study06_base_glist,
      study06_config$chr, study06_sample_ids,
      study06_filtered_markers$marker_ids, study06_split)),
  targets::tar_target(study06_working_glist,
    .study02_set_training_af(.study01_set_rsids_ld(
      study06_base_glist, study06_config$chr,
      study06_filtered_markers$marker_ids),
      study06_config$chr, study06_filtered_markers$marker_ids,
      study06_scaled_genotypes$allele_frequency)),
  targets::tar_target(study06_ld_bundle, {
    started <- proc.time()[["elapsed"]]
    x <- .study02_make_ld(study06_working_glist, study06_split,
      study06_filtered_markers$marker_ids, study06_config,
      study06_paths$genotype_output_dir)
    list(Glist = x,
      elapsed_seconds = proc.time()[["elapsed"]] - started)
  }),
  targets::tar_target(study06_simulation_specs,
    .study06_simulation_specs(study06_config),
    iteration = "list"),
  targets::tar_target(study06_simulation_spec,
    study06_simulation_specs, pattern = map(study06_simulation_specs),
    iteration = "list"),
  targets::tar_target(study06_simulation, {
    x <- .study06_simulate(study06_simulation_spec,
      study06_scaled_genotypes$all, study06_config)
    .study06_validate_simulation(x, study06_scaled_genotypes$all,
      study06_config)
    x
  }, pattern = map(study06_simulation_spec), iteration = "list"),
  targets::tar_target(study06_stats_bundle,
    .study06_summary_stats(study06_simulation,
      study06_ld_bundle$Glist, study06_split, study06_config),
    pattern = map(study06_simulation), iteration = "list"),

  # Deterministic operator and block-selection gate.
  targets::tar_target(study06_operator_simulation, {
    spec <- .study06_simulation_specs(study06_config)[[1L]]
    x <- .study06_simulate(spec, study06_scaled_genotypes$all,
      study06_config)
    .study06_validate_simulation(x, study06_scaled_genotypes$all,
      study06_config)
    x
  }),
  targets::tar_target(study06_operator_stats,
    .study06_summary_stats(study06_operator_simulation,
      study06_ld_bundle$Glist, study06_split,
      study06_config)$stats),
  targets::tar_target(study06_full_csr,
    sblr::sparseLD_read_CSR(
      study06_ld_bundle$Glist$sparseLD$prefix,
      one_based = FALSE)),
  targets::tar_target(study06_block_candidate_table, {
    base <- .study06_block_candidates(
      study06_filtered_markers$marker_ids, study06_config)
    cross <- do.call(rbind, lapply(study06_config$block_candidates,
      function(size) .study06_cross_block_summary(study06_full_csr,
        length(study06_filtered_markers$marker_ids), size)))
    construction <- do.call(rbind, lapply(
      study06_config$block_candidates, function(size) {
        blocks <- .study06_blocks(
          study06_filtered_markers$marker_ids, size)
        started <- proc.time()[["elapsed"]]
        inspect <- .study06_inspect_operator(
          study06_ld_bundle$Glist, study06_operator_stats, blocks,
          filter = "ridge_fixed", eta = 0)
        dense <- .study06_dense_blocks(inspect)
        eigenvalues <- unlist(lapply(dense, function(A)
          eigen(A, symmetric = TRUE, only.values = TRUE)$values),
          use.names = FALSE)
        data.frame(target_size = size,
          operator_construction_and_eigendecomposition_seconds =
            proc.time()[["elapsed"]] - started,
          total_positive_eigenvalue_mass =
            sum(pmax(eigenvalues, 0)),
          negative_eigenvalue_count =
            sum(eigenvalues < -1e-10),
          negative_eigenvalue_mass =
            sum(abs(eigenvalues[eigenvalues < -1e-10])),
          runtime_diagonal_minimum = min(inspect$diagonal),
          runtime_diagonal_maximum = max(inspect$diagonal),
          gate_construction_status = "passed",
          stringsAsFactors = FALSE)
      }))
    Reduce(function(x, y) merge(x, y, by = "target_size",
      sort = FALSE), list(base, cross, construction))
  }),
  targets::tar_target(study06_selected_block_design,
    .study06_select_block_design(study06_block_candidate_table,
      study06_config)),
  targets::tar_target(study06_blocks, {
    b <- .study06_blocks(study06_filtered_markers$marker_ids,
      study06_selected_block_design$target_size)
    .study06_validate_blocks(b,
      study06_filtered_markers$marker_ids)
    b
  }),
  targets::tar_target(study06_unfiltered_operator,
    .study06_inspect_operator(study06_ld_bundle$Glist,
      study06_operator_stats, study06_blocks,
      filter = "ridge_fixed", eta = 0)),
  targets::tar_target(study06_runtime_csr,
    .study06_write_runtime_csr(study06_unfiltered_operator,
      file.path(study06_paths$operator_dir,
        paste0("runtime_matched_block_", 
          study06_selected_block_design$target_size)),
      study06_filtered_markers$marker_ids)),
  targets::tar_target(study06_runtime_glist,
    .study06_runtime_glist(study06_ld_bundle$Glist,
      study06_operator_stats, study06_runtime_csr)),
  targets::tar_target(study06_runtime_csr_inspection,
    .study06_inspect_from_csr(study06_runtime_csr$prefix,
      study06_unfiltered_operator$diagonal, study06_blocks,
      study06_unfiltered_operator)),
  targets::tar_target(study06_equivalence,
    .study06_equivalence_gate(study06_runtime_csr_inspection,
      study06_unfiltered_operator, study06_blocks,
      study06_config, action_vectors = list(
        summary_score = unlist(study06_operator_stats$wy,
          use.names = FALSE),
        true_effect = as.numeric(study06_operator_simulation$effects)))),

  targets::tar_target(study06_synthetic_filter_validation,
    .study06_synthetic_filter_validation()),

  # Hard-filter and LW pilot. No samplers run in these targets.
  targets::tar_target(study06_hard_filter_pilot, {
    rows <- list()
    for (tau in study06_config$hard_tau_candidates) {
      hard <- .study06_inspect_operator(study06_ld_bundle$Glist,
        study06_operator_stats, study06_blocks,
        filter = "hard_truncate", tau = tau)
      rows[[length(rows) + 1L]] <- .study06_filter_candidate(
        study06_unfiltered_operator, hard, study06_blocks,
        study06_config, tau)
    }
    do.call(rbind, rows)
  }),
  targets::tar_target(study06_selected_hard_filter,
    .study06_select_hard_filter(study06_hard_filter_pilot,
      study06_config$selected_hard_tau)),
  targets::tar_target(study06_selected_hard_operator,
    .study06_inspect_operator(study06_ld_bundle$Glist,
      study06_operator_stats, study06_blocks,
      filter = "hard_truncate",
      tau = study06_selected_hard_filter$eigen_tau)),
  targets::tar_target(study06_ridge_lw_operator,
    .study06_inspect_operator(study06_ld_bundle$Glist,
      study06_operator_stats, study06_blocks,
      filter = "ridge_lw")),
  targets::tar_target(study06_ridge_lw_audit,
    .study06_lw_audit(study06_unfiltered_operator,
      study06_ridge_lw_operator, study06_blocks, study06_config,
      length(study06_split$train_ids))),
  targets::tar_target(study06_fixed_ridge_pilot, {
    rows <- lapply(study06_config$ridge_fixed_shrinkage_candidates,
      function(a) {
        candidate <- .study06_inspect_operator(study06_ld_bundle$Glist,
          study06_operator_stats, study06_blocks, filter = "ridge_fixed",
          eta = a / (1 - a))
        .study06_fixed_ridge_candidate(study06_unfiltered_operator,
          candidate, study06_blocks, study06_config, a)
      })
    do.call(rbind, rows)
  }),
  targets::tar_target(study06_selected_ridge_filter,
    .study06_select_fixed_ridge(study06_fixed_ridge_pilot,
      study06_config$selected_ridge_shrinkage)),
  targets::tar_target(study06_selected_ridge_operator,
    .study06_inspect_operator(study06_ld_bundle$Glist,
      study06_operator_stats, study06_blocks, filter = "ridge_fixed",
      eta = study06_selected_ridge_filter$eigen_eta)),
  targets::tar_target(study06_operator_perturbation_metrics, {
    vectors <- list(summary_score = unlist(study06_operator_stats$wy,
      use.names = FALSE),
      true_effect = as.numeric(study06_operator_simulation$effects))
    candidates <- list(block_eigen_unfiltered = study06_unfiltered_operator,
      block_eigen_hard = study06_selected_hard_operator,
      block_eigen_ridge_fixed = study06_selected_ridge_operator,
      block_eigen_ridge_lw_sensitivity = study06_ridge_lw_operator)
    block_rows <- do.call(rbind, lapply(names(candidates), function(id)
      .study06_operator_metrics(study06_unfiltered_operator,
        candidates[[id]], study06_blocks, study06_config,
        operator_id = id, reference_id = "block_eigen_unfiltered",
        action_vectors = vectors)))
    full_actions <- do.call(rbind, lapply(names(vectors), function(id) {
      x <- vectors[[id]]
      full <- .study06_apply_csr_crossproduct(study06_full_csr,
        unlist(study06_operator_stats$ww, use.names = FALSE), x)
      block <- .study06_apply_blocks(
        .study06_dense_blocks(study06_unfiltered_operator), x)
      .study06_action_comparison(full, block, id, "full_csr",
        "runtime_matched_block_csr", x)
    }))
    list(block_metrics = block_rows, full_csr_actions = full_actions)
  }),
  targets::tar_target(study06_operator_output_dir,
    file.path(study06_paths$local_dir, "operator_output")),
  targets::tar_target(study06_operator_files, {
    tables <- list(
      block_design_candidates.csv = study06_block_candidate_table,
      selected_block_design.csv = study06_selected_block_design,
      selected_block_definitions.csv = study06_blocks,
      operator_equivalence_summary.csv =
        study06_equivalence$summary,
      block_level_operator_diagnostics.csv =
        study06_equivalence$block_metrics,
      hard_filter_candidates.csv = study06_hard_filter_pilot,
      selected_hard_filter.csv = study06_selected_hard_filter,
      synthetic_filter_validation.csv = study06_synthetic_filter_validation,
      ridge_lw_diagnostics.csv = study06_ridge_lw_operator$diagnostics,
      ridge_lw_audit.csv = study06_ridge_lw_audit,
      fixed_ridge_candidates.csv = study06_fixed_ridge_pilot,
      selected_ridge_filter.csv = study06_selected_ridge_filter,
      operator_perturbation_metrics.csv =
        study06_operator_perturbation_metrics$block_metrics,
      full_csr_operator_action_metrics.csv =
        study06_operator_perturbation_metrics$full_csr_actions)
    vapply(names(tables), function(name)
      .study06_write_csv(tables[[name]],
        file.path(study06_operator_output_dir, name)), "")
  }, format = "file"),

  # Ten operational-setting one-replicate operator pilot fits.
  targets::tar_target(study06_operator_pilot_runs, {
    grid <- .study06_method_grid(study06_config)
    grid <- grid[grid$configuration != "bed", , drop = FALSE]
    runs <- list()
    for (architecture in study06_config$architectures) {
      simulation <- study06_operator_simulation
      if (!identical(simulation$architecture, architecture)) {
        spec <- Filter(function(x) x$architecture == architecture &&
          x$replicate == 1L,
          .study06_simulation_specs(study06_config))[[1L]]
        simulation <- .study06_simulate(spec,
          study06_scaled_genotypes$all, study06_config)
      }
      stats <- .study06_summary_stats(simulation,
        study06_ld_bundle$Glist, study06_split, study06_config)$stats
      methods <- grid[grid$architecture == architecture, , drop = FALSE]
      for (i in seq_len(nrow(methods)))
        runs[[length(runs) + 1L]] <- .study06_run_cached_fit(
          as.list(methods[i, , drop = FALSE]), simulation, stats,
          study06_ld_bundle$Glist, study06_split, study06_blocks,
          study06_runtime_glist, study06_config, "pilot", NULL,
          study06_selected_hard_filter$eigen_tau,
          study06_selected_ridge_filter$eigen_eta)
    }
    runs
  }, iteration = "list"),
  targets::tar_target(study06_operator_pilot_summary,
    .study06_operator_pilot_summary(study06_operator_pilot_runs,
      study06_scaled_genotypes$test, study06_split,
      study06_scaled_genotypes$all,
      lapply(study06_config$architectures, function(a) {
        spec <- Filter(function(x) x$architecture == a && x$replicate == 1L,
          .study06_simulation_specs(study06_config))[[1L]]
        .study06_simulate(spec, study06_scaled_genotypes$all, study06_config)
      }))),
  targets::tar_target(study06_operator_pilot_file,
    .study06_write_csv(study06_operator_pilot_summary,
      file.path(study06_operator_output_dir,
        "one_replicate_operator_pilot.csv")), format = "file"),

  # Eight maximum-history block-operator convergence fits.
  targets::tar_target(study06_convergence_runs, {
    grid <- .study06_method_grid(study06_config)
    grid <- grid[grid$configuration %in% c("block_csr",
      "block_eigen_unfiltered", "block_eigen_hard",
      "block_eigen_ridge_fixed"), , drop = FALSE]
    runs <- list()
    for (architecture in study06_config$architectures) {
      spec <- Filter(function(x) x$architecture == architecture &&
        x$replicate == 1L,
        .study06_simulation_specs(study06_config))[[1L]]
      simulation <- .study06_simulate(spec,
        study06_scaled_genotypes$all, study06_config)
      stats <- .study06_summary_stats(simulation,
        study06_ld_bundle$Glist, study06_split,
        study06_config)$stats
      methods <- grid[grid$architecture == architecture, ,
        drop = FALSE]
      for (i in seq_len(nrow(methods))) {
        method <- as.list(methods[i, , drop = FALSE])
        runs[[length(runs) + 1L]] <- .study06_run_cached_fit(
          method, simulation, stats, study06_ld_bundle$Glist,
          study06_split, study06_blocks, study06_runtime_glist,
          study06_config, "convergence", NULL,
          study06_selected_hard_filter$eigen_tau,
          study06_selected_ridge_filter$eigen_eta)
      }
    }
    runs
  }, iteration = "list"),
  targets::tar_target(study06_convergence_selections, {
    lapply(study06_convergence_runs, function(run) {
      draws <- .study06_extract_draws(run)
      .study06_select_recommendation(draws, study06_config)
    })
  }, iteration = "list"),
  targets::tar_target(study06_convergence_recommendations,
    do.call(rbind, lapply(study06_convergence_selections,
      `[[`, "recommendation"))),
  targets::tar_target(study06_convergence_diagnostics,
    do.call(rbind, lapply(study06_convergence_selections,
      `[[`, "diagnostics"))),
  targets::tar_target(study06_convergence_candidates,
    do.call(rbind, lapply(study06_convergence_selections,
      `[[`, "candidates"))),
  targets::tar_target(study06_convergence_status,
    do.call(rbind, lapply(study06_convergence_runs, function(x)
      data.frame(architecture = x$architecture,
        replicate = x$replicate,
        configuration = x$method$configuration,
        method = x$method$native_method,
        status = x$status, error_message = x$reason,
        chain_count = if (x$status == "ok") length(x$fit$chains) else 0L,
        elapsed_seconds = x$runtime, warnings = x$warnings,
        stringsAsFactors = FALSE)))),
  targets::tar_target(study06_convergence_output_dir,
    file.path(study06_paths$local_dir, "convergence_output")),
  targets::tar_target(study06_convergence_files, {
    seed <- .study06_seed_registry(study06_config)
    seed <- seed[seed$replicate == 1L &
      seed$configuration %in% c("block_csr",
        "block_eigen_unfiltered", "block_eigen_hard",
        "block_eigen_ridge_fixed"), ]
    tables <- list(
      fit_status.csv = study06_convergence_status,
      convergence_diagnostics.csv =
        study06_convergence_diagnostics,
      candidate_settings.csv = study06_convergence_candidates,
      method_recommendations.csv =
        study06_convergence_recommendations,
      seed_registry.csv = seed,
      computational_summary.csv = study06_convergence_status)
    vapply(names(tables), function(name)
      .study06_write_csv(tables[[name]],
        file.path(study06_convergence_output_dir, name)), "")
  }, format = "file"),

  # Main 60-fit benchmark, cached per validated fit.
  targets::tar_target(study06_block_recommendations, {
    path <- file.path(study06_config$convergence_capsule,
      "method_recommendations.csv")
    if (!file.exists(path))
      stop("Validated Study 06 convergence capsule is required.")
    utils::read.csv(path, stringsAsFactors = FALSE)
  }),
  targets::tar_target(study06_benchmark_run_group, {
    methods <- .study06_method_grid(study06_config)
    methods <- methods[methods$architecture ==
      study06_simulation$architecture, , drop = FALSE]
    lapply(seq_len(nrow(methods)), function(i)
      .study06_run_cached_fit(
        as.list(methods[i, , drop = FALSE]), study06_simulation,
        study06_stats_bundle$stats, study06_ld_bundle$Glist,
        study06_split, study06_blocks, study06_runtime_glist,
        study06_config, "benchmark",
        study06_block_recommendations,
        study06_selected_hard_filter$eigen_tau,
        study06_selected_ridge_filter$eigen_eta))
  }, pattern = map(study06_simulation, study06_stats_bundle),
  iteration = "list"),
  targets::tar_target(study06_benchmark_runs,
    unlist(study06_benchmark_run_group, recursive = FALSE),
    iteration = "list"),
  targets::tar_target(study06_benchmark_outputs, {
    runs <- study06_benchmark_runs
    simulations <- study06_simulation
    find_sim <- function(run) Filter(function(x)
      x$architecture == run$architecture &&
        x$replicate == run$replicate, simulations)[[1L]]
    draws <- lapply(runs, .study06_extract_draws)
    prediction <- do.call(rbind, Map(function(run, sim)
      .study06_prediction_metrics(run, sim,
        study06_scaled_genotypes$test, study06_split),
      runs, lapply(runs, find_sim)))
    parameter <- do.call(rbind, Map(function(x, sim)
      .study06_parameter_estimates(x, sim), draws,
      lapply(runs, find_sim)))
    marker <- do.call(rbind, Map(function(run, sim)
      .study06_marker_metrics(run, sim,
        study06_scaled_genotypes$all),
      runs, lapply(runs, find_sim)))
    diagnostics <- do.call(rbind, Map(function(x, run)
      .study06_diagnostics(x, 0L, run$controls$nit,
        study06_config), draws, runs))
    status <- do.call(rbind, lapply(runs, function(x)
      data.frame(architecture = x$architecture,
        replicate = x$replicate,
        configuration = x$method$configuration,
        method = x$method$native_method,
        status = x$status, error_message = x$reason,
        chain_count = length(x$fit$chains),
        warnings = x$warnings, stringsAsFactors = FALSE)))
    computational <- do.call(rbind, lapply(runs, function(x)
      data.frame(architecture = x$architecture,
        replicate = x$replicate,
        configuration = x$method$configuration,
        method = x$method$native_method,
        preprocessing_time = study06_ld_bundle$elapsed_seconds,
        summary_statistic_time = {
          z <- Filter(function(s) s$architecture == x$architecture &&
            s$replicate == x$replicate, simulations)[[1L]]
          0
        },
        ld_construction_time = study06_ld_bundle$elapsed_seconds,
        model_fitting_time = x$runtime,
        prediction_time = NA_real_, metric_generation_time = NA_real_,
        total_time = x$runtime,
        output_size_bytes = as.numeric(object.size(x$fit)),
        memory_bytes_estimate = if (is.null(
          x$fit$memory_estimate$estimated_total_bytes)) NA_real_ else
            as.numeric(x$fit$memory_estimate$estimated_total_bytes),
        warnings = x$warnings, failures = x$reason,
        stringsAsFactors = FALSE)))
    sbayesr <- do.call(rbind, lapply(runs[
      vapply(runs, function(x)
        x$architecture == "sparse_mixture", logical(1))],
      .study06_sbayesr_evidence))
    marker_agreement <- .study06_marker_agreement(runs)
    combined_metrics <- rbind(
      prediction[c("architecture", "replicate", "configuration",
        "method", "metric", "value")],
      marker[c("architecture", "replicate", "configuration",
        "method", "metric", "value")])
    paired <- .study06_paired_differences(combined_metrics)
    list(runs = runs, prediction = prediction,
      parameter = parameter, marker = marker,
      diagnostics = diagnostics, status = status,
      computational = computational, sbayesr = sbayesr,
      marker_agreement = marker_agreement,
      paired = paired, paired_summary = .study06_paired_summary(paired),
      recovery_summary = .study06_recovery_summary(parameter),
      convergence_validation =
        .study06_convergence_validation_summary(diagnostics))
  }, iteration = "list"),
  targets::tar_target(study06_benchmark_output_dir,
    file.path(study06_paths$local_dir, "benchmark_output")),
  targets::tar_target(study06_benchmark_files, {
    x <- study06_benchmark_outputs
    simulations <- study06_simulation
    availability <- expand.grid(
      configuration = study06_config$configurations,
      output = c("marker_effects", "nonnull_probabilities",
        "component_probabilities", "vgs", "ves", "vbs",
        "vld", "vle", "log_cpo", "mean_log_cpo"),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    availability$available <- availability$output %in% c(
      "marker_effects", "nonnull_probabilities", "vgs", "ves",
      "vbs", "vld", "vle")
    availability$reason <- ifelse(availability$available, "",
      "not returned by the installed public fit contract")
    operator_summary <- rbind(
      transform(study06_equivalence$summary,
        operator_id = "block_eigen_unfiltered"),
      data.frame(gate = "hard_filter_selected",
        pass = study06_selected_hard_filter$pass,
        marker_mapping_pass = TRUE, block_coverage_pass = TRUE,
        diagonal_maximum_absolute_error = NA_real_,
        reconstruction_maximum_absolute_error =
          study06_selected_hard_filter$operator_frobenius_maximum_error,
        reconstruction_maximum_relative_error = NA_real_,
        matrix_vector_maximum_error =
          study06_selected_hard_filter$matrix_vector_maximum_error,
        quadratic_form_maximum_error =
          study06_selected_hard_filter$quadratic_form_maximum_error,
        matrix_error_block = NA_character_,
        relative_error_block = NA_character_,
        matrix_vector_error_block = NA_character_,
        quadratic_form_error_block = NA_character_,
        float32_reconstruction_contribution =
          study06_selected_hard_filter$operator_frobenius_maximum_error,
        absolute_tolerance = study06_config$operator_tolerance$absolute,
        relative_tolerance = study06_config$operator_tolerance$relative,
        product_absolute_tolerance =
          study06_config$operator_tolerance$product_absolute,
        quadratic_absolute_tolerance =
          study06_config$operator_tolerance$quadratic_absolute,
        rationale = "selected hard filter relative to unfiltered",
        operator_id = "block_eigen_hard"))
    tables <- list(
      selected_block_definitions.csv = study06_blocks,
      selected_block_design.csv = study06_selected_block_design,
      filter_specification.csv = rbind(
        data.frame(configuration = "block_eigen_unfiltered",
          policy = "ridge_fixed", tau = NA_real_, eta = 0),
        data.frame(configuration = "block_eigen_hard",
          policy = "hard_truncate",
          tau = study06_selected_hard_filter$eigen_tau,
          eta = NA_real_),
        data.frame(configuration = "block_eigen_ridge_fixed",
          policy = "ridge_fixed", tau = NA_real_,
          eta = study06_selected_ridge_filter$eigen_eta)),
      seed_registry.csv = .study06_seed_registry(study06_config),
      simulation_summary.csv = do.call(rbind,
        lapply(simulations, .study06_simulation_summary)),
      operator_summaries.csv = operator_summary,
      operator_perturbation_metrics.csv =
        study06_operator_perturbation_metrics$block_metrics,
      full_csr_operator_action_metrics.csv =
        study06_operator_perturbation_metrics$full_csr_actions,
      synthetic_filter_validation.csv = study06_synthetic_filter_validation,
      block_level_operator_diagnostics.csv =
        study06_equivalence$block_metrics,
      fit_status.csv = x$status,
      prediction_metrics.csv = x$prediction,
      parameter_estimates.csv = x$parameter,
      parameter_recovery_summary.csv = x$recovery_summary,
      marker_effect_metrics.csv = x$marker,
      marker_effect_agreement.csv = x$marker_agreement,
      sbayesr_diagnostic_evidence.csv = x$sbayesr,
      paired_replicate_differences.csv = x$paired,
      paired_comparison_summary.csv = x$paired_summary,
      convergence_diagnostics.csv = x$diagnostics,
      convergence_validation_summary.csv =
        x$convergence_validation,
      computational_summary.csv = x$computational,
      method_output_availability.csv = availability)
    vapply(names(tables), function(name)
      .study06_write_csv(tables[[name]],
        file.path(study06_benchmark_output_dir, name)), "")
  }, format = "file")
)
