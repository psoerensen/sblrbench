source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
for (f in c("diagnostic_registry.R", "methods.R", "chain_extraction.R", "diagnostics.R", "pilot.R"))
  source(file.path("studies", "04_convergence", f), local = TRUE)
source(file.path("studies", "five_replicate_helpers.R"), local = TRUE)

.study04_validation_status <- function(successful, passing) {
  if (successful < 5L) return("indeterminate")
  if (passing == 5L) return("supported_in_all_replicates")
  if (passing >= 3L) return("supported_in_most_replicates")
  if (passing >= 1L) return("mixed")
  "not_supported"
}

.study04_validation_method_summary <- function(replicates, diagnostics) {
  keys <- interaction(replicates$architecture, replicates$method, drop = TRUE)
  do.call(rbind, lapply(split(replicates, keys), function(z) {
    d <- diagnostics[diagnostics$architecture == z$architecture[1L] &
      diagnostics$method == z$method[1L], , drop = FALSE]
    successful <- sum(z$status == "ok")
    passing <- sum(z$all_core_estimands_pass & z$status == "ok")
    limiting <- d$estimand[!d$overall_pass]
    limiting <- if (!length(limiting)) "none" else names(sort(table(limiting), decreasing = TRUE))[1L]
    data.frame(architecture = z$architecture[1L], method = z$method[1L],
      replicate_count = 5L, successful_replicates = successful,
      replicates_passing_all_core_estimands = passing, pass_proportion = passing / 5,
      maximum_rhat = if (nrow(d)) max(d$rhat, na.rm = TRUE) else NA_real_,
      minimum_bulk_ess = if (nrow(d)) min(d$ess_bulk, na.rm = TRUE) else NA_real_,
      minimum_tail_ess = if (nrow(d)) min(d$ess_tail, na.rm = TRUE) else NA_real_,
      maximum_relative_mcse = if (nrow(d)) max(d$relative_mcse, na.rm = TRUE) else NA_real_,
      most_frequent_limiting_estimand = limiting,
      recommendation_validation_status = .study04_validation_status(successful, passing),
      validation_rule = "indeterminate if fewer than five successful replicates; otherwise all=5, most=3-4, mixed=1-2, not_supported=0",
      stringsAsFactors = FALSE)
  }))
}

list(
  targets::tar_target(validation_config_file,
    file.path("studies", "04_convergence", "config.R"), format = "file"),
  targets::tar_target(validation_config, {
    x <- source(validation_config_file, local = TRUE)$value
    if (!identical(x$profile, "five_replicate_validation"))
      stop("Study 04 validation requires the five_replicate_validation profile.", call. = FALSE)
    .five_replicate_recommendations(x$profiles$five_replicate_validation$recommendation_source)
    x
  }),
  targets::tar_target(validation_paths, list(
    glist_path = Sys.getenv("SBLR_BENCH_GLIST", ""),
    data_dir = Sys.getenv("SBLR_BENCH_DATA_DIR",
      file.path("results", "local", "03_parameter_estimation", "data")),
    output_dir = Sys.getenv("SBLR_BENCH_LD_DIR",
      file.path("results", "local", "03_parameter_estimation", "genotype_setup")))),
  targets::tar_target(validation_example_files,
    .study01_example_files(validation_paths$data_dir, validation_config$example_data)),
  targets::tar_target(validation_base_glist,
    .study01_load_glist(validation_paths, validation_example_files)),
  targets::tar_target(validation_markers,
    .study01_run_qc(validation_base_glist, validation_config)),
  targets::tar_target(validation_ids,
    .study01_selected_ids(validation_base_glist, validation_config$sample_limit)),
  targets::tar_target(validation_working_glist, .study01_set_rsids_ld(
    validation_base_glist, validation_config$chr, validation_markers$marker_ids)),
  targets::tar_target(validation_genotypes, .study01_extract_genotypes(
    validation_working_glist, validation_config$chr, validation_ids,
    validation_markers$marker_ids)),
  targets::tar_target(validation_sparse_ld_glist, .study01_make_sparse_ld(
    validation_working_glist, validation_markers, validation_config,
    validation_paths$output_dir)),
  targets::tar_target(validation_sim_specs,
    .study04_sim_specs(validation_config), iteration = "list"),
  targets::tar_target(validation_sim_spec, validation_sim_specs,
    pattern = map(validation_sim_specs), iteration = "list"),
  targets::tar_target(validation_sim_bundle, {
    sim <- sblrbench:::simulate_prediction_architecture(list(
      scenario=validation_sim_spec$architecture,
      replicate=validation_sim_spec$replicate,
      simulation_seed=validation_sim_spec$simulation_seed),
      validation_genotypes,validation_config$parameter_spec)
    oracle <- sblrbench::check_oracle_genetic_values(sim,
      tolerance = validation_config$oracle_tolerance)
    if (!oracle$ok) stop("Study 04 validation oracle failed.", call. = FALSE)
    stats <- sblr::make_summary_stats(Glist = validation_sparse_ld_glist,
      y = sim$truth$phenotypes, chr = validation_config$chr,
      rows = seq_len(nrow(sim$truth$phenotypes)), scale = TRUE, nthreads = 1L)
    list(simulation = sim, stats = stats,
      truth = sblrbench:::parameter_estimand_truth(sim,
        validation_config$parameter_spec))
  }, pattern = map(validation_sim_spec), iteration = "list"),
  targets::tar_target(validation_specs,
    .study04_specs(validation_config), iteration = "list"),
  targets::tar_target(validation_spec, validation_specs,
    pattern = map(validation_specs), iteration = "list"),
  targets::tar_target(validation_method_run, {
    hit <- which(vapply(validation_sim_bundle, function(x)
      x$simulation$scenario$architecture == validation_spec$architecture &&
        x$simulation$scenario$replicate == validation_spec$replicate, logical(1)))
    if (length(hit) != 1L) stop("No unique paired simulation bundle.", call. = FALSE)
    bundle <- validation_sim_bundle[[hit]]
    run <- .study04_fit(validation_spec, bundle$simulation, bundle$stats,
      validation_sparse_ld_glist, validation_config)
    draws <- if (run$status == "ok") tryCatch(.study04_extract_chain_draws(
      run$fit, validation_spec$architecture, validation_spec$method,
      validation_spec$replicate), error = function(e) {
        run$status <<- "failed"; run$reason <<- conditionMessage(e); data.frame()
      }) else data.frame()
    agreement <- .study04_marker_agreement(run, validation_spec$architecture)
    agreement$replicate <- validation_spec$replicate
    list(run = run, draws = draws, agreement = agreement, spec = validation_spec)
  }, pattern = map(validation_spec), iteration = "list"),
  targets::tar_target(validation_draws_branch,
    validation_method_run$draws,
    pattern = map(validation_method_run), iteration = "list"),
  targets::tar_target(validation_agreement_branch,
    validation_method_run$agreement,
    pattern = map(validation_method_run), iteration = "list"),
  targets::tar_target(validation_run_metadata_branch, {
    x <- validation_method_run
    list(spec = x$spec, status = x$run$status, reason = x$run$reason,
      runtime = x$run$runtime)
  }, pattern = map(validation_method_run), iteration = "list"),
  targets::tar_target(validation_scalar_chain_draws,
    do.call(rbind, validation_draws_branch)),
  targets::tar_target(validation_convergence_diagnostics, {
    groups <- split(validation_scalar_chain_draws, interaction(
      validation_scalar_chain_draws$architecture,
      validation_scalar_chain_draws$replicate,
      validation_scalar_chain_draws$method, drop = TRUE))
    do.call(rbind, lapply(groups, function(z) .study04_diagnostics(z, 0L,
      length(unique(z$raw_iteration)), .study04_registry(), validation_config$thresholds)))
  }),
  targets::tar_target(validation_replicate_diagnostic_summary, {
    expected <- do.call(rbind, lapply(validation_specs, as.data.frame))
    rows <- lapply(seq_len(nrow(expected)), function(i) {
      d <- validation_convergence_diagnostics[
        validation_convergence_diagnostics$architecture == expected$architecture[i] &
        validation_convergence_diagnostics$replicate == expected$replicate[i] &
        validation_convergence_diagnostics$method == expected$method[i], , drop = FALSE]
      run <- validation_run_metadata_branch[[i]]
      limiting <- d$estimand[!d$overall_pass]
      data.frame(expected[i, ], status = run$status,
        all_core_estimands_pass = nrow(d) == 4L && all(d$overall_pass),
        limiting_estimands = if (length(limiting)) paste(limiting, collapse = ";") else "none",
        maximum_rhat = if (nrow(d)) max(d$rhat) else NA_real_,
        minimum_bulk_ess = if (nrow(d)) min(d$ess_bulk) else NA_real_,
        minimum_tail_ess = if (nrow(d)) min(d$ess_tail) else NA_real_,
        maximum_relative_mcse = if (nrow(d)) max(d$relative_mcse) else NA_real_,
        reason = run$reason, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }),
  targets::tar_target(validation_method_summary,
    .study04_validation_method_summary(validation_replicate_diagnostic_summary,
      validation_convergence_diagnostics)),
  targets::tar_target(validation_chain_status, do.call(rbind,
    lapply(validation_run_metadata_branch, function(x) do.call(rbind, lapply(seq_len(4L), function(ch)
      data.frame(architecture = x$spec$architecture, replicate = x$spec$replicate,
        method = x$spec$method, chain = ch, status = x$status,
        error_class = if (x$status == "ok") "" else "fit_or_extraction_error",
        error_message = x$reason, stringsAsFactors = FALSE)))))),
  targets::tar_target(validation_computational_summary, do.call(rbind,
    lapply(validation_run_metadata_branch, function(x) {
      p <- .five_replicate_mcmc(x$spec$method)
      data.frame(architecture = x$spec$architecture, replicate = x$spec$replicate,
        method = x$spec$method, status = x$status, runtime = x$runtime,
        nit = p$nit, nburn = p$nburn, nthin = p$nthin, nchains = p$nchains,
        ncores = p$ncores, error_message = x$reason,
        OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
        OMP_THREAD_LIMIT = Sys.getenv("OMP_THREAD_LIMIT"),
        OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS"),
        MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS"), stringsAsFactors = FALSE)
    }))),
  targets::tar_target(validation_chain_summaries, {
    ids <- c("effect_variance", "genetic_variance", "residual_variance", "heritability")
    groups <- split(validation_scalar_chain_draws, interaction(
      validation_scalar_chain_draws$architecture, validation_scalar_chain_draws$replicate,
      validation_scalar_chain_draws$method, validation_scalar_chain_draws$chain, drop = TRUE))
    do.call(rbind, lapply(groups, function(z) do.call(rbind, lapply(ids, function(id)
      data.frame(architecture = z$architecture[1], replicate = z$replicate[1],
        method = z$method[1], chain = z$chain[1], estimand = id,
        draw_count = nrow(z), mean = mean(z[[id]]), sd = stats::sd(z[[id]]),
        minimum = min(z[[id]]), maximum = max(z[[id]]), stringsAsFactors = FALSE)))))
  }),
  targets::tar_target(validation_chain_agreement,
    do.call(rbind, validation_agreement_branch)),
  targets::tar_target(validation_seed_registry, do.call(rbind,
    lapply(validation_specs, function(x) {
      fit_seed <- .study04_chain_seeds(x$architecture, x$method, validation_config,
        x$replicate)
      sim <- validation_sim_specs[[which(vapply(validation_sim_specs, function(s)
        s$architecture == x$architecture && s$replicate == x$replicate, logical(1)))]]
      data.frame(architecture = x$architecture, replicate = x$replicate,
        method = x$method, data_selection_seed = 3301L,
        architecture_seed = validation_config$simulation$base_seed +
          match(x$architecture, names(validation_config$simulation$architectures)) * 1000L,
        simulation_seed = sim$simulation_seed, fit_seed = fit_seed[1L],
        chain = seq_len(4L), chain_seed = fit_seed, stringsAsFactors = FALSE)
    }))),
  targets::tar_target(validation_output_dir, Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
    file.path("results", "local", "five_replicate_overnight", "04_convergence"))),
  targets::tar_target(validation_draws_file, .study04_write(validation_scalar_chain_draws,
    file.path(validation_output_dir, "scalar_chain_draws.csv")), format = "file"),
  targets::tar_target(validation_diagnostics_file, .study04_write(validation_convergence_diagnostics,
    file.path(validation_output_dir, "convergence_diagnostics.csv")), format = "file"),
  targets::tar_target(validation_replicate_summary_file, .study04_write(validation_replicate_diagnostic_summary,
    file.path(validation_output_dir, "replicate_diagnostic_summary.csv")), format = "file"),
  targets::tar_target(validation_method_summary_file, .study04_write(validation_method_summary,
    file.path(validation_output_dir, "method_validation_summary.csv")), format = "file"),
  targets::tar_target(validation_chain_summary_file, .study04_write(validation_chain_summaries,
    file.path(validation_output_dir, "chain_summaries.csv")), format = "file"),
  targets::tar_target(validation_agreement_file, .study04_write(validation_chain_agreement,
    file.path(validation_output_dir, "chain_agreement.csv")), format = "file"),
  targets::tar_target(validation_computational_file, .study04_write(validation_computational_summary,
    file.path(validation_output_dir, "computational_summary.csv")), format = "file"),
  targets::tar_target(validation_status_file, .study04_write(validation_chain_status,
    file.path(validation_output_dir, "chain_status.csv")), format = "file"),
  targets::tar_target(validation_seed_file, .study04_write(validation_seed_registry,
    file.path(validation_output_dir, "seed_registry.csv")), format = "file"),
  targets::tar_target(validation_manifest_file, {
    sblr <- .five_replicate_sblr_provenance()
    complete <- nrow(validation_chain_status) == 80L &&
      all(validation_chain_status$status == "ok")
    path <- file.path(validation_output_dir, "benchmark_manifest.json")
    jsonlite::write_json(list(study_id = "04_convergence_validation",
      task = "single_trait_multichain_convergence_validation",
      benchmark_scope = "five_replicate_fixed_setting_validation",
      benchmark_status = if (complete) "complete" else "incomplete",
      validation_scope = "five simulations; not universal validation",
      replicate_count = 5L, expected_method_fit_count = 20L,
      successful_method_fit_count = sum(validation_computational_summary$status == "ok"),
      failed_method_fit_count = sum(validation_computational_summary$status != "ok"),
      expected_chain_count = 80L,
      successful_chain_count = sum(validation_chain_status$status == "ok"),
      failed_chain_count = sum(validation_chain_status$status != "ok"),
      active_methods = validation_config$methods,
      architecture_grid = validation_config$matched_grid,
      mcmc = .five_replicate_recommendations(), chain_count = 4L,
      seeds = validation_seed_registry, package_versions = list(
        sblr = sblr$version, sblrbench = as.character(packageVersion("sblrbench")),
        posterior = as.character(packageVersion("posterior"))),
      repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
      sblr_commit = sblr$commit, sblr_source_status = sblr$source_status,
      qgdata_commit = validation_config$example_data$commit,
      source_status = "uncommitted five-replicate validation source; source inventory retained in capsule",
      completion_counts = list(method_fits = nrow(validation_computational_summary),
        chains = nrow(validation_chain_status)),
      failures = validation_chain_status[validation_chain_status$status != "ok", ],
      provenance = validation_config$example_data,
      validation_status = validation_method_summary,
      diagnostic_thresholds = validation_config$thresholds,
      thread_settings = .five_replicate_thread_settings()), path, pretty = TRUE,
      auto_unbox = TRUE)
    path
  }, format = "file")
)
