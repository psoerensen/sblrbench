source(file.path("studies", "01_finemapping", "setup_example_data.R"), local = TRUE)
for (f in c("diagnostic_registry.R", "methods.R", "chain_extraction.R", "diagnostics.R", "recommendations.R", "pilot.R")) source(file.path("studies", "04_convergence", f), local = TRUE)

list(
  targets::tar_target(convergence_config_file, file.path("studies", "04_convergence", "config.R"), format="file"),
  targets::tar_target(convergence_config, source(convergence_config_file, local=TRUE)$value),
  targets::tar_target(convergence_paths, list(
    glist_path=Sys.getenv("SBLR_BENCH_GLIST",""),
    data_dir=Sys.getenv("SBLR_BENCH_DATA_DIR",file.path("results","local","03_parameter_estimation","data")),
    output_dir=Sys.getenv("SBLR_BENCH_LD_DIR",file.path("results","local","03_parameter_estimation","genotype_setup")))),
  targets::tar_target(convergence_example_files, .study01_example_files(convergence_paths$data_dir, convergence_config$example_data)),
  targets::tar_target(convergence_base_glist, .study01_load_glist(convergence_paths, convergence_example_files)),
  targets::tar_target(convergence_markers, .study01_run_qc(convergence_base_glist, convergence_config)),
  targets::tar_target(convergence_ids, .study01_selected_ids(convergence_base_glist, convergence_config$sample_limit)),
  targets::tar_target(convergence_working_glist, .study01_set_rsids_ld(convergence_base_glist, convergence_config$chr, convergence_markers$marker_ids)),
  targets::tar_target(convergence_genotypes, .study01_extract_genotypes(convergence_working_glist, convergence_config$chr, convergence_ids, convergence_markers$marker_ids)),
  targets::tar_target(convergence_sparse_ld_glist, .study01_make_sparse_ld(convergence_working_glist, convergence_markers, convergence_config, convergence_paths$output_dir)),
  targets::tar_target(convergence_sim_specs, .study04_sim_specs(convergence_config), iteration="list"),
  targets::tar_target(convergence_sim_spec, convergence_sim_specs, pattern=map(convergence_sim_specs), iteration="list"),
  targets::tar_target(convergence_sim_bundle, { sim <- sblrbench:::simulate_prediction_architecture(list(scenario=convergence_sim_spec$architecture,replicate=convergence_sim_spec$replicate,simulation_seed=convergence_sim_spec$simulation_seed),convergence_genotypes,convergence_config$parameter_spec); oracle <- sblrbench::check_oracle_genetic_values(sim, tolerance=convergence_config$oracle_tolerance); if(!oracle$ok) stop("oracle failed"); stats <- sblrbench:::parameter_summary_stats(sim,convergence_sparse_ld_glist,convergence_config$parameter_spec); list(simulation=sim,stats=stats,truth=sblrbench:::parameter_estimand_truth(sim,convergence_config$parameter_spec),oracle=oracle)}, pattern=map(convergence_sim_spec), iteration="list"),
  targets::tar_target(convergence_specs, .study04_specs(convergence_config), iteration="list"),
  targets::tar_target(convergence_spec, convergence_specs, pattern=map(convergence_specs), iteration="list"),
  targets::tar_target(convergence_method_run, { bundle <- convergence_sim_bundle[[match(convergence_spec$architecture,vapply(convergence_sim_bundle,function(x)x$simulation$scenario$architecture,character(1)))]]; run <- .study04_fit(convergence_spec,bundle$simulation,bundle$stats,convergence_sparse_ld_glist,convergence_config); draws <- if(run$status=="ok") .study04_extract_chain_draws(run$fit,convergence_spec$architecture,convergence_spec$method) else data.frame(); list(run=run,draws=draws,agreement=.study04_marker_agreement(run,convergence_spec$architecture),spec=convergence_spec)}, pattern=map(convergence_spec), iteration="list"),
  targets::tar_target(convergence_scalar_chain_draws, do.call(rbind,lapply(convergence_method_run,`[[`,"draws"))),
  targets::tar_target(convergence_diagnostics, do.call(rbind,lapply(split(convergence_scalar_chain_draws,interaction(convergence_scalar_chain_draws$architecture,convergence_scalar_chain_draws$method,drop=TRUE)),.study04_diagnostic_grid,config=convergence_config))),
  targets::tar_target(convergence_burnin_stability, do.call(rbind,lapply(split(convergence_scalar_chain_draws,interaction(convergence_scalar_chain_draws$architecture,convergence_scalar_chain_draws$method,drop=TRUE)),.study04_burnin_stability,diagnostics=convergence_diagnostics,config=convergence_config))),
  targets::tar_target(convergence_recommendations, .study04_recommend(convergence_diagnostics,convergence_burnin_stability,convergence_config)),
  targets::tar_target(convergence_chain_agreement, do.call(rbind,lapply(convergence_method_run,`[[`,"agreement"))),
  targets::tar_target(convergence_computational, do.call(rbind,lapply(convergence_method_run,function(x) do.call(rbind,lapply(1:4,function(ch)data.frame(architecture=x$spec$architecture,method=x$spec$method,chain=ch,status=x$run$status,runtime=NA_real_,method_total_runtime=x$run$runtime,error_class=if(x$run$status=="ok")"" else "fit_error",error_message=x$run$reason,raw_draw_count=if(x$run$status=="ok")3000L else NA,nit_argument=3000L,nburn_argument=0L,nthin=1L,nchains_argument=4L,ncores_argument=4L,fit_seed=x$run$chain_seeds[ch],OMP_NUM_THREADS=Sys.getenv("OMP_NUM_THREADS",NA),OMP_THREAD_LIMIT=Sys.getenv("OMP_THREAD_LIMIT",NA),OPENBLAS_NUM_THREADS=Sys.getenv("OPENBLAS_NUM_THREADS",NA),MKL_NUM_THREADS=Sys.getenv("MKL_NUM_THREADS",NA))))))),
  targets::tar_target(convergence_chain_status, convergence_computational[,c("architecture","method","chain","status","error_class","error_message")]),
  targets::tar_target(convergence_chain_summaries, aggregate(convergence_scalar_chain_draws[5:8], convergence_scalar_chain_draws[1:3], function(x) c(mean=mean(x),sd=sd(x),minimum=min(x),maximum=max(x)))),
  targets::tar_target(convergence_checkpoint_summary, aggregate(overall_pass ~ architecture + method + burnin_candidate + retained_draw_candidate, convergence_diagnostics, function(x) all(x))),
  targets::tar_target(convergence_seed_registry, do.call(rbind,lapply(convergence_specs,function(x)data.frame(architecture=x$architecture,method=x$method,chain=1:4,simulation_seed=.study04_sim_specs(convergence_config)[[match(x$architecture,names(convergence_config$simulation$architectures))]]$simulation_seed,fit_seed=.study04_chain_seeds(x$architecture,x$method,convergence_config))))),
  targets::tar_target(convergence_simulation_truth, do.call(rbind,lapply(convergence_sim_bundle,`[[`,"truth"))),
  targets::tar_target(convergence_registry, .study04_registry()),
  targets::tar_target(convergence_thresholds, .study04_thresholds(convergence_config)),
  targets::tar_target(convergence_output_dir, Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
    file.path("results", "local", "04_convergence"))),
  targets::tar_target(convergence_draws_file,.study04_write(convergence_scalar_chain_draws,file.path(convergence_output_dir,"scalar_chain_draws.csv")),format="file"),
  targets::tar_target(convergence_diagnostics_file,.study04_write(convergence_diagnostics,file.path(convergence_output_dir,"convergence_diagnostics.csv")),format="file"),
  targets::tar_target(convergence_stability_file,.study04_write(convergence_burnin_stability,file.path(convergence_output_dir,"burnin_stability.csv")),format="file"),
  targets::tar_target(convergence_recommendations_file,.study04_write(convergence_recommendations,file.path(convergence_output_dir,"method_recommendations.csv")),format="file"),
  targets::tar_target(convergence_computational_file,.study04_write(convergence_computational,file.path(convergence_output_dir,"computational_summary.csv")),format="file")
  ,targets::tar_target(convergence_checkpoint_file,.study04_write(convergence_checkpoint_summary,file.path(convergence_output_dir,"checkpoint_summary.csv")),format="file")
  ,targets::tar_target(convergence_chain_summary_file,.study04_write(convergence_chain_summaries,file.path(convergence_output_dir,"chain_summaries.csv")),format="file")
  ,targets::tar_target(convergence_agreement_file,.study04_write(convergence_chain_agreement,file.path(convergence_output_dir,"chain_agreement.csv")),format="file")
  ,targets::tar_target(convergence_status_file,.study04_write(convergence_chain_status,file.path(convergence_output_dir,"chain_status.csv")),format="file")
  ,targets::tar_target(convergence_seed_file,.study04_write(convergence_seed_registry,file.path(convergence_output_dir,"seed_registry.csv")),format="file")
  ,targets::tar_target(convergence_truth_file,.study04_write(convergence_simulation_truth,file.path(convergence_output_dir,"simulation_truth.csv")),format="file")
  ,targets::tar_target(convergence_registry_file,.study04_write(convergence_registry,file.path(convergence_output_dir,"diagnostic_registry.csv")),format="file")
  ,targets::tar_target(convergence_threshold_file,.study04_write(convergence_thresholds,file.path(convergence_output_dir,"diagnostic_thresholds.csv")),format="file")
  ,targets::tar_target(convergence_manifest_file, {
    path <- file.path(convergence_output_dir,"benchmark_manifest.json")
    jsonlite::write_json(list(study_id=convergence_config$study,task=convergence_config$task,
      benchmark_scope="one_replicate_matched_architecture_development",
      benchmark_status=if(nrow(convergence_chain_status)==16L && all(convergence_chain_status$status=="ok")) "complete" else "incomplete",
      development_settings=TRUE,created_at=format(Sys.time(),tz="UTC",usetz=TRUE),
      repository_commit=system2("git",c("rev-parse","HEAD"),stdout=TRUE),
      source_tree_clean=!length(system2("git",c("status","--porcelain"),stdout=TRUE)),
      source_state="uncommitted Study 04 snapshot; source_files.csv records content hashes",
      sblr_commit=system2("git",c("-C","../sblr","rev-parse","HEAD"),stdout=TRUE),
      qgdata_commit=convergence_config$example_data$commit,active_methods=convergence_config$methods,
      architectures=names(convergence_config$simulation$architectures),matched_method_architecture_grid=convergence_config$matched_grid,
      replicate_count=1L,expected_method_fit_count=4L,successful_method_fit_count=sum(vapply(convergence_method_run,function(x)x$run$status=="ok",logical(1))),
      expected_chain_count=16L,successful_chain_count=sum(convergence_chain_status$status=="ok"),failed_chain_count=sum(convergence_chain_status$status!="ok"),
      analysis_sample_count=length(convergence_ids),canonical_marker_count=length(convergence_markers$marker_ids),causal_marker_count=convergence_config$simulation$n_causal,target_h2=convergence_config$simulation$h2,
      chain_execution_strategy="native multichain with retained chain identities",native_multichain_used=TRUE,
      chain_retention_contract="convergence_traces: iteration x chain x quantity; fit$chains retained",
      iteration_semantics="nit post-warmup draws; benchmark nburn=0 retains all raw history; candidate burn-in applied downstream once",
      mcmc_maximum_run=convergence_config$profiles$development,burnin_candidates=convergence_config$burnin_candidates,
      retained_draw_candidates=convergence_config$retained_draw_candidates,diagnostic_thresholds=convergence_config$thresholds,
      seeds=convergence_seed_registry,core_estimands=convergence_registry$estimand[convergence_registry$classification=="core"],
      extended_estimands=convergence_registry$estimand[convergence_registry$classification=="extended"],recommendations=convergence_recommendations,
      data_provenance=convergence_config$example_data),path,pretty=TRUE,auto_unbox=TRUE)
    path
  },format="file")
)
