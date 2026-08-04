# Common execution path for completed Studies 01--05.

benchmark_output_paths <- function(output_dir) {
  .benchmark_scalar_string(output_dir, "output_dir")
  list(root = output_dir, checkpoints = file.path(output_dir, "checkpoints"),
    qualification = file.path(output_dir, "qualification"),
    tables = file.path(output_dir, "tables"),
    figures = file.path(output_dir, "figures"),
    manifest = file.path(output_dir, "manifest.json"),
    session_info = file.path(output_dir, "session_info.txt"))
}

parse_benchmark_cli_arguments <- function(args) {
  allowed <- c("--study", "--profile", "--output-dir", "--resume",
    "--validate-only", "--mode")
  if (!is.character(args) || length(args) %% 2L != 0L)
    stop("Every command-line option requires a value.", call. = FALSE)
  options <- list()
  if (length(args)) for (i in seq(1L, length(args), by = 2L)) {
    option <- args[[i]]
    if (!option %in% allowed)
      stop("Unknown command-line option: ", option, call. = FALSE)
    if (!is.null(options[[option]]))
      stop("Command-line option supplied more than once: ", option,
        call. = FALSE)
    options[[option]] <- args[[i + 1L]]
  }
  missing <- setdiff(c("--study", "--profile", "--output-dir"),
    names(options))
  if (length(missing))
    stop("Missing required command-line options: ",
      paste(missing, collapse = ", "), ".", call. = FALSE)
  supported_studies <- c("01_finemapping", "02_prediction",
    "03_parameter_estimation", "04_convergence", "05_ld_operator",
    "06_annotation_models")
  if (identical(options[["--study"]], "06_ld_operator"))
    stop("Study ID `06_ld_operator` is retired; use `05_ld_operator`.",
      call. = FALSE)
  if (identical(options[["--study"]], "05_annotation_models"))
    stop("Study ID `05_annotation_models` is retired; use ",
      "`06_annotation_models` (in development).", call. = FALSE)
  if (identical(options[["--study"]], "07_mtblr_validation"))
    stop("Study ID `07_mtblr_validation` is retired; use `07_mt_validation` ",
      "(in development).", call. = FALSE)
  if (identical(options[["--study"]], "07_mt_validation"))
    stop("Study `07_mt_validation` is in development and is not yet ",
      "supported by run_benchmark().", call. = FALSE)
  if (!options[["--study"]] %in% supported_studies)
    stop("Unsupported --study; choose a completed Study 01--05 ID or `06_annotation_models`.",
      call. = FALSE)
  if (!options[["--profile"]] %in% c("workshop", "benchmark"))
    stop("--profile must be `workshop` or `benchmark`.", call. = FALSE)
  parse_boolean <- function(value, option) {
    normalized <- tolower(value)
    if (!normalized %in% c("true", "false"))
      stop(option, " must be `true` or `false`.", call. = FALSE)
    identical(normalized, "true")
  }
  mode <- if (is.null(options[["--mode"]])) NULL else options[["--mode"]]
  if (!is.null(mode) && !mode %in% c("validate_only", "qualification", "final"))
    stop("--mode must be `validate_only`, `qualification`, or `final`.",
      call. = FALSE)
  list(study = options[["--study"]], profile = options[["--profile"]],
    output_dir = options[["--output-dir"]],
    resume = parse_boolean(if (is.null(options[["--resume"]])) "true" else
      options[["--resume"]], "--resume"),
    validate_only = parse_boolean(
      if (is.null(options[["--validate-only"]])) "false" else
        options[["--validate-only"]], "--validate-only"), mode = mode)
}

.run_finemapping <- function(spec,profile,coordinates,methods,paths,resume) {
  data <- prepare_finemapping_data(spec,paths$root)
  bundles <- prepare_finemapping_simulations(spec,profile,data)
  names(bundles) <- vapply(bundles,function(x) .prediction_bundle_key(
    x$coordinate$scenario,x$coordinate$replicate),character(1))
  dispatch <- getOption("sblrbench.finemapping_fit_dispatch",
    fit_finemapping_method)
  status_rows <- runtime_rows <- marker_rows <- set_rows <- list()
  for(i in seq_len(nrow(coordinates))) {
    coordinate <- coordinates[i,,drop=FALSE]
    bundle <- bundles[[.prediction_bundle_key(coordinate$scenario,
      coordinate$replicate)]]
    method <- methods[[coordinate$method]]
    controls <- benchmark_method_controls(spec,coordinate$method,profile,
      coordinate$fit_seed,coordinate$chain_seeds[[1L]])
    identity <- benchmark_semantic_checkpoint_identity(
      diagnostic_id="study01-finemapping-fit",scientific_inputs=list(
        study=spec$study,scenario=coordinate$scenario,
        replicate=coordinate$replicate,method=coordinate$method,
        sample_ids=data$sample_ids,marker_ids=data$markers$marker_ids,
        causal_markers=bundle$selection$marker_ids,
        simulation_seed=coordinate$simulation_seed,controls=controls,
        sblr_sha=spec$packages$sblr$sha,qgdata_sha=spec$packages$qgdata$sha))
    input_hash <- benchmark_semantic_checkpoint_hash(identity)
    checkpoint <- .prediction_checkpoint_path(paths,as.list(coordinate))
    result <- NULL; reused <- FALSE; reason <- ""
    if(isTRUE(resume) && file.exists(checkpoint)) {
      loaded <- benchmark_load_semantic_checkpoint(checkpoint,input_hash,
        validator=function(x) inherits(x$result,"sblrbench_result"))
      result <- loaded$value$result; reused <- TRUE
    } else {
      result <- tryCatch(dispatch(method=method,controls=controls,
        simulation=bundle$simulation,stats=bundle$stats,glist=data$ld_glist),
        error=function(e) {reason <<- conditionMessage(e); NULL})
      if(!is.null(result)) benchmark_atomic_save_rds(list(
        checkpoint_schema="sblrbench-semantic-v2",identity_payload=identity,
        semantic_hash=input_hash,result=result),checkpoint,compress=FALSE,
        temporary_prefix=".finemapping-fit-")
    }
    ok <- !is.null(result)
    status_rows[[i]] <- data.frame(study=spec$study,scenario=coordinate$scenario,
      replicate=coordinate$replicate,method=coordinate$method,
      status=if(ok)"ok" else "failed",reason=reason,reused=reused,
      input_hash=input_hash,stringsAsFactors=FALSE)
    runtime_rows[[i]] <- data.frame(study=spec$study,scenario=coordinate$scenario,
      replicate=coordinate$replicate,method=coordinate$method,
      elapsed_seconds=if(ok)extract_runtime(result) else NA_real_,
      status=if(ok)"ok" else "failed",reason=reason,reused=reused,
      stringsAsFactors=FALSE)
    if(ok) {
      marker_rows[[i]] <- finemapping_marker_table(result,bundle$simulation,
        coordinate$method)
      cs <- spec$locus_design
      native <- sblr::make_credible_sets(fit=result$native_fit,
        Glist=data$ld_glist,trait=1L,coverage=cs$credible_set_target,
        min_r2=cs$min_r2,pip_cutoff=cs$pip_cutoff,
        locus_pip_cutoff=cs$locus_pip_cutoff,
        max_locus_distance=cs$max_locus_distance)
      member_ids <- unique(finemapping_credible_set_members(native)$marker)
      ids <- intersect(unique(c(bundle$selection$marker_ids,member_ids)),
        colnames(data$scaled))
      ld <- if(length(ids)) stats::cor(data$scaled[,ids,drop=FALSE]) else
        matrix(numeric(),0L,0L)
      index <- match(data$markers$marker_ids,
        data$working_glist$rsids[[spec$data$chromosome]])
      positions <- stats::setNames(
        data$working_glist$pos[[spec$data$chromosome]][index],
        data$markers$marker_ids)
      set_rows[[i]] <- evaluate_finemapping_credible_sets(native,
        bundle$simulation,positions,ld,coordinate$method,cs$min_r2,
        cs$max_locus_distance)
    }
  }
  status <- .benchmark_bind_rows(status_rows)
  runtime <- .benchmark_bind_rows(runtime_rows)
  marker_results <- .benchmark_bind_rows(marker_rows)
  credible_sets <- .benchmark_bind_rows(set_rows)
  metrics <- if(is.null(marker_results)) NULL else
    summarise_finemapping_metrics(marker_results,credible_sets)
  truth <- .benchmark_bind_rows(lapply(bundles,.prediction_truth_table))
  loci <- .benchmark_bind_rows(lapply(bundles,function(x)
    transform(x$loci,scenario=x$coordinate$scenario,
      replicate=x$coordinate$replicate)))
  oracle <- .benchmark_oracle_table(bundles,spec$study)
  files <- list(
    fit_status=.benchmark_write_csv(status,file.path(paths$tables,"fit_status.csv")),
    simulation_truth=.benchmark_write_csv(truth,file.path(paths$tables,"simulation_truth.csv")),
    loci=.benchmark_write_csv(loci,file.path(paths$tables,"loci.csv")),
    marker_results=.benchmark_write_csv(marker_results,file.path(paths$tables,"marker_results.csv")),
    credible_sets=.benchmark_write_csv(credible_sets,file.path(paths$tables,"credible_sets.csv")),
    finemapping_metrics=.benchmark_write_csv(metrics,file.path(paths$tables,"finemapping_metrics.csv")),
    runtime=.benchmark_write_csv(runtime,file.path(paths$tables,"runtime.csv")),
    simulation_oracle=.benchmark_write_csv(oracle,file.path(paths$tables,"simulation_oracle.csv")))
  manifest <- .prediction_manifest(spec,profile,paths,coordinates,FALSE,status,data)
  .write_prediction_manifest(manifest,paths$manifest)
  writeLines(benchmark_session_information(),paths$session_info)
  list(spec=spec,paths=c(paths,files),status=status,truth=truth,loci=loci,
    marker_results=marker_results,credible_sets=credible_sets,metrics=metrics,
    convergence=NULL,runtime=runtime,oracle=oracle)
}

.benchmark_create_output_dirs <- function(paths) {
  for (path in paths[c("root", "checkpoints", "qualification", "tables",
      "figures")])
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(paths)
}

.benchmark_write_csv <- function(x, path) {
  if (!is.data.frame(x) || !nrow(x)) return(NULL)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  path
}

.benchmark_bind_rows <- function(x) {
  x <- Filter(function(value) is.data.frame(value) && nrow(value), x)
  if (!length(x)) NULL else {
    out <- do.call(rbind, x)
    rownames(out) <- NULL
    out
  }
}

.run_ld_operator <- function(spec, profile, coordinates, paths, resume) {
  runner <- getOption("sblrbench.ld_operator_runner")
  if (!is.function(runner))
    stop("Study 05 execution requires its study-specific operator runner. ",
      "Source `studies/05_ld_operator/operator-design.R` before calling ",
      "run_benchmark() outside validate-only mode.", call. = FALSE)
  result <- runner(spec = spec, profile = profile, coordinates = coordinates,
    paths = paths, resume = resume)
  required <- c("status", "operator_summary", "operator_comparisons",
    "eigenvalue_summary", "convergence", "recovery_metrics", "runtime",
    "sbayesr_fit_status", "sbayesr_scheduler", "sbayesr_variance",
    "sbayesr_conditionals", "sbayesr_quadratics", "sbayesr_recovery",
    "sbayesr_eigenvalues")
  missing <- setdiff(required, names(result))
  if (length(missing))
    stop("The Study 05 operator runner omitted required result components: ",
      paste(missing, collapse = ", "), ".", call. = FALSE)
  result$spec <- spec
  result$paths <- c(paths, result$paths)
  result
}

.prediction_validation_status <- function(coordinates, study = "02_prediction") {
  out <- data.frame(study = study, scenario = coordinates$scenario,
    replicate = coordinates$replicate, method = coordinates$method,
    status = "not_run_validate_only",
    reason = "Specification, profile, coordinates, seeds, methods, and package provenance validated; fitting disabled.",
    stringsAsFactors = FALSE)
  if ("stage" %in% names(coordinates))
    out <- cbind(out["study"], stage = coordinates$stage, out[-1L])
  out
}

.prediction_manifest <- function(spec, profile, paths, coordinates,
                                 validate_only, status = NULL,
                                 data = NULL) {
  package <- benchmark_package_provenance("sblr")
  list(schema_version = 1L, study = spec$study, task = spec$task,
    profile = profile, validate_only = isTRUE(validate_only),
    coordinate_count = nrow(coordinates),
    scenarios = if (is.null(spec$scenarios)) unique(coordinates$scenario) else
      names(spec$scenarios),
    methods = if (is.null(spec$methods)) unique(coordinates$method) else
      names(spec$methods),
    replicate_count = resolve_benchmark_profile(spec, profile)$replicate_count,
    split = if(is.null(spec$split)) NULL else spec$split,
    marker_policy = if (is.null(spec$markers)) NULL else spec$markers,
    scientific_controls = spec$controls, seeds = spec$seeds,
    qgdata = if (is.null(spec$data)) spec$packages$qgdata else
      spec$data$example_data,
    packages = list(sblr = package, expected = spec$packages),
    data_summary = if (is.null(data)) NULL else list(
      sample_count = length(data$sample_ids),
      training_sample_count = if(is.null(data$split)) NULL else
        length(data$split$train_ids),
      test_sample_count = if(is.null(data$split)) NULL else
        length(data$split$test_ids),
      marker_count = length(data$markers$marker_ids),
      sample_order_hash = benchmark_hash_object(data$sample_ids),
      marker_order_hash = benchmark_hash_object(data$markers$marker_ids)),
    completion = if (is.null(status)) NULL else as.list(table(status$status)),
    output_layout = paths,
    source_commit = benchmark_git_sha(".", warn = FALSE),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

.write_prediction_manifest <- function(manifest, path) {
  jsonlite::write_json(manifest, path, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null")
  path
}

.prediction_bundle_key <- function(scenario, replicate) {
  paste(scenario, as.integer(replicate), sep = "::")
}

.benchmark_oracle_table <- function(bundles, study) {
  .benchmark_bind_rows(lapply(bundles, function(bundle) data.frame(
    study = study, scenario = bundle$coordinate$scenario,
    replicate = bundle$coordinate$replicate,
    status = if (isTRUE(bundle$oracle$ok)) "passed" else "failed",
    max_abs_error = bundle$oracle$max_abs_error,
    tolerance = bundle$oracle$tolerance, stringsAsFactors = FALSE)))
}

.run_parameter_estimation <- function(spec,profile,resolved,coordinates,
                                      methods,paths,resume) {
  data <- prepare_parameter_estimation_data(spec,paths$root)
  bundles <- prepare_parameter_estimation_simulations(spec,profile,data)
  keys <- vapply(bundles,function(x) .prediction_bundle_key(
    x$coordinate$scenario,x$coordinate$replicate),character(1))
  names(bundles) <- keys
  dispatch <- getOption("sblrbench.parameter_fit_dispatch",
    fit_parameter_estimation_method)
  if(!is.function(dispatch)) stop("Parameter fit dispatch must be a function.",
    call.=FALSE)
  status_rows <- runtime_rows <- estimate_rows <- marker_rows <- list()
  convergence_rows <- list()
  for(i in seq_len(nrow(coordinates))) {
    coordinate <- coordinates[i,,drop=FALSE]; cl <- as.list(coordinate)
    bundle <- bundles[[.prediction_bundle_key(coordinate$scenario,
      coordinate$replicate)]]; method <- methods[[coordinate$method]]
    controls <- benchmark_method_controls(spec,coordinate$method,profile,
      coordinate$fit_seed,coordinate$chain_seeds[[1L]])
    identity <- list(schema_version=1L,study=spec$study,coordinate=cl,
      sample_ids=data$sample_ids,marker_ids=data$markers$marker_ids,
      simulation_effects=bundle$simulation$truth$effects,
      simulation_phenotypes=bundle$simulation$truth$phenotypes,
      controls=controls,package_sha=spec$packages$sblr$sha,
      qgdata_commit=spec$data$example_data$commit)
    input_hash <- benchmark_hash_object(identity)
    checkpoint <- .prediction_checkpoint_path(paths,cl)
    result <- NULL; reused <- FALSE; reason <- ""
    if(isTRUE(resume) && file.exists(checkpoint)) {
      loaded <- benchmark_load_checkpoint(checkpoint,input_hash,
        validator=function(x) identical(x$package_sha,spec$packages$sblr$sha) &&
          identical(x$coordinate,cl) && inherits(x$result,"sblrbench_result"))
      result <- loaded$value$result; reused <- TRUE
    } else {
      result <- tryCatch(dispatch(method=method,controls=controls,
        simulation=bundle$simulation,stats=bundle$stats,glist=data$ld_glist),
        error=function(e) {reason <<- conditionMessage(e); NULL})
      if(!is.null(result)) benchmark_atomic_save_rds(list(input_hash=input_hash,
        package_sha=spec$packages$sblr$sha,coordinate=cl,result=result),checkpoint,
        compress=FALSE,temporary_prefix=".parameter-fit-")
    }
    ok <- !is.null(result)
    status_rows[[i]] <- data.frame(study=spec$study,
      scenario=coordinate$scenario,replicate=coordinate$replicate,
      method=coordinate$method,status=if(ok)"ok" else "failed",reason=reason,
      reused=reused,input_hash=input_hash,stringsAsFactors=FALSE)
    runtime_rows[[i]] <- data.frame(study=spec$study,
      scenario=coordinate$scenario,replicate=coordinate$replicate,
      method=coordinate$method,elapsed_seconds=if(ok)extract_runtime(result) else NA_real_,
      status=if(ok)"ok" else "failed",reason=reason,reused=reused,
      stringsAsFactors=FALSE)
    if(ok) {
      draws <- extract_parameter_draws(result,coordinate$method,spec$estimands,
        length(data$markers$marker_ids),expected_chains=controls$nchains)
      summary <- complete_parameter_summary(summarise_parameter_draws(draws),
        spec$estimands,coordinate$method,controls$nchains)
      estimate_rows[[i]] <- parameter_recovery_metrics(summary,bundle$truth,
        method,spec$estimands,spec$validation$relative_error_tolerance)
      marker_rows[[i]] <- prediction_marker_table(result,coordinate$scenario,
        coordinate$replicate,coordinate$method)
      convergence_rows[[i]] <- .prediction_convergence_row(result,cl)
    }
  }
  status <- .benchmark_bind_rows(status_rows); runtime <- .benchmark_bind_rows(runtime_rows)
  estimates <- .benchmark_bind_rows(estimate_rows)
  marker_results <- .benchmark_bind_rows(marker_rows)
  convergence <- .benchmark_bind_rows(convergence_rows)
  truth <- .benchmark_bind_rows(lapply(bundles,`[[`,"truth"))
  oracle <- .benchmark_oracle_table(bundles, spec$study)
  metrics <- if(is.null(estimates)) NULL else parameter_recovery_summary(estimates)
  paired <- if(is.null(estimates)) NULL else parameter_paired_differences(estimates)
  paired_summary <- if(is.null(paired)) NULL else parameter_paired_summary(paired)
  files <- list(
    fit_status=.benchmark_write_csv(status,file.path(paths$tables,"fit_status.csv")),
    simulation_truth=.benchmark_write_csv(truth,file.path(paths$tables,"simulation_truth.csv")),
    simulation_oracle=.benchmark_write_csv(oracle,
      file.path(paths$tables,"simulation_oracle.csv")),
    estimates=.benchmark_write_csv(estimates,file.path(paths$tables,"estimates.csv")),
    parameter_metrics=.benchmark_write_csv(metrics,file.path(paths$tables,"parameter_metrics.csv")),
    marker_results=.benchmark_write_csv(marker_results,file.path(paths$tables,"marker_results.csv")),
    convergence=.benchmark_write_csv(convergence,file.path(paths$tables,"convergence.csv")),
    runtime=.benchmark_write_csv(runtime,file.path(paths$tables,"runtime.csv")),
    paired=.benchmark_write_csv(paired,file.path(paths$tables,"paired_parameter_differences.csv")),
    paired_summary=.benchmark_write_csv(paired_summary,file.path(paths$tables,"paired_comparison_summary.csv")))
  manifest <- .prediction_manifest(spec,profile,paths,coordinates,FALSE,status,data)
  .write_prediction_manifest(manifest,paths$manifest)
  writeLines(benchmark_session_information(),paths$session_info)
  list(spec=spec,paths=c(paths,files),status=status,truth=truth,
    estimates=estimates,marker_results=marker_results,metrics=metrics,
    convergence=convergence,runtime=runtime,oracle=oracle)
}

.prediction_checkpoint_identity <- function(spec, coordinate, controls,
                                            simulation, data) {
  list(schema_version = 1L, study = spec$study, task = spec$task,
    scenario = coordinate$scenario, replicate = as.integer(coordinate$replicate),
    method = coordinate$method, split_seed = data$split$split_seed,
    train_ids = data$split$train_ids, test_ids = data$split$test_ids,
    marker_ids = data$markers$marker_ids,
    allele_frequency = data$scaled$allele_frequency,
    simulation_seed = as.integer(coordinate$simulation_seed),
    simulation_effects = simulation$truth$effects,
    simulation_phenotypes = simulation$truth$phenotypes,
    controls = controls, package_sha = spec$packages$sblr$sha,
    qgdata_commit = spec$data$example_data$commit)
}

.prediction_checkpoint_path <- function(paths, coordinate) {
  file.path(paths$checkpoints, "fits", coordinate$scenario,
    paste0("replicate-", coordinate$replicate),
    paste0(coordinate$method, ".rds"))
}

.prediction_truth_table <- function(bundle) {
  effects <- bundle$simulation$truth$effects
  data.frame(study = bundle$simulation$scenario$study,
    scenario = bundle$simulation$scenario$architecture,
    replicate = bundle$simulation$scenario$replicate,
    marker = rep(rownames(effects), times = ncol(effects)),
    trait = rep(colnames(effects), each = nrow(effects)),
    true_effect = as.numeric(effects),
    causal = rep(rownames(effects) %in% bundle$simulation$truth$causal$all,
      times = ncol(effects)),
    simulation_seed = bundle$simulation$provenance$seed,
    stringsAsFactors = FALSE)
}

.prediction_convergence_row <- function(result, coordinate) {
  information <- extract_chain_information(result)
  compact <- information$compact_summaries
  data.frame(study = "02_prediction", scenario = coordinate$scenario,
    replicate = coordinate$replicate, method = coordinate$method,
    true_traces_available = !is.null(information$true_traces),
    compact_summary_available = !is.null(compact),
    final_states_available = length(information$final_states) > 0L,
    compact_summary = if (is.null(compact)) NA_character_ else
      jsonlite::toJSON(compact, auto_unbox = TRUE, null = "null"),
    stringsAsFactors = FALSE)
}

.convergence_checkpoint_path <- function(paths, coordinate) {
  file.path(paths$checkpoints, "fits", coordinate$stage,
    coordinate$scenario, paste0("replicate-", coordinate$replicate),
    paste0(coordinate$method, ".rds"))
}

.convergence_checkpoint_identity <- function(spec, parameter_spec,
                                             coordinate, controls,
                                             simulation, data) {
  benchmark_semantic_checkpoint_identity(
    diagnostic_id = "study04-convergence-fit",
    scientific_inputs = list(
      study = spec$study, task = spec$task, stage = coordinate$stage,
      scenario = coordinate$scenario,
      replicate = as.integer(coordinate$replicate), method = coordinate$method,
      sample_ids = data$sample_ids, marker_ids = data$markers$marker_ids,
      simulation_seed = as.integer(coordinate$simulation_seed),
      simulation_effects = simulation$truth$effects,
      simulation_phenotypes = simulation$truth$phenotypes,
      controls = controls, diagnostic_registry = spec$diagnostics$registry,
      thresholds = spec$diagnostics$thresholds,
      candidates = spec$diagnostics[c("burnin_candidates",
        "retained_draw_candidates")],
      parameter_priors = parameter_spec$controls$priors,
      sblr_sha = spec$packages$sblr$sha,
      qgdata_sha = spec$packages$qgdata$sha))
}

.convergence_fit <- function(spec, parameter_spec, coordinate, bundle, data,
                             controls, paths, resume, dispatch) {
  coordinate_list <- as.list(coordinate)
  identity <- .convergence_checkpoint_identity(spec, parameter_spec,
    coordinate_list, controls, bundle$simulation, data)
  input_hash <- benchmark_semantic_checkpoint_hash(identity)
  checkpoint <- .convergence_checkpoint_path(paths, coordinate_list)
  result <- NULL
  reused <- FALSE
  reason <- ""
  if (isTRUE(resume) && file.exists(checkpoint)) {
    loaded <- benchmark_load_semantic_checkpoint(checkpoint, input_hash,
      validator = function(x) inherits(x$result, "sblrbench_result"))
    result <- loaded$value$result
    reused <- TRUE
  } else {
    method <- parameter_spec$methods[[coordinate$method]]
    method$id <- coordinate$method
    result <- tryCatch(dispatch(method = method, controls = controls,
      simulation = bundle$simulation, stats = bundle$stats,
      glist = data$ld_glist), error = function(error) {
        reason <<- conditionMessage(error)
        NULL
      })
    if (!is.null(result))
      benchmark_atomic_save_rds(list(
        checkpoint_schema = "sblrbench-semantic-v2",
        identity_payload = identity, semantic_hash = input_hash,
        result = result), checkpoint, compress = FALSE,
        temporary_prefix = ".convergence-fit-")
  }
  list(result = result, reused = reused, reason = reason,
    input_hash = input_hash, checkpoint = checkpoint)
}

.convergence_status_row <- function(spec, coordinate, fit) {
  data.frame(study = spec$study, stage = coordinate$stage,
    scenario = coordinate$scenario, replicate = coordinate$replicate,
    method = coordinate$method,
    status = if (is.null(fit$result)) "failed" else "ok",
    reason = fit$reason, reused = fit$reused, input_hash = fit$input_hash,
    stringsAsFactors = FALSE)
}

.convergence_runtime_row <- function(spec, coordinate, fit) {
  data.frame(study = spec$study, stage = coordinate$stage,
    scenario = coordinate$scenario, replicate = coordinate$replicate,
    method = coordinate$method,
    elapsed_seconds = if (is.null(fit$result)) NA_real_ else
      extract_runtime(fit$result),
    status = if (is.null(fit$result)) "failed" else "ok",
    reason = fit$reason, reused = fit$reused, stringsAsFactors = FALSE)
}

.convergence_bundle <- function(bundles, scenario, replicate) {
  key <- .prediction_bundle_key(scenario, replicate)
  bundle <- bundles[[key]]
  if (is.null(bundle)) stop("No matched Study 03 simulation bundle for ", key,
    ".", call. = FALSE)
  bundle
}

.run_convergence <- function(spec, profile, resolved, coordinates, paths,
                             resume) {
  parameter_spec <- benchmark_matched_spec(spec)
  data <- prepare_parameter_estimation_data(parameter_spec, paths$root)
  parameter_profile <- if (resolved$replicate_count == 1L) "workshop" else
    "benchmark"
  bundles <- prepare_parameter_estimation_simulations(parameter_spec,
    parameter_profile, data)
  names(bundles) <- vapply(bundles, function(x) .prediction_bundle_key(
    x$coordinate$scenario, x$coordinate$replicate), character(1))
  dispatch <- getOption("sblrbench.convergence_fit_dispatch",
    fit_parameter_estimation_method)
  if (!is.function(dispatch))
    stop("Convergence fit dispatch must be a function.", call. = FALSE)

  selection <- coordinates[coordinates$stage == "selection", , drop = FALSE]
  status_rows <- runtime_rows <- selection_draws <- diagnostics <- stability <-
    list()
  for (i in seq_len(nrow(selection))) {
    coordinate <- selection[i, , drop = FALSE]
    controls <- convergence_method_controls(spec, parameter_spec,
      coordinate)
    bundle <- .convergence_bundle(bundles, coordinate$scenario,
      coordinate$replicate)
    fit <- .convergence_fit(spec, parameter_spec, coordinate, bundle, data,
      controls, paths, resume, dispatch)
    status_rows[[length(status_rows) + 1L]] <-
      .convergence_status_row(spec, coordinate, fit)
    runtime_rows[[length(runtime_rows) + 1L]] <-
      .convergence_runtime_row(spec, coordinate, fit)
    if (!is.null(fit$result)) {
      draws <- extract_convergence_traces(fit$result, coordinate,
        spec$diagnostics$registry,
        spec$diagnostics$thresholds$chain_count)
      selection_draws[[length(selection_draws) + 1L]] <- draws
      diagnostics[[length(diagnostics) + 1L]] <-
        benchmark_convergence_candidate_grid(draws, spec)
      stability[[length(stability) + 1L]] <-
        benchmark_burnin_stability(draws, spec)
    }
  }
  selection_diagnostics <- .benchmark_bind_rows(diagnostics)
  burnin_stability <- .benchmark_bind_rows(stability)
  recommendations <- if (!is.null(selection_diagnostics) &&
      !is.null(burnin_stability))
    benchmark_convergence_recommendations(selection_diagnostics,
      burnin_stability, spec) else NULL
  if (nrow(selection) && (is.null(recommendations) ||
      nrow(recommendations) != nrow(selection)))
    stop("Selection fits did not produce a complete recommendation grid.",
      call. = FALSE)

  validation <- coordinates[coordinates$stage == "validation", , drop = FALSE]
  validation_diagnostics <- validation_status <- list()
  if (nrow(validation)) for (i in seq_len(nrow(validation))) {
    coordinate <- validation[i, , drop = FALSE]
    controls <- convergence_method_controls(spec, parameter_spec,
      coordinate, recommendations)
    bundle <- .convergence_bundle(bundles, coordinate$scenario,
      coordinate$replicate)
    fit <- .convergence_fit(spec, parameter_spec, coordinate, bundle, data,
      controls, paths, resume, dispatch)
    status_row <- .convergence_status_row(spec, coordinate, fit)
    status_rows[[length(status_rows) + 1L]] <- status_row
    runtime_rows[[length(runtime_rows) + 1L]] <-
      .convergence_runtime_row(spec, coordinate, fit)
    diagnostic <- NULL
    if (!is.null(fit$result)) {
      draws <- extract_convergence_traces(fit$result, coordinate,
        spec$diagnostics$registry,
        spec$diagnostics$thresholds$chain_count)
      diagnostic <- benchmark_convergence_diagnostics(draws, 0L,
        as.integer(controls$nit), spec$diagnostics$registry,
        spec$diagnostics$thresholds)
      validation_diagnostics[[length(validation_diagnostics) + 1L]] <-
        diagnostic
    }
    limiting <- if (is.null(diagnostic)) character() else
      diagnostic$quantity[!diagnostic$overall_pass]
    validation_status[[length(validation_status) + 1L]] <- data.frame(
      scenario = coordinate$scenario, replicate = coordinate$replicate,
      method = coordinate$method, status = status_row$status,
      all_core_quantities_pass = !is.null(diagnostic) &&
        nrow(diagnostic) == sum(spec$diagnostics$registry$required) &&
        all(diagnostic$overall_pass),
      limiting_quantities = if (length(limiting))
        paste(limiting, collapse = ";") else "none",
      maximum_rhat = if (is.null(diagnostic)) NA_real_ else
        max(diagnostic$rhat),
      minimum_bulk_ess = if (is.null(diagnostic)) NA_real_ else
        min(diagnostic$ess_bulk),
      minimum_tail_ess = if (is.null(diagnostic)) NA_real_ else
        min(diagnostic$ess_tail),
      maximum_relative_mcse = if (is.null(diagnostic)) NA_real_ else
        max(diagnostic$relative_mcse), reason = fit$reason,
      stringsAsFactors = FALSE)
  }
  validation_diagnostics <- .benchmark_bind_rows(validation_diagnostics)
  validation_status <- .benchmark_bind_rows(validation_status)
  validation_summary <- if (is.null(validation_status)) NULL else
    benchmark_convergence_validation_summary(validation_status,
      validation_diagnostics,
      spec$validation$required_successful_replicates)
  convergence <- .benchmark_bind_rows(list(selection_diagnostics,
    validation_diagnostics))
  status <- .benchmark_bind_rows(status_rows)
  runtime <- .benchmark_bind_rows(runtime_rows)
  candidate_summary <- if (is.null(selection_diagnostics)) NULL else
    aggregate(overall_pass ~ scenario + method + burnin_candidate +
      retained_draw_candidate, selection_diagnostics, all)
  files <- list(
    fit_status = .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv")),
    coordinate_grid = .benchmark_write_csv(coordinates,
      file.path(paths$tables, "coordinate_grid.csv")),
    convergence = .benchmark_write_csv(convergence,
      file.path(paths$tables, "convergence.csv")),
    candidate_summary = .benchmark_write_csv(candidate_summary,
      file.path(paths$tables, "candidate_summary.csv")),
    burnin_stability = .benchmark_write_csv(burnin_stability,
      file.path(paths$tables, "burnin_stability.csv")),
    recommendations = .benchmark_write_csv(recommendations,
      file.path(paths$tables, "recommendations.csv")),
    validation_replicates = .benchmark_write_csv(validation_status,
      file.path(paths$tables, "validation_replicates.csv")),
    validation_summary = .benchmark_write_csv(validation_summary,
      file.path(paths$tables, "validation_summary.csv")),
    runtime = .benchmark_write_csv(runtime,
      file.path(paths$tables, "runtime.csv")))
  manifest <- .prediction_manifest(spec, profile, paths, coordinates, FALSE,
    status, data)
  .write_prediction_manifest(manifest, paths$manifest)
  writeLines(benchmark_session_information(), paths$session_info)
  list(spec = spec, paths = c(paths, files), status = status,
    coordinate_grid = coordinates, convergence = convergence,
    candidate_summary = candidate_summary,
    burnin_stability = burnin_stability,
    recommendations = recommendations,
    validation_replicates = validation_status,
    validation_summary = validation_summary, runtime = runtime)
}

.annotation_logic <- function(spec) {
  path <- benchmark_spec_path(spec, spec$annotation_design$implementation)
  if (!file.exists(path))
    stop("Study 06 annotation implementation is missing: ", path,
      call. = FALSE)
  environment <- new.env(parent = environment())
  sys.source(path, envir = environment)
  environment
}

benchmark_annotation_spec_hash <- function(spec) {
  x <- spec
  attributes(x) <- NULL
  benchmark_hash_object(x)
}

annotation_qualification_artifact_schema <- function() list(
  schema = "sblrbench-annotation-qualification-v1",
  required_fields = c("study", "spec_hash", "sblr_sha", "qgdata_sha",
    "entries", "overall_decision", "created_at"),
  entry_fields = c("scenario", "replicate", "method",
    "available_history", "selected_burnin", "selected_retained", "rhat",
    "ess_bulk", "ess_tail", "relative_mcse", "all_quantities_pass",
    "quantity_decisions", "semantic_checkpoint_hash",
    "reusable_history_hash"),
  quantity_fields = c("quantity", "rhat", "ess_bulk", "ess_tail",
    "relative_mcse", "pass"))

validate_annotation_qualification_decision <- function(decision, spec) {
  schema <- annotation_qualification_artifact_schema()
  if (!is.list(decision) ||
      !all(schema$required_fields %in% names(decision)) ||
      !identical(decision$schema, schema$schema))
    stop("Study 06 qualification decision artifact is incomplete.",
      call. = FALSE)
  if (!identical(decision$study, spec$study) ||
      !identical(decision$spec_hash, benchmark_annotation_spec_hash(spec)) ||
      !identical(decision$sblr_sha, spec$packages$sblr$sha) ||
      !identical(decision$qgdata_sha, spec$packages$qgdata$sha))
    stop("Study 06 qualification decision identity is stale or incompatible.",
      call. = FALSE)
  if (!is.list(decision$entries) || length(decision$entries) != 4L ||
      any(!vapply(decision$entries, function(x)
        all(schema$entry_fields %in% names(x)), logical(1))))
    stop("Study 06 qualification decision must contain four complete entries.",
      call. = FALSE)
  if (any(!vapply(decision$entries, function(entry)
      is.list(entry$quantity_decisions) && length(entry$quantity_decisions) &&
        all(vapply(entry$quantity_decisions, function(quantity)
          all(schema$quantity_fields %in% names(quantity)), logical(1))),
      logical(1))))
    stop("Study 06 qualification decision lacks per-quantity diagnostics.",
      call. = FALSE)
  expected <- spec$qualification$entries
  observed <- do.call(rbind, lapply(decision$entries, function(x)
    data.frame(scenario = x$scenario, replicate = as.integer(x$replicate),
      method = x$method, stringsAsFactors = FALSE)))
  rownames(observed) <- NULL
  if (!identical(observed, expected))
    stop("Study 06 qualification entries differ from the specification.",
      call. = FALSE)
  if (!identical(decision$overall_decision, "passed") ||
      any(!vapply(decision$entries, function(x)
        isTRUE(x$all_quantities_pass), logical(1))))
    stop("Study 06 qualification did not pass every required entry.",
      call. = FALSE)
  invisible(decision)
}

.read_annotation_qualification_decision <- function(spec, paths) {
  path <- file.path(paths$root, spec$qualification$decision_path)
  if (!file.exists(path))
    stop("Final Study 06 mode requires a passing qualification decision: ",
      path, call. = FALSE)
  decision <- jsonlite::read_json(path, simplifyVector = FALSE)
  validate_annotation_qualification_decision(decision, spec)
  decision
}

.annotation_checkpoint_identities <- function(spec, coordinate, controls,
                                              data, bundle, mode) {
  annotation_hash <- benchmark_hash_object(bundle$annotations)
  common <- list(study = spec$study, task = spec$task,
    scenario = as.character(coordinate$scenario),
    replicate = as.integer(coordinate$replicate),
    method = as.character(coordinate$method),
    train_sample_order_hash = benchmark_hash_object(data$split$train_ids),
    test_sample_order_hash = benchmark_hash_object(data$split$test_ids),
    marker_order_hash = benchmark_hash_object(data$markers$marker_ids),
    annotation_matrix_hash = annotation_hash,
    annotation_columns = colnames(bundle$annotations),
    phenotype_hash = benchmark_hash_object(bundle$simulation$truth$phenotypes),
    truth_effect_hash = benchmark_hash_object(bundle$simulation$truth$effects),
    priors = spec$controls$priors, controls = controls,
    mcmc_history = controls[c("nit", "nburn", "nthin", "nchains")],
    seeds = list(component = as.integer(coordinate$component_seed),
      effect = as.integer(coordinate$effect_seed),
      residual = as.integer(coordinate$residual_seed),
      fit = as.integer(coordinate$fit_seed),
      chains = as.integer(coordinate$chain_seeds[[1L]])),
    ld_settings = spec$data$sparse_ld,
    sblr_sha = spec$packages$sblr$sha,
    qgdata_sha = spec$packages$qgdata$sha)
  checkpoint <- benchmark_semantic_checkpoint_identity(
    "study06-annotation-fit", c(list(mode = mode), common))
  history <- benchmark_semantic_checkpoint_identity(
    "study06-annotation-history", common)
  list(checkpoint = checkpoint,
    checkpoint_hash = benchmark_semantic_checkpoint_hash(checkpoint),
    history = history,
    history_hash = benchmark_semantic_checkpoint_hash(history))
}

.annotation_checkpoint_path <- function(paths, coordinate, mode) {
  directory <- if (identical(mode, "qualification"))
    file.path(paths$qualification, "checkpoints") else paths$checkpoints
  file.path(directory, paste0(coordinate$scenario, "--r",
    coordinate$replicate, "--", coordinate$method, ".rds"))
}

.annotation_reusable_qualification_history <- function(spec, coordinate,
    controls, data, bundle, paths, decision, resume) {
  if (!isTRUE(resume) || as.integer(coordinate$replicate) != 1L)
    return(NULL)
  method <- spec$methods[[as.character(coordinate$method)]]
  if (!isTRUE(method$annotation_aware)) return(NULL)
  hit <- vapply(decision$entries, function(x)
    identical(x$scenario, as.character(coordinate$scenario)) &&
      identical(x$method, as.character(coordinate$method)), logical(1))
  if (sum(hit) != 1L) return(NULL)
  selected <- decision$entries[[which(hit)]]
  if (selected$selected_burnin + selected$selected_retained >
      selected$available_history) return(NULL)
  qualification_controls <- annotation_method_controls(spec, coordinate,
    profile = "benchmark", mode = "qualification")
  identity <- .annotation_checkpoint_identities(spec, coordinate,
    qualification_controls, data, bundle, "qualification")
  if (!identical(selected$semantic_checkpoint_hash,
      identity$checkpoint_hash) ||
      !identical(selected$reusable_history_hash, identity$history_hash))
    return(NULL)
  path <- .annotation_checkpoint_path(paths, coordinate, "qualification")
  if (!file.exists(path)) return(NULL)
  payload <- readRDS(path)
  if (!is.list(payload) ||
      !identical(payload$checkpoint_schema, "sblrbench-semantic-v2") ||
      !identical(payload$semantic_hash, identity$checkpoint_hash) ||
      !identical(payload$reusable_history_hash, identity$history_hash))
    return(NULL)
  key <- paste0("burnin_", selected$selected_burnin, "--retained_",
    selected$selected_retained)
  window <- payload$qualified_windows[[key]]
  if (!is.list(window) ||
      !identical(as.integer(window$burnin),
        as.integer(selected$selected_burnin)) ||
      !identical(as.integer(window$retained),
        as.integer(selected$selected_retained)) ||
      !inherits(window$result, "sblrbench_result")) return(NULL)
  # A reusable window must contain the same final-output contract, not merely
  # a maximum-history fit or convergence summaries.
  valid <- tryCatch({
    extract_marker_effects(window$result)
    extract_marker_probabilities(window$result)
    traces <- extract_annotation_coefficient_traces(window$result,
      expected_chains = controls$nchains)
    identical(traces$status[[1L]], "ok")
  }, error = function(...) FALSE)
  if (!isTRUE(valid)) return(NULL)
  list(result = window$result, reused = TRUE, reason = "",
    checkpoint = path,
    identities = .annotation_checkpoint_identities(spec, coordinate,
      controls, data, bundle, "final"))
}

.annotation_fit <- function(spec, coordinate, method, controls, data, bundle,
                            paths, mode, resume, dispatch) {
  identities <- .annotation_checkpoint_identities(spec, coordinate, controls,
    data, bundle, mode)
  checkpoint <- .annotation_checkpoint_path(paths, coordinate, mode)
  result <- NULL
  reused <- FALSE
  reason <- ""
  if (isTRUE(resume) && file.exists(checkpoint)) {
    loaded <- benchmark_load_semantic_checkpoint(checkpoint,
      identities$checkpoint_hash, validator = function(x)
        inherits(x$result, "sblrbench_result") &&
          identical(x$reusable_history_hash, identities$history_hash))
    result <- loaded$value$result
    reused <- TRUE
  } else {
    result <- tryCatch(dispatch(method = method, controls = controls,
      simulation = bundle$simulation, stats = bundle$stats,
      glist = data$ld_glist, split = data$split,
      annotations = bundle$annotations,
      annotation_truth = bundle$annotation_truth), error = function(error) {
        reason <<- conditionMessage(error)
        NULL
      })
    if (!is.null(result)) benchmark_atomic_save_rds(list(
      checkpoint_schema = "sblrbench-semantic-v2",
      identity_payload = identities$checkpoint,
      semantic_hash = identities$checkpoint_hash,
      reusable_history_identity = identities$history,
      reusable_history_hash = identities$history_hash,
      result = result), checkpoint, compress = FALSE,
      temporary_prefix = ".annotation-fit-")
  }
  list(result = result, reused = reused, reason = reason,
    checkpoint = checkpoint, identities = identities)
}

.annotation_required_traces <- function(result, coordinate, bundle,
                                        expected_chains) {
  alpha <- extract_annotation_coefficient_traces(result, expected_chains)
  if (any(alpha$status != "ok"))
    stop(alpha$reason[[1L]], call. = FALSE)
  alpha$quantity <- ifelse(alpha$parameter == "alpha",
    paste("alpha", alpha$annotation, alpha$stick, sep = ":"),
    paste("sigmaSqAlpha", alpha$stick, sep = ":"))
  out <- alpha[c("iteration", "chain", "quantity", "value")]
  native <- .benchmark_native_fit(result)
  values <- native$convergence_traces$values
  descriptors <- native$convergence_traces$quantities
  groups <- as.character(descriptors$group)
  index <- match(c("vbs", "vgs", "ves"), groups)
  if (anyNA(index))
    stop("Study 06 qualification requires true vbs, vgs, and ves traces.",
      call. = FALSE)
  base <- expand.grid(iteration = seq_len(dim(values)[1L]),
    chain = seq_len(dim(values)[2L]))
  core <- do.call(rbind, lapply(seq_along(index), function(i) data.frame(base,
    quantity = c("effect_variance", "genetic_variance",
      "residual_variance")[i], value = as.vector(values[, , index[i]]),
    stringsAsFactors = FALSE)))
  h2 <- values[, , index[2L]] /
    (values[, , index[2L]] + values[, , index[3L]])
  core <- rbind(core, data.frame(base, quantity = "heritability",
    value = as.vector(h2), stringsAsFactors = FALSE))
  prior <- summarise_drawwise_annotation_prior(alpha, bundle$annotations,
    bundle$spec$controls$simulation$mixture_var,
    bundle$simulation$extras$true_marker_prior, bundle$marker_truth)
  if (!identical(prior$status, "ok"))
    stop(prior$reason, call. = FALSE)
  derived <- reshape(prior$draws, varying = c("expected_active",
    "mean_prior_enriched", "mean_prior_unannotated",
    "enriched_prior_contrast"), v.names = "value",
    timevar = "quantity", times = c("prior_expected_active",
      "prior_nonnull_mean_enriched", "prior_nonnull_mean_unannotated",
      "prior_enriched_contrast"), direction = "long")
  derived <- derived[c("iteration", "chain", "quantity", "value")]
  out <- rbind(out, core, derived)
  out$scenario <- as.character(coordinate$scenario)
  out$replicate <- as.integer(coordinate$replicate)
  out$method <- as.character(coordinate$method)
  out
}

.annotation_candidate_diagnostics <- function(draws, spec) {
  rows <- list()
  for (burnin in spec$qualification$burnin_candidates)
    for (retained in spec$qualification$retained_candidates)
      if (burnin + retained <= spec$qualification$maximum_history) {
        z <- draws[draws$iteration > burnin &
          draws$iteration <= burnin + retained, , drop = FALSE]
        counts <- table(z$quantity, z$chain)
        if (!nrow(counts) || ncol(counts) != 4L || any(counts != retained))
          stop("Qualification candidate lacks equal draws in four chains.",
            call. = FALSE)
        for (quantity in unique(z$quantity)) {
          x <- z[z$quantity == quantity, , drop = FALSE]
          diagnostic <- benchmark_scalar_diagnostics(x,
            spec$qualification$thresholds)
          rows[[length(rows) + 1L]] <- cbind(
            x[1L, c("scenario", "replicate", "method"), drop = FALSE],
            burnin = burnin, retained = retained, quantity = quantity,
            diagnostic, stringsAsFactors = FALSE)
        }
      }
  do.call(rbind, rows)
}

.annotation_select_candidate <- function(diagnostics) {
  key <- interaction(diagnostics$burnin, diagnostics$retained, drop = TRUE)
  summary <- do.call(rbind, lapply(split(diagnostics, key), function(x)
    data.frame(burnin = x$burnin[1L], retained = x$retained[1L],
      all_quantities_pass = all(x$overall_pass),
      maximum_rhat = max(x$rhat), minimum_bulk_ess = min(x$ess_bulk),
      minimum_tail_ess = min(x$ess_tail),
      maximum_relative_mcse = max(x$relative_mcse),
      stringsAsFactors = FALSE)))
  summary <- summary[order(summary$burnin + summary$retained,
    summary$burnin, summary$retained), , drop = FALSE]
  pass <- summary[summary$all_quantities_pass, , drop = FALSE]
  selected <- if (nrow(pass)) pass[1L, , drop = FALSE] else
    summary[nrow(summary), , drop = FALSE]
  list(summary = summary, selected = selected)
}

.annotation_validation_design <- function(spec) {
  logic <- .annotation_logic(spec)
  marker_ids <- sprintf("marker_%05d", seq_len(spec$data$expected_marker_count))
  annotations <- logic$construct_annotation_design(marker_ids, spec)
  truth <- logic$construct_annotation_truth(annotations, spec)
  list(annotations = annotations, truth = truth,
    summary = logic$annotation_design_summary(annotations, truth, spec))
}

.annotation_bundle_map <- function(bundles) {
  keys <- vapply(bundles, function(x) .prediction_bundle_key(
    x$coordinate$scenario, x$coordinate$replicate), character(1))
  if (anyDuplicated(keys))
    stop("Study 06 simulation bundle coordinates are duplicated.",
      call. = FALSE)
  names(bundles) <- keys
  bundles
}

.run_annotation_qualification <- function(spec, profile, coordinates, paths,
                                          resume) {
  data <- prepare_prediction_data(spec, paths$root)
  if (length(data$sample_ids) != spec$validation$expected_sample_count ||
      length(data$markers$marker_ids) != spec$validation$expected_marker_count)
    stop("Study 06 prepared data counts differ from the specification.",
      call. = FALSE)
  logic <- .annotation_logic(spec)
  bundles <- logic$prepare_annotation_simulations(spec, profile, data)
  bundles <- .annotation_bundle_map(bundles)
  methods <- resolve_benchmark_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")
  dispatch <- getOption("sblrbench.annotation_fit_dispatch",
    fit_annotation_method)
  status_rows <- runtime_rows <- diagnostic_rows <- candidate_rows <- list()
  decision_entries <- vector("list", nrow(coordinates))
  for (i in seq_len(nrow(coordinates))) {
    coordinate <- coordinates[i, , drop = FALSE]
    bundle <- bundles[[.prediction_bundle_key(coordinate$scenario,
      coordinate$replicate)]]
    bundle$spec <- spec
    controls <- annotation_method_controls(spec, coordinate, profile,
      mode = "qualification")
    fit <- .annotation_fit(spec, coordinate, methods[[coordinate$method]],
      controls, data, bundle, paths, "qualification", resume, dispatch)
    ok <- !is.null(fit$result)
    status_rows[[i]] <- data.frame(study = spec$study,
      mode = "qualification", scenario = coordinate$scenario,
      replicate = coordinate$replicate, method = coordinate$method,
      status = if (ok) "ok" else "failed", reason = fit$reason,
      reused = fit$reused, input_hash = fit$identities$checkpoint_hash,
      stringsAsFactors = FALSE)
    runtime_rows[[i]] <- data.frame(study = spec$study,
      mode = "qualification", scenario = coordinate$scenario,
      replicate = coordinate$replicate, method = coordinate$method,
      elapsed_seconds = if (ok) extract_runtime(fit$result) else NA_real_,
      status = if (ok) "ok" else "failed", reason = fit$reason,
      reused = fit$reused, stringsAsFactors = FALSE)
    selected <- selected_diagnostics <- NULL
    if (ok) {
      draws <- .annotation_required_traces(fit$result, coordinate, bundle,
        spec$qualification$nchains)
      diagnostics <- .annotation_candidate_diagnostics(draws, spec)
      choice <- .annotation_select_candidate(diagnostics)
      diagnostic_rows[[i]] <- diagnostics
      candidate_rows[[i]] <- cbind(
        coordinate[c("scenario", "replicate", "method")], choice$summary)
      selected <- choice$selected
      selected_diagnostics <- diagnostics[
        diagnostics$burnin == selected$burnin &
          diagnostics$retained == selected$retained, , drop = FALSE]
    }
    quantity_decisions <- if (is.null(selected_diagnostics)) list() else
      lapply(seq_len(nrow(selected_diagnostics)), function(j) list(
        quantity = as.character(selected_diagnostics$quantity[j]),
        rhat = selected_diagnostics$rhat[j],
        ess_bulk = selected_diagnostics$ess_bulk[j],
        ess_tail = selected_diagnostics$ess_tail[j],
        relative_mcse = selected_diagnostics$relative_mcse[j],
        pass = isTRUE(selected_diagnostics$overall_pass[j])))
    decision_entries[[i]] <- list(scenario = as.character(coordinate$scenario),
      replicate = as.integer(coordinate$replicate),
      method = as.character(coordinate$method),
      available_history = if (ok) spec$qualification$maximum_history else 0L,
      selected_burnin = if (is.null(selected)) NA_integer_ else
        as.integer(selected$burnin),
      selected_retained = if (is.null(selected)) NA_integer_ else
        as.integer(selected$retained),
      rhat = if (is.null(selected)) NA_real_ else selected$maximum_rhat,
      ess_bulk = if (is.null(selected)) NA_real_ else selected$minimum_bulk_ess,
      ess_tail = if (is.null(selected)) NA_real_ else selected$minimum_tail_ess,
      relative_mcse = if (is.null(selected)) NA_real_ else
        selected$maximum_relative_mcse,
      all_quantities_pass = !is.null(selected) &&
        isTRUE(selected$all_quantities_pass),
      quantity_decisions = quantity_decisions,
      semantic_checkpoint_hash = fit$identities$checkpoint_hash,
      reusable_history_hash = fit$identities$history_hash)
  }
  status <- .benchmark_bind_rows(status_rows)
  runtime <- .benchmark_bind_rows(runtime_rows)
  convergence <- .benchmark_bind_rows(diagnostic_rows)
  candidates <- .benchmark_bind_rows(candidate_rows)
  passed <- all(vapply(decision_entries, function(x)
    isTRUE(x$all_quantities_pass), logical(1)))
  decision <- list(schema = "sblrbench-annotation-qualification-v1",
    study = spec$study, spec_hash = benchmark_annotation_spec_hash(spec),
    sblr_sha = spec$packages$sblr$sha,
    qgdata_sha = spec$packages$qgdata$sha, entries = decision_entries,
    overall_decision = if (passed) "passed" else "failed",
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  decision_path <- file.path(paths$root, spec$qualification$decision_path)
  dir.create(dirname(decision_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(decision, decision_path, pretty = TRUE,
    auto_unbox = TRUE, null = "null", na = "null")
  files <- list(fit_status = .benchmark_write_csv(status,
    file.path(paths$qualification, "fit_status.csv")),
    convergence = .benchmark_write_csv(convergence,
      file.path(paths$qualification, "convergence.csv")),
    candidate_windows = .benchmark_write_csv(candidates,
      file.path(paths$qualification, "candidate_windows.csv")),
    qualification_decision = decision_path,
    runtime = .benchmark_write_csv(runtime,
      file.path(paths$qualification, "runtime.csv")))
  if (!passed)
    stop("Study 06 qualification failed; final benchmark remains blocked.",
      call. = FALSE)
  list(spec = spec, paths = c(paths, files), status = status,
    qualification_coordinates = coordinates, qualification = decision,
    convergence = convergence, runtime = runtime)
}

.annotation_marker_table <- function(result, bundle, coordinate) {
  effects <- extract_marker_effects(result)
  if (is.matrix(effects)) effects <- effects[, 1L]
  ids <- names(effects)
  if (is.null(ids)) ids <- rownames(extract_marker_effects(result))
  marker_ids <- bundle$simulation$data$marker_ids
  effects <- as.numeric(effects[match(marker_ids, ids)])
  probabilities <- extract_marker_probabilities(result)
  pip <- probabilities$posterior_inclusion
  if (is.matrix(pip)) pip <- pip[, 1L]
  if (is.null(pip))
    stop("Study 06 requires marker posterior inclusion probabilities.",
      call. = FALSE)
  pip_ids <- names(pip)
  if (is.null(pip_ids) && !is.null(rownames(probabilities$posterior_inclusion)))
    pip_ids <- rownames(probabilities$posterior_inclusion)
  pip <- as.numeric(pip[match(marker_ids, pip_ids)])
  if (any(!is.finite(effects)) || any(!is.finite(pip)))
    stop("Study 06 marker output alignment failed.", call. = FALSE)
  data.frame(study = "06_annotation_models",
    scenario = coordinate$scenario, replicate = coordinate$replicate,
    method = coordinate$method, marker_id = marker_ids,
    posterior_mean_effect = effects,
    posterior_inclusion_probability = pip, stringsAsFactors = FALSE)
}

.annotation_parameter_estimates <- function(result, coordinate, bundle) {
  native <- .benchmark_native_fit(result)
  bundle <- native$convergence_traces
  if (is.null(bundle$values)) return(NULL)
  groups <- as.character(bundle$quantities$group)
  index <- match(c("vgs", "ves"), groups)
  if (anyNA(index)) return(NULL)
  vg <- as.numeric(bundle$values[, , index[1L]])
  ve <- as.numeric(bundle$values[, , index[2L]])
  values <- list(genetic_variance = vg, residual_variance = ve,
    heritability = vg / (vg + ve))
  train <- match(bundle$simulation$data$train_ids,
    bundle$simulation$data$sample_ids)
  true_vg <- stats::var(bundle$simulation$truth$genetic_values[train, 1L])
  true_ve <- stats::var(bundle$simulation$truth$residuals[train, 1L])
  truth <- c(genetic_variance = true_vg, residual_variance = true_ve,
    heritability = true_vg / (true_vg + true_ve))
  do.call(rbind, lapply(names(values), function(id) data.frame(
    study = "06_annotation_models", scenario = coordinate$scenario,
    replicate = coordinate$replicate, method = coordinate$method,
    parameter = id, posterior_mean = mean(values[[id]]),
    posterior_sd = stats::sd(values[[id]]),
    lower_95 = unname(stats::quantile(values[[id]], .025)),
    upper_95 = unname(stats::quantile(values[[id]], .975)),
    truth = truth[[id]], status = "ok",
    reason = "", stringsAsFactors = FALSE)))
}

.run_annotation_final <- function(spec, profile, coordinates, paths, resume) {
  decision <- .read_annotation_qualification_decision(spec, paths)
  data <- prepare_prediction_data(spec, paths$root)
  logic <- .annotation_logic(spec)
  bundles <- .annotation_bundle_map(
    logic$prepare_annotation_simulations(spec, profile, data))
  methods <- resolve_benchmark_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")
  dispatch <- getOption("sblrbench.annotation_fit_dispatch",
    fit_annotation_method)
  status_rows <- runtime_rows <- marker_rows <- prior_rows <- list()
  alpha_rows <- metric_rows <- prediction_rows <- parameter_rows <- list()
  convergence_rows <- list()
  for (i in seq_len(nrow(coordinates))) {
    coordinate <- coordinates[i, , drop = FALSE]
    bundle <- bundles[[.prediction_bundle_key(coordinate$scenario,
      coordinate$replicate)]]
    bundle$spec <- spec
    controls <- annotation_method_controls(spec, coordinate, profile,
      mode = "final", qualification_decision = decision)
    fit <- .annotation_reusable_qualification_history(spec, coordinate,
      controls, data, bundle, paths, decision, resume)
    if (is.null(fit)) fit <- .annotation_fit(spec, coordinate,
      methods[[coordinate$method]], controls, data, bundle, paths, "final",
      resume, dispatch)
    ok <- !is.null(fit$result)
    status_rows[[i]] <- data.frame(study = spec$study, mode = "final",
      scenario = coordinate$scenario, replicate = coordinate$replicate,
      method = coordinate$method, status = if (ok) "ok" else "failed",
      reason = fit$reason, reused = fit$reused,
      input_hash = fit$identities$checkpoint_hash,
      stringsAsFactors = FALSE)
    runtime_rows[[i]] <- data.frame(study = spec$study,
      scenario = coordinate$scenario, replicate = coordinate$replicate,
      method = coordinate$method,
      elapsed_seconds = if (ok) extract_runtime(fit$result) else NA_real_,
      status = if (ok) "ok" else "failed", reason = fit$reason,
      reused = fit$reused, stringsAsFactors = FALSE)
    if (ok) {
      fit$result <- predict_prediction_result(fit$result, bundle$simulation,
        bundle$test_simulation, data$scaled$test)
      marker <- .annotation_marker_table(fit$result, bundle, coordinate)
      marker_rows[[i]] <- marker
      metadata <- list(study = spec$study,
        scenario = as.character(coordinate$scenario),
        replicate = as.integer(coordinate$replicate),
        method = as.character(coordinate$method))
      metric_rows[[i]] <- annotation_marker_recovery(marker,
        bundle$marker_truth, metadata)
      prediction_rows[[i]] <- rbind(
        prediction_effect_recovery(bundle$simulation, fit$result),
        prediction_genetic_value_recovery(bundle$test_simulation, fit$result))
      parameter_rows[[i]] <- .annotation_parameter_estimates(fit$result,
        coordinate, bundle)
      metric_rows[[i]] <- rbind(metric_rows[[i]],
        annotation_parameter_recovery(parameter_rows[[i]], metadata))
      if (isTRUE(methods[[coordinate$method]]$annotation_aware)) {
        traces <- extract_annotation_coefficient_traces(fit$result,
          expected_chains = controls$nchains)
        prior <- summarise_drawwise_annotation_prior(traces,
          bundle$annotations, spec$controls$simulation$mixture_var,
          bundle$simulation$extras$true_marker_prior, bundle$marker_truth)
        if (!identical(prior$status, "ok"))
          stop(prior$reason, call. = FALSE)
        prior$marker$study <- spec$study
        prior$marker$scenario <- coordinate$scenario
        prior$marker$replicate <- coordinate$replicate
        prior$marker$method <- coordinate$method
        prior_rows[[i]] <- prior$marker
        alpha_rows[[i]] <- annotation_alpha_recovery(traces,
          bundle$simulation$extras$true_alpha, metadata)
        metric_rows[[i]] <- rbind(metric_rows[[i]],
          annotation_prior_recovery(prior,
            bundle$simulation$extras$true_marker_prior,
            bundle$annotations, bundle$marker_truth, metadata))
        required <- .annotation_required_traces(fit$result, coordinate,
          bundle, controls$nchains)
        hit <- vapply(decision$entries, function(x)
          identical(x$scenario, as.character(coordinate$scenario)) &&
            identical(x$method, as.character(coordinate$method)), logical(1))
        selected <- decision$entries[[which(hit)]]
        z <- required[required$iteration > selected$selected_burnin &
          required$iteration <= selected$selected_burnin +
            selected$selected_retained, , drop = FALSE]
        diagnostic_spec <- spec
        diagnostic_spec$qualification$burnin_candidates <- 0L
        diagnostic_spec$qualification$retained_candidates <-
          as.integer(selected$selected_retained)
        diagnostic_spec$qualification$maximum_history <-
          as.integer(selected$selected_retained)
        z$iteration <- z$iteration - as.integer(selected$selected_burnin)
        convergence_rows[[i]] <- .annotation_candidate_diagnostics(z,
          diagnostic_spec)
      }
    }
  }
  status <- .benchmark_bind_rows(status_rows)
  if (any(status$status != "ok"))
    stop("Study 06 final benchmark has failed or missing coordinates.",
      call. = FALSE)
  marker_results <- .benchmark_bind_rows(marker_rows)
  metrics <- .benchmark_bind_rows(metric_rows)
  paired <- annotation_paired_advantages(metrics)
  truth <- .benchmark_bind_rows(lapply(bundles, .prediction_truth_table))
  annotations <- attr(bundles, "annotations")
  annotation_truth <- attr(bundles, "annotation_truth")
  annotation_truth_table <- do.call(rbind, lapply(names(spec$scenarios),
    function(scenario) {
      x <- annotation_truth[[scenario]]
      data.frame(scenario = scenario,
        annotation = rep(rownames(x), times = ncol(x)),
        stick = rep(colnames(x), each = nrow(x)), alpha = as.numeric(x),
        stringsAsFactors = FALSE)
    }))
  runtime <- .benchmark_bind_rows(runtime_rows)
  files <- list(
    fit_status = .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv")),
    simulation_truth = .benchmark_write_csv(truth,
      file.path(paths$tables, "simulation_truth.csv")),
    annotation_truth = .benchmark_write_csv(annotation_truth_table,
      file.path(paths$tables, "annotation_truth.csv")),
    annotation_prior_summary = .benchmark_write_csv(
      .benchmark_bind_rows(prior_rows),
      file.path(paths$tables, "annotation_prior_summary.csv")),
    marker_results = .benchmark_write_csv(marker_results,
      file.path(paths$tables, "marker_results.csv")),
    parameter_estimates = .benchmark_write_csv(
      .benchmark_bind_rows(parameter_rows),
      file.path(paths$tables, "parameter_estimates.csv")),
    annotation_alpha = .benchmark_write_csv(.benchmark_bind_rows(alpha_rows),
      file.path(paths$tables, "annotation_alpha.csv")),
    annotation_metrics = .benchmark_write_csv(metrics,
      file.path(paths$tables, "annotation_metrics.csv")),
    paired_comparisons = .benchmark_write_csv(paired,
      file.path(paths$tables, "paired_comparisons.csv")),
    prediction_metrics = .benchmark_write_csv(
      .benchmark_bind_rows(prediction_rows),
      file.path(paths$tables, "prediction_metrics.csv")),
    convergence = .benchmark_write_csv(
      .benchmark_bind_rows(convergence_rows),
      file.path(paths$tables, "convergence.csv")),
    runtime = .benchmark_write_csv(runtime,
      file.path(paths$tables, "runtime.csv")))
  manifest <- .prediction_manifest(spec, profile, paths, coordinates, FALSE,
    status, data)
  manifest$qualification_decision <- file.path(paths$root,
    spec$qualification$decision_path)
  .write_prediction_manifest(manifest, paths$manifest)
  writeLines(benchmark_session_information(), paths$session_info)
  parameter_estimates <- .benchmark_bind_rows(parameter_rows)
  annotation_prior_summary <- .benchmark_bind_rows(prior_rows)
  list(spec = spec, paths = c(paths, files), status = status, truth = truth,
    annotation_truth = annotation_truth_table,
    annotation_prior = annotation_prior_summary,
    annotation_prior_summary = annotation_prior_summary,
    marker_results = marker_results,
    parameter_estimates = parameter_estimates, estimates = parameter_estimates,
    alpha = .benchmark_bind_rows(alpha_rows), metrics = metrics,
    paired = paired, prediction_metrics = .benchmark_bind_rows(prediction_rows),
    convergence = .benchmark_bind_rows(convergence_rows), runtime = runtime,
    qualification = decision)
}

#' Run a benchmark
#'
#' Runs the ordinary-R execution path for Study 02 prediction. Large native fit
#' objects remain in science-identity checkpoints and are not returned.
#'
#' @param spec A specification list or a path accepted by
#'   `read_benchmark_spec()`.
#' @param output_dir Local output directory.
#' @param profile Either `"workshop"` or `"benchmark"`.
#' @param resume Reuse checkpoints only after strict identity validation.
#' @param validate_only Validate coordinates, controls, seeds, and provenance
#'   without preparing data or calling a fit function.
#' @param mode Study 06 execution mode: `"validate_only"`,
#'   `"qualification"`, or `"final"`. Other studies support validation-only
#'   and final execution, but not qualification mode.
#' @return A compact result index containing output tables and paths.
#' @export
run_benchmark <- function(spec, output_dir, profile = "benchmark",
                          resume = TRUE, validate_only = FALSE, mode = NULL) {
  if (is.character(spec) && length(spec) == 1L) spec <- read_benchmark_spec(spec)
  validate_benchmark_spec(spec)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume))
    stop("resume must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(validate_only) || length(validate_only) != 1L ||
      is.na(validate_only))
    stop("validate_only must be TRUE or FALSE.", call. = FALSE)
  if (is.null(mode)) mode <- if (isTRUE(validate_only)) "validate_only" else
    "final"
  mode <- match.arg(mode, c("validate_only", "qualification", "final"))
  if (isTRUE(validate_only)) mode <- "validate_only"
  if (!identical(spec$task, "annotation_models") && identical(mode, "qualification"))
    stop("Qualification mode is supported only for Study 06 annotation models.",
      call. = FALSE)
  resolved <- resolve_benchmark_profile(spec, profile)
  coordinates <- if (identical(spec$task, "annotation_models") &&
      identical(mode, "qualification"))
    benchmark_annotation_seeds(spec, profile, mode = "qualification") else
      benchmark_seeds(spec, profile)
  benchmark_assert_package_sha("sblr", spec$packages$sblr$sha)
  if (!identical(as.character(utils::packageVersion("sblr")),
      spec$packages$sblr$version))
    stop("Installed sblr version does not match spec$packages$sblr$version.",
      call. = FALSE)
  paths <- benchmark_output_paths(output_dir)
  .benchmark_create_output_dirs(paths)
  if (isTRUE(validate_only) || identical(mode, "validate_only")) {
    if (identical(spec$task, "convergence")) {
      benchmark_matched_spec(spec)
      benchmark_convergence_design(spec)
    }
    status <- .prediction_validation_status(coordinates, spec$study)
    status_path <- .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv"))
    coordinate_path <- if (spec$task %in% c("convergence", "finemapping",
        "ld_operator", "annotation_models"))
      .benchmark_write_csv(coordinates,
        file.path(paths$tables, "coordinate_grid.csv")) else NULL
    qualification_coordinates <- annotation_summary <- qualification_path <-
      annotation_path <- NULL
    qualification_schema <- NULL
    if (identical(spec$task, "annotation_models")) {
      qualification_coordinates <- benchmark_annotation_seeds(spec, profile,
        mode = "qualification")
      design <- .annotation_validation_design(spec)
      annotation_summary <- design$summary
      qualification_path <- .benchmark_write_csv(qualification_coordinates,
        file.path(paths$qualification, "qualification_coordinates.csv"))
      annotation_path <- .benchmark_write_csv(annotation_summary,
        file.path(paths$tables, "annotation_design_summary.csv"))
      qualification_schema <- annotation_qualification_artifact_schema()
    }
    manifest <- .prediction_manifest(spec, profile, paths, coordinates, TRUE,
      status)
    .write_prediction_manifest(manifest, paths$manifest)
    writeLines(benchmark_session_information(), paths$session_info)
    return(list(spec = spec, paths = c(paths, fit_status = status_path,
      coordinate_grid = coordinate_path,
      qualification_coordinates = qualification_path,
      annotation_design_summary = annotation_path),
      status = status, truth = NULL, estimates = NULL,
      marker_results = NULL, metrics = NULL, convergence = NULL,
      runtime = NULL, oracle = NULL,
      coordinate_grid = if (spec$task %in% c("convergence", "finemapping",
          "ld_operator", "annotation_models"))
        coordinates else NULL,
      qualification_coordinates = qualification_coordinates,
      annotation_summary = annotation_summary,
      qualification_schema = qualification_schema,
      candidate_summary = NULL, recommendations = NULL))
  }

  if (identical(spec$task, "annotation_models")) {
    if (identical(mode, "qualification"))
      return(.run_annotation_qualification(spec, profile, coordinates, paths,
        resume))
    return(.run_annotation_final(spec, profile, coordinates, paths, resume))
  }

  if (identical(spec$task, "convergence"))
    return(.run_convergence(spec, profile, resolved, coordinates, paths,
      resume))

  if (identical(spec$task, "ld_operator"))
    return(.run_ld_operator(spec, profile, coordinates, paths, resume))

  methods <- resolve_benchmark_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")

  if(identical(spec$task,"finemapping"))
    return(.run_finemapping(spec,profile,coordinates,methods,paths,resume))

  if(identical(spec$task,"parameter_estimation"))
    return(.run_parameter_estimation(spec,profile,resolved,coordinates,methods,
      paths,resume))

  data <- prepare_prediction_data(spec, output_dir)
  bundles <- prepare_prediction_simulations(spec, profile, data)
  bundle_keys <- vapply(bundles, function(x) .prediction_bundle_key(
    x$coordinate$scenario, x$coordinate$replicate), character(1))
  if (anyDuplicated(bundle_keys))
    stop("Prediction simulation coordinates are not unique.", call. = FALSE)
  names(bundles) <- bundle_keys
  fit_dispatch <- getOption("sblrbench.fit_dispatch", fit_prediction_method)
  if (!is.function(fit_dispatch))
    stop("The configured prediction fit dispatch must be a function.",
      call. = FALSE)

  status_rows <- runtime_rows <- estimate_rows <- marker_rows <- list()
  metric_rows <- convergence_rows <- list()
  for (i in seq_len(nrow(coordinates))) {
    coordinate <- coordinates[i, , drop = FALSE]
    coordinate_list <- as.list(coordinate)
    bundle <- bundles[[.prediction_bundle_key(coordinate$scenario,
      coordinate$replicate)]]
    method <- methods[[coordinate$method]]
    controls <- benchmark_method_controls(spec, coordinate$method, profile,
      coordinate$fit_seed, coordinate$chain_seeds[[1L]])
    identity <- .prediction_checkpoint_identity(spec, coordinate_list,
      controls, bundle$simulation, data)
    input_hash <- benchmark_hash_object(identity)
    checkpoint <- .prediction_checkpoint_path(paths, coordinate_list)
    result <- NULL
    reused <- FALSE
    reason <- ""
    if (isTRUE(resume) && file.exists(checkpoint)) {
      loaded <- benchmark_load_checkpoint(checkpoint, input_hash,
        validator = function(x) identical(x$package_sha,
          spec$packages$sblr$sha) && identical(x$coordinate,
          coordinate_list) && inherits(x$result, "sblrbench_result"))
      result <- loaded$value$result
      reused <- TRUE
    } else {
      result <- tryCatch(fit_dispatch(method = method, controls = controls,
        simulation = bundle$simulation, stats = bundle$stats,
        glist = data$ld_glist, split = data$split),
        error = function(error) {
          reason <<- conditionMessage(error)
          NULL
        })
      if (!is.null(result)) {
        result <- predict_prediction_result(result, bundle$simulation,
          bundle$test_simulation, data$scaled$test)
        benchmark_atomic_save_rds(list(input_hash = input_hash,
          package_sha = spec$packages$sblr$sha,
          coordinate = coordinate_list, result = result), checkpoint,
          compress = FALSE, temporary_prefix = ".prediction-fit-")
      }
    }
    ok <- !is.null(result)
    status_rows[[i]] <- data.frame(study = spec$study,
      scenario = coordinate$scenario, replicate = coordinate$replicate,
      method = coordinate$method, status = if (ok) "ok" else "failed",
      reason = reason, reused = reused, input_hash = input_hash,
      stringsAsFactors = FALSE)
    runtime_rows[[i]] <- data.frame(study = spec$study,
      scenario = coordinate$scenario, replicate = coordinate$replicate,
      method = coordinate$method,
      elapsed_seconds = if (ok) extract_runtime(result) else NA_real_,
      status = if (ok) "ok" else "failed", reason = reason, reused = reused,
      stringsAsFactors = FALSE)
    if (ok) {
      estimate_rows[[i]] <- prediction_estimate_table(result,
        coordinate$scenario, coordinate$replicate, coordinate$method)
      marker_rows[[i]] <- prediction_marker_table(result,
        coordinate$scenario, coordinate$replicate, coordinate$method)
      metric_rows[[i]] <- prediction_metric_table(bundle$test_simulation,
        result, spec$metrics)
      convergence_rows[[i]] <- .prediction_convergence_row(result,
        coordinate_list)
    }
  }
  status <- .benchmark_bind_rows(status_rows)
  runtime <- .benchmark_bind_rows(runtime_rows)
  estimates <- .benchmark_bind_rows(estimate_rows)
  marker_results <- .benchmark_bind_rows(marker_rows)
  metrics <- .benchmark_bind_rows(metric_rows)
  convergence <- .benchmark_bind_rows(convergence_rows)
  truth <- .benchmark_bind_rows(lapply(bundles, .prediction_truth_table))
  oracle <- .benchmark_oracle_table(bundles, spec$study)
  paired <- if (is.null(metrics)) NULL else prediction_paired_metrics(metrics)
  files <- list(
    fit_status = .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv")),
    simulation_truth = .benchmark_write_csv(truth,
      file.path(paths$tables, "simulation_truth.csv")),
    simulation_oracle = .benchmark_write_csv(oracle,
      file.path(paths$tables, "simulation_oracle.csv")),
    estimates = .benchmark_write_csv(estimates,
      file.path(paths$tables, "estimates.csv")),
    marker_results = .benchmark_write_csv(marker_results,
      file.path(paths$tables, "marker_results.csv")),
    metrics = .benchmark_write_csv(metrics,
      file.path(paths$tables, "metrics.csv")),
    convergence = .benchmark_write_csv(convergence,
      file.path(paths$tables, "convergence.csv")),
    runtime = .benchmark_write_csv(runtime,
      file.path(paths$tables, "runtime.csv")),
    paired = .benchmark_write_csv(paired,
      file.path(paths$tables, "paired_method_differences.csv")))
  manifest <- .prediction_manifest(spec, profile, paths, coordinates, FALSE,
    status, data)
  .write_prediction_manifest(manifest, paths$manifest)
  writeLines(benchmark_session_information(), paths$session_info)
  list(spec = spec, paths = c(paths, files), status = status, truth = truth,
    estimates = estimates, marker_results = marker_results, metrics = metrics,
    convergence = convergence, runtime = runtime, oracle = oracle)
}
