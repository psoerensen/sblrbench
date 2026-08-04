# Minimal common execution path for migrated Studies 02--04.

benchmark_output_paths <- function(output_dir) {
  .benchmark_scalar_string(output_dir, "output_dir")
  list(root = output_dir, checkpoints = file.path(output_dir, "checkpoints"),
    tables = file.path(output_dir, "tables"),
    figures = file.path(output_dir, "figures"),
    manifest = file.path(output_dir, "manifest.json"),
    session_info = file.path(output_dir, "session_info.txt"))
}

parse_benchmark_cli_arguments <- function(args) {
  allowed <- c("--study", "--profile", "--output-dir", "--resume",
    "--validate-only")
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
  supported_studies <- c("01_finemapping", "02_prediction", "03_parameter_estimation",
    "04_convergence", "05_ld_operator")
  if (identical(options[["--study"]], "06_ld_operator"))
    stop("Study ID `06_ld_operator` is retired; use `05_ld_operator`.",
      call. = FALSE)
  if (identical(options[["--study"]], "05_annotation_models"))
    stop("Study ID `05_annotation_models` is retired; use ",
      "`06_annotation_models` (in development).", call. = FALSE)
  if (identical(options[["--study"]], "06_annotation_models"))
    stop("Study `06_annotation_models` is in development and is not yet ",
      "supported by run_benchmark().", call. = FALSE)
  if (!options[["--study"]] %in% supported_studies)
    stop("Unsupported --study; choose `01_finemapping`, `02_prediction`, `03_parameter_estimation`, `04_convergence`, or `05_ld_operator`.",
      call. = FALSE)
  if (!options[["--profile"]] %in% c("workshop", "benchmark"))
    stop("--profile must be `workshop` or `benchmark`.", call. = FALSE)
  parse_boolean <- function(value, option) {
    normalized <- tolower(value)
    if (!normalized %in% c("true", "false"))
      stop(option, " must be `true` or `false`.", call. = FALSE)
    identical(normalized, "true")
  }
  list(study = options[["--study"]], profile = options[["--profile"]],
    output_dir = options[["--output-dir"]],
    resume = parse_boolean(if (is.null(options[["--resume"]])) "true" else
      options[["--resume"]], "--resume"),
    validate_only = parse_boolean(
      if (is.null(options[["--validate-only"]])) "false" else
        options[["--validate-only"]], "--validate-only"))
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
  for (path in paths[c("root", "checkpoints", "tables", "figures")])
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
#' @return A compact result index containing output tables and paths.
#' @export
run_benchmark <- function(spec, output_dir, profile = "benchmark",
                          resume = TRUE, validate_only = FALSE) {
  if (is.character(spec) && length(spec) == 1L) spec <- read_benchmark_spec(spec)
  validate_benchmark_spec(spec)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume))
    stop("resume must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(validate_only) || length(validate_only) != 1L ||
      is.na(validate_only))
    stop("validate_only must be TRUE or FALSE.", call. = FALSE)
  resolved <- resolve_benchmark_profile(spec, profile)
  coordinates <- benchmark_seeds(spec, profile)
  benchmark_assert_package_sha("sblr", spec$packages$sblr$sha)
  if (!identical(as.character(utils::packageVersion("sblr")),
      spec$packages$sblr$version))
    stop("Installed sblr version does not match spec$packages$sblr$version.",
      call. = FALSE)
  paths <- benchmark_output_paths(output_dir)
  .benchmark_create_output_dirs(paths)
  if (isTRUE(validate_only)) {
    if (identical(spec$task, "convergence")) {
      benchmark_matched_spec(spec)
      benchmark_convergence_design(spec)
    }
    status <- .prediction_validation_status(coordinates, spec$study)
    status_path <- .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv"))
    coordinate_path <- if (spec$task %in% c("convergence", "finemapping",
        "ld_operator"))
      .benchmark_write_csv(coordinates,
        file.path(paths$tables, "coordinate_grid.csv")) else NULL
    manifest <- .prediction_manifest(spec, profile, paths, coordinates, TRUE,
      status)
    .write_prediction_manifest(manifest, paths$manifest)
    writeLines(benchmark_session_information(), paths$session_info)
    return(list(spec = spec, paths = c(paths, fit_status = status_path,
      coordinate_grid = coordinate_path),
      status = status, truth = NULL, estimates = NULL,
      marker_results = NULL, metrics = NULL, convergence = NULL,
      runtime = NULL, oracle = NULL,
      coordinate_grid = if (spec$task %in% c("convergence", "finemapping",
          "ld_operator"))
        coordinates else NULL,
      candidate_summary = NULL, recommendations = NULL))
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
