# Minimal common execution path. Only Study 02 single-trait prediction is
# implemented; unsupported tasks fail before data preparation or fitting.

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
  if (!identical(options[["--study"]], "02_prediction"))
    stop("Unsupported --study; only `02_prediction` is implemented.",
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

.prediction_validation_status <- function(coordinates) {
  data.frame(study = "02_prediction", scenario = coordinates$scenario,
    replicate = coordinates$replicate, method = coordinates$method,
    status = "not_run_validate_only",
    reason = "Specification, profile, coordinates, seeds, methods, and package provenance validated; fitting disabled.",
    stringsAsFactors = FALSE)
}

.prediction_manifest <- function(spec, profile, paths, coordinates,
                                 validate_only, status = NULL,
                                 data = NULL) {
  package <- benchmark_package_provenance("sblr")
  list(schema_version = 1L, study = spec$study, task = spec$task,
    profile = profile, validate_only = isTRUE(validate_only),
    coordinate_count = nrow(coordinates),
    scenarios = names(spec$scenarios), methods = names(spec$methods),
    replicate_count = resolve_benchmark_profile(spec, profile)$replicate_count,
    split = spec$split, marker_policy = spec$markers,
    scientific_controls = spec$controls, seeds = spec$seeds,
    qgdata = spec$data$example_data,
    packages = list(sblr = package, expected = spec$packages),
    data_summary = if (is.null(data)) NULL else list(
      sample_count = length(data$sample_ids),
      training_sample_count = length(data$split$train_ids),
      test_sample_count = length(data$split$test_ids),
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
  methods <- resolve_prediction_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")
  benchmark_assert_package_sha("sblr", spec$packages$sblr$sha)
  if (!identical(as.character(utils::packageVersion("sblr")),
      spec$packages$sblr$version))
    stop("Installed sblr version does not match spec$packages$sblr$version.",
      call. = FALSE)
  paths <- benchmark_output_paths(output_dir)
  .benchmark_create_output_dirs(paths)
  if (isTRUE(validate_only)) {
    status <- .prediction_validation_status(coordinates)
    status_path <- .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv"))
    manifest <- .prediction_manifest(spec, profile, paths, coordinates, TRUE,
      status)
    .write_prediction_manifest(manifest, paths$manifest)
    writeLines(benchmark_session_information(), paths$session_info)
    return(list(spec = spec, paths = c(paths, fit_status = status_path),
      status = status, truth = NULL, estimates = NULL,
      marker_results = NULL, metrics = NULL, convergence = NULL,
      runtime = NULL))
  }

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
    controls <- prediction_method_controls(spec, coordinate$method, profile,
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
  paired <- if (is.null(metrics)) NULL else prediction_paired_metrics(metrics)
  files <- list(
    fit_status = .benchmark_write_csv(status,
      file.path(paths$tables, "fit_status.csv")),
    simulation_truth = .benchmark_write_csv(truth,
      file.path(paths$tables, "simulation_truth.csv")),
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
    convergence = convergence, runtime = runtime)
}
