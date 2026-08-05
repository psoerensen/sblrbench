#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
resume <- !any(args == "--no-resume")
validate_only <- any(args == "--validate-only")
trace_smoke <- any(args == "--trace-smoke")
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION")))
  stop("Run this script from the sblrbench repository root.", call. = FALSE)
isolated_library <- file.path(root, "results", "local",
  "current_benchmark_refresh", "rlib")
if (!dir.exists(isolated_library))
  stop("The repository's isolated benchmark R library is unavailable.",
    call. = FALSE)
.libPaths(c(isolated_library, .libPaths()))
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "studies/06_annotation_models/power-isolation.R"),
  local = FALSE)
spec <- read_benchmark_spec(file.path(root,
  "studies/06_annotation_models/spec.R"))
profile <- study06_power_profile()
output <- file.path(root, profile$output)
if (trace_smoke) output <- file.path(tempdir(), "study06-trace-smoke")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

glist_cache <- file.path(root, "results/local/06_annotation_models",
  "checkpoints/data/human_glist.rds")
if (!nzchar(Sys.getenv("SBLR_BENCH_GLIST", "")) && file.exists(glist_cache))
  Sys.setenv(SBLR_BENCH_GLIST = glist_cache)

if (!identical(spec$packages$sblr$sha,
    benchmark_package_provenance("sblr")$sha))
  stop("The loaded sblr package does not match the Study 06 pin.", call. = FALSE)
data <- prepare_prediction_data(spec, output)
logic <- .annotation_logic(spec)
annotations <- logic$construct_annotation_design(data$markers$marker_ids, spec)
annotation_truth <- logic$construct_annotation_truth(annotations, spec)
seed_grid <- benchmark_annotation_seeds(spec, "benchmark", "qualification")
seed_row <- seed_grid[seed_grid$scenario == "informative_annotations" &
  seed_grid$replicate == 1L & seed_grid$method == "st_bed_bayesrc", ,
  drop = FALSE]
simulation_coordinate <- as.list(seed_row[1L, c("scenario", "replicate",
  "component_seed", "effect_seed", "residual_seed")])
simulation <- logic$simulate_annotation_architecture(simulation_coordinate,
  data$scaled$all, data$split$train_rows, annotations, annotation_truth, spec)
test_simulation <- subset_sblrbench_simulation_samples(simulation,
  data$split$test_ids)
stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
  spec$data)
bundle <- list(spec = spec, simulation = simulation,
  test_simulation = test_simulation, stats = stats, annotations = annotations,
  annotation_truth = annotation_truth,
  marker_truth = simulation$extras$marker_truth)
trace_markers <- study06_trace_marker_set(bundle$marker_truth)
truth_identity <- study06_power_truth_identity(data, bundle)
truth_hash <- benchmark_hash_object(truth_identity)
shuffle <- study06_shuffle_annotations(annotations, profile$shuffle_seed)
registry <- study06_power_registry(spec)
registry$fit_seed <- as.integer(seed_row$fit_seed)
registry$chain_seeds <- paste(as.integer(seed_row$chain_seeds[[1L]]),
  collapse = ";")
registry$truth_hash <- truth_hash
registry$annotation_hash <- vapply(registry$annotation_treatment, function(x)
  if (x == "none") "none" else if (x == "shuffled")
    shuffle$annotation_hash else benchmark_hash_object(annotations), character(1))

causal <- bundle$marker_truth$true_nonnull
shuffle_enriched <- shuffle$annotations[, "enriched_binary"] == 1
shuffle_audit <- data.frame(seed = profile$shuffle_seed,
  permutation_hash = shuffle$permutation_hash,
  annotation_hash = shuffle$annotation_hash,
  intercept_exact = identical(shuffle$annotations[, "Intercept"],
    annotations[, "Intercept"]),
  marginal_distributions_exact = all(vapply(seq_len(ncol(annotations)),
    function(j) identical(unname(sort(shuffle$annotations[, j])),
      unname(sort(annotations[, j]))), logical(1))),
  maximum_correlation_difference = max(abs(
    stats::cor(shuffle$annotations[, -1L]) -
      stats::cor(annotations[, -1L]))),
  causal_rate_enriched = mean(causal[shuffle_enriched]),
  causal_rate_unenriched = mean(causal[!shuffle_enriched]),
  causal_enrichment_ratio = mean(causal[shuffle_enriched]) /
    mean(causal[!shuffle_enriched]), stringsAsFactors = FALSE)
truth_audit <- data.frame(profile = profile$id, truth_hash = truth_hash,
  individuals = length(data$sample_ids), training = length(data$split$train_ids),
  validation = length(data$split$test_ids), markers = length(data$markers$marker_ids),
  blocks = length(data$block_start), realized_component_0 =
    simulation$extras$component_counts[1L], realized_component_1 =
    simulation$extras$component_counts[2L], realized_component_2 =
    simulation$extras$component_counts[3L], realized_component_3 =
    simulation$extras$component_counts[4L], realized_heritability =
    simulation$extras$realized_heritability,
  marker_order_hash = benchmark_hash_object(data$markers$marker_ids),
  phenotype_hash = benchmark_hash_object(simulation$truth$phenotypes),
  effect_hash = benchmark_hash_object(simulation$truth$effects),
  summary_statistics_hash = benchmark_hash_object(stats),
  block_hash = benchmark_hash_object(list(data$block_start, data$marker_panel)),
  traced_marker_count = trace_markers$marker_count,
  traced_causal_count = trace_markers$causal_count,
  traced_noncausal_count = trace_markers$noncausal_count,
  trace_estimated_extended_gib = trace_markers$estimated_extended_gib,
  complete_genomewide_occupancy = trace_markers$complete_genomewide_occupancy,
  actual_component_trace_marker_count = 0L,
  component_trace_status = "unavailable_native_abort_with_one_marker",
  stringsAsFactors = FALSE)
study06_write_csv(registry, file.path(output, "registry.csv"))
study06_write_csv(truth_audit, file.path(output, "shared_truth_audit.csv"))
study06_write_csv(shuffle_audit, file.path(output, "shuffle_audit.csv"))
if (trace_smoke) {
  methods <- resolve_benchmark_methods(spec)
  names(methods) <- vapply(methods, `[[`, character(1), "id")
  controls <- study06_power_controls(spec, data$markers$marker_ids[1L],
    seed_row$fit_seed, seed_row$chain_seeds[[1L]], FALSE, FALSE)
  controls$nit <- 20L
  result <- fit_annotation_method(methods[["st_bed_bayesr"]], controls,
    simulation, stats, data$ld_glist, data$split, annotations,
    annotation_truth, data$block_start)
  component <- study06_component_trace(result, data$markers$marker_ids[1L])
  if (!identical(dim(component), c(20L, 4L, 1L)))
    stop("Selected-component trace smoke dimensions are invalid.",
      call. = FALSE)
  message("NON-INFERENTIAL TRACE SMOKE PASSED; no result retained.")
  quit(save = "no", status = 0L)
}
if (validate_only) {
  message("Validated eight-fit paired power isolation registry; no fits run.")
  quit(save = "no", status = 0L)
}

methods <- resolve_benchmark_methods(spec)
names(methods) <- vapply(methods, `[[`, character(1), "id")
checkpoints <- file.path(output, "checkpoints")
dir.create(checkpoints, recursive = TRUE, showWarnings = FALSE)
status_rows <- runtime_rows <- convergence_rows <- occupancy_rows <- list()
occupancy_summary_rows <- metric_rows <- component_rows <- variance_rows <- list()
fixed_rows <- alpha_rows <- prior_rows <- list()
pip_chain_rows <- list()
records <- list()

for (i in seq_len(nrow(registry))) {
  row <- registry[i, , drop = FALSE]
  message("Starting ", i, "/", nrow(registry), ": ", row$fit_id)
  method <- methods[[row$method]]
  fit_annotations <- if (row$annotation_treatment == "shuffled")
    shuffle$annotations else annotations
  controls <- study06_power_controls(spec, character(),
    seed_row$fit_seed, seed_row$chain_seeds[[1L]],
    method$annotation_aware, row$update_alpha)
  identity <- list(schema = "sblrbench-study06-power-isolation-v1",
    profile = profile$id, sblrbench_sha = system("git rev-parse HEAD",
      intern = TRUE), sblr_sha = spec$packages$sblr$sha,
    spec_hash = benchmark_annotation_spec_hash(spec), truth_hash = truth_hash,
    fit = as.list(row), method = method, controls = controls,
    annotation_hash = row$annotation_hash,
    trace_policy = profile$trace_policy,
    trace_marker_hash = "unavailable")
  semantic_hash <- benchmark_hash_object(identity)
  checkpoint <- file.path(checkpoints, paste0(row$fit_id, ".rds"))
  reused <- FALSE
  started <- Sys.time()
  if (resume && file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    saved_identity <- saved$identity
    saved_identity$diagnostic_source_hash <- NULL
    if (!identical(saved_identity, identity))
      stop("Diagnostic checkpoint identity mismatch for ", row$fit_id,
        call. = FALSE)
    result <- saved$result
    semantic_hash <- saved$semantic_hash
    reused <- TRUE
  } else {
    result <- fit_annotation_method(method, controls, simulation, stats,
      data$ld_glist, data$split, fit_annotations, annotation_truth,
      data$block_start)
    benchmark_atomic_save_rds(list(identity = identity,
      semantic_hash = semantic_hash, result = result), checkpoint,
      compress = FALSE, temporary_prefix = ".power-isolation-")
  }
  elapsed <- if (reused) 0 else as.numeric(difftime(Sys.time(), started,
    units = "secs"))
  result <- predict_prediction_result(result, simulation, test_simulation,
    data$scaled$test)
  component <- NULL
  occupancy <- NULL
  draws <- study06_power_required_traces(result, component, occupancy,
    fit_annotations, bundle$marker_truth, row, spec)
  selected <- study06_selected_diagnostics(draws, spec)
  probabilities <- extract_marker_probabilities(result)$posterior_inclusion
  if (is.matrix(probabilities)) probabilities <- probabilities[, 1L]
  pip <- as.numeric(probabilities[match(data$markers$marker_ids,
    names(probabilities))])
  if (any(!is.finite(pip))) stop("Full-marker PIPs are unavailable.",
    call. = FALSE)
  effect <- extract_marker_effects(result)[, 1L]
  effect <- effect[match(data$markers$marker_ids, names(effect))]
  prediction <- as.numeric(data$scaled$test %*% effect)
  genetic_truth <- as.numeric(simulation$truth$genetic_values[
    data$split$test_ids, 1L])
  phenotype <- as.numeric(simulation$truth$phenotypes[data$split$test_ids, 1L])
  metrics <- study06_power_metrics(pip, effect, bundle$marker_truth,
    prediction, genetic_truth, phenotype, row$fit_id)
  metric_rows[[i]] <- metrics
  pip_chain_rows[[i]] <- do.call(rbind, lapply(seq_along(
    result$native_fit$chains),
    function(chain) data.frame(fit_id = row$fit_id, chain = chain,
      causal_mean_pip = mean(result$native_fit$chains[[chain]]$dm[
        bundle$marker_truth$marker_id[causal]]),
      noncausal_mean_pip = mean(result$native_fit$chains[[chain]]$dm[
        bundle$marker_truth$marker_id[!causal]]),
      expected_active = sum(result$native_fit$chains[[chain]]$dm),
      stringsAsFactors = FALSE)))
  component_rows[[i]] <- study06_component_recovery(pip,
    bundle$marker_truth, row$fit_id)
  convergence_rows[[i]] <- selected$rows
  occupancy_rows[[i]] <- do.call(rbind, lapply(seq_along(
    result$native_fit$chains), function(chain) {
      state <- result$native_fit$chains[[chain]]$component
      count <- tabulate(as.integer(state) + 1L, nbins = 4L)
      data.frame(fit_id = row$fit_id, chain = chain,
        component_0 = count[1L], component_1 = count[2L],
        component_2 = count[3L], component_3 = count[4L],
        active_count = sum(count[-1L]),
        source = "final_state_diagnostic_not_convergence_trace",
        stringsAsFactors = FALSE)
    }))
  keep <- draws$iteration > selected$selected$burnin &
    draws$iteration <= selected$selected$burnin + selected$selected$retained
  scalar <- draws[keep & draws$quantity %in% c("effect_variance",
    "genetic_variance", "residual_variance", "heritability"), ]
  variance_rows[[i]] <- aggregate(value ~ quantity, scalar,
    function(x) c(mean = mean(x), sd = stats::sd(x),
      lower = unname(stats::quantile(x, .025)),
      upper = unname(stats::quantile(x, .975))))
  variance_value <- variance_rows[[i]]$value
  variance_rows[[i]] <- data.frame(fit_id = row$fit_id,
    quantity = variance_rows[[i]]$quantity,
    mean = variance_value[, "mean"], sd = variance_value[, "sd"],
    lower = variance_value[, "lower"], upper = variance_value[, "upper"],
    row.names = NULL)
  if (isTRUE(row$fixed_true_alpha)) fixed_rows[[i]] <-
    study06_fixed_alpha_audit(result,
      simulation$extras$true_alpha,
      simulation$extras$true_marker_prior,
      spec$controls$simulation$mixture_var, row$fit_id)
  if (isTRUE(row$update_alpha)) {
    alpha <- extract_annotation_coefficient_traces(result, 4L)
    alpha <- alpha[alpha$iteration > selected$selected$burnin &
      alpha$iteration <= selected$selected$burnin +
        selected$selected$retained, ]
    alpha_rows[[i]] <- aggregate(value ~ parameter + annotation + stick,
      alpha, function(x) c(mean = mean(x), sd = stats::sd(x),
        lower = unname(stats::quantile(x, .025)),
        upper = unname(stats::quantile(x, .975))))
    alpha_value <- alpha_rows[[i]]$value
    alpha_rows[[i]] <- data.frame(fit_id = row$fit_id,
      alpha_rows[[i]][c("parameter", "annotation", "stick")],
      mean = alpha_value[, "mean"], sd = alpha_value[, "sd"],
      lower = alpha_value[, "lower"], upper = alpha_value[, "upper"],
      row.names = NULL)
    prior <- draws[keep & grepl("^prior_", draws$quantity), ]
    prior_rows[[i]] <- aggregate(value ~ quantity, prior,
      function(x) c(mean = mean(x), sd = stats::sd(x),
        lower = unname(stats::quantile(x, .025)),
        upper = unname(stats::quantile(x, .975))))
    prior_value <- prior_rows[[i]]$value
    prior_rows[[i]] <- data.frame(fit_id = row$fit_id,
      quantity = prior_rows[[i]]$quantity,
      mean = prior_value[, "mean"], sd = prior_value[, "sd"],
      lower = prior_value[, "lower"], upper = prior_value[, "upper"],
      row.names = NULL)
  }
  status_rows[[i]] <- data.frame(fit_id = row$fit_id, status = "ok",
    reused = reused, semantic_hash = semantic_hash,
    checkpoint_sha256 = digest::digest(file = checkpoint, algo = "sha256"),
    selected_burnin = selected$selected$burnin,
    selected_retained = selected$selected$retained,
    convergence_pass = all(selected$rows$overall_pass), stringsAsFactors = FALSE)
  runtime_rows[[i]] <- data.frame(fit_id = row$fit_id,
    fit_reported_seconds = extract_runtime(result),
    wall_seconds = elapsed, reused = reused, stringsAsFactors = FALSE)
  records[[row$fit_id]] <- list(effect = effect, prediction = prediction,
    pip = pip, variance = variance_rows[[i]],
    mean_active = mean(vapply(result$native_fit$chains, function(x)
      sum(x$dm), numeric(1))))
  message("Completed ", row$fit_id, "; convergence_pass=",
    status_rows[[i]]$convergence_pass, "; reused=", reused)
  rm(result, component, draws); invisible(gc())
}

status <- do.call(rbind, status_rows)
metrics <- do.call(rbind, metric_rows)
convergence <- do.call(rbind, convergence_rows)
convergence_status <- setNames(status$convergence_pass, status$fit_id)
contrasts <- study06_metric_contrasts(metrics, registry, convergence_status)
route_rows <- list()
for (condition in unique(registry$condition)) {
  bed <- registry$fit_id[registry$condition == condition & registry$route == "bed"]
  block <- registry$fit_id[registry$condition == condition &
    registry$route == "block_eigen"]
  vb <- records[[bed]]$variance; ve <- records[[block]]$variance
  value <- function(x, q) x$mean[x$quantity == q]
  route_rows[[condition]] <- data.frame(condition = condition,
    heritability_difference = abs(value(vb, "heritability") -
      value(ve, "heritability")),
    genetic_variance_difference = abs(value(vb, "genetic_variance") -
      value(ve, "genetic_variance")),
    residual_variance_difference = abs(value(vb, "residual_variance") -
      value(ve, "residual_variance")),
    validation_genetic_value_correlation = stats::cor(records[[bed]]$prediction,
      records[[block]]$prediction),
    posterior_effect_correlation = stats::cor(records[[bed]]$effect,
      records[[block]]$effect), pip_correlation = stats::cor(records[[bed]]$pip,
      records[[block]]$pip),
    pip_auprc_difference = abs(metrics$value[metrics$fit_id == bed &
      metrics$metric == "pip_auprc"] - metrics$value[metrics$fit_id == block &
      metrics$metric == "pip_auprc"]),
    expected_active_difference = abs(records[[bed]]$mean_active -
      records[[block]]$mean_active), stringsAsFactors = FALSE)
}
route <- do.call(rbind, route_rows)
baseline <- status$convergence_pass[match(c("baseline--bed",
  "baseline--block_eigen"), status$fit_id)]
fixed <- status$convergence_pass[match(c("fixed_true_alpha--bed",
  "fixed_true_alpha--block_eigen"), status$fit_id)]
learned <- status$convergence_pass[match(c("learned_informative--bed",
  "learned_informative--block_eigen"), status$fit_id)]
interpretation <- if (all(baseline) && all(fixed) && !all(learned)) {
  "alpha hierarchy and feedback dominate"
} else if (all(baseline) && !all(fixed)) {
  "marker-specific-prior allocation mixing dominates"
} else if (!all(baseline)) {
  if (baseline[1L] && !baseline[2L]) "mixed outcome" else
    "general BayesR allocation or variance mixing dominates"
} else if (all(c(baseline[1L], fixed[1L], learned[1L])) &&
    !all(c(baseline[2L], fixed[2L], learned[2L]))) {
  "block-eigen route issue dominates"
} else "mixed outcome"
decision <- list(schema = "sblrbench-study06-power-isolation-v1",
  profile = profile$id, formal_qualification_unchanged = TRUE,
  final_benchmark_authorized = FALSE, spec_hash = benchmark_annotation_spec_hash(spec),
  truth_hash = truth_hash, shuffle_seed = profile$shuffle_seed,
  shuffle_hash = shuffle$annotation_hash, sblrbench_sha = system(
    "git rev-parse HEAD", intern = TRUE), sblr_sha = spec$packages$sblr$sha,
  all_fits_complete = all(status$status == "ok"), interpretation = interpretation,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))

study06_write_csv(status, file.path(output, "fit_status.csv"))
study06_write_csv(do.call(rbind, runtime_rows), file.path(output, "runtime.csv"))
study06_write_csv(convergence, file.path(output, "selected_convergence.csv"))
study06_write_csv(do.call(rbind, occupancy_rows),
  file.path(output, "final_component_states.csv"))
study06_write_csv(data.frame(
  requested_public_control = "selected_marker_quantities=component",
  smallest_smoke_marker_count = 1L, smallest_smoke_iterations = 20L,
  smallest_smoke_chains = 4L, status = "unavailable_native_abort",
  formal_occupancy_convergence_used = FALSE,
  reason = paste("Pinned public BED BayesR aborted before returning when",
    "selected-component tracing was requested; final states are diagnostic only.")),
  file.path(output, "occupancy_trace_availability.csv"))
study06_write_csv(do.call(rbind, variance_rows),
  file.path(output, "variance_summary.csv"))
study06_write_csv(metrics, file.path(output, "causal_identification_metrics.csv"))
study06_write_csv(do.call(rbind, component_rows),
  file.path(output, "component_stratified_recovery.csv"))
study06_write_csv(do.call(rbind, pip_chain_rows),
  file.path(output, "chain_pip_summary.csv"))
study06_write_csv(contrasts, file.path(output, "paired_contrasts.csv"))
study06_write_csv(route, file.path(output, "route_comparisons.csv"))
if (length(Filter(Negate(is.null), fixed_rows))) study06_write_csv(
  do.call(rbind, fixed_rows), file.path(output, "fixed_alpha_audit.csv"))
if (length(Filter(Negate(is.null), alpha_rows))) study06_write_csv(
  do.call(rbind, alpha_rows), file.path(output, "learned_alpha_summary.csv"))
if (length(Filter(Negate(is.null), prior_rows))) study06_write_csv(
  do.call(rbind, prior_rows), file.path(output, "learned_prior_summary.csv"))
jsonlite::write_json(decision, file.path(output, "diagnostic_decision.json"),
  pretty = TRUE, auto_unbox = TRUE)
message("Study 06 paired power isolation complete: ", interpretation)
