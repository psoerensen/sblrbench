# Study 06 exact analysis: annotation-informed models

## 1. Repository setup ------------------------------------------------------
find_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        dir.exists(file.path(path, "studies"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot locate the sblrbench root.")
    path <- parent
  }
}

root <- find_root()
library(sblr)
library(sblrbench)
source(file.path(root, "studies/06_annotation_models/spec.R"), local = TRUE)
source(file.path(root, "studies/06_annotation_models/annotation-design.R"),
  local = TRUE)

profile <- Sys.getenv("SBLR_BENCH_PROFILE", "benchmark")
mode <- Sys.getenv("SBLR_BENCH_MODE", "validate_only")
mode <- match.arg(mode, c("validate_only", "qualification", "final"))
output_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path(root, "results/local/06_annotation_models"))
validate_benchmark_spec(spec)

## 2. Package and data provenance -----------------------------------------
provenance <- data.frame(
  item = c("study", "status", "profile", "mode", "expected sblr SHA",
    "installed sblr version", "installed sblr SHA", "qgdata SHA",
    "output directory", "partial capsule"),
  value = c(spec$study, spec$status, profile, mode, spec$packages$sblr$sha,
    as.character(utils::packageVersion("sblr")),
    benchmark_package_provenance("sblr")$sha, spec$data$qgdata_sha,
    normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    spec$frozen_capsule$current_stop),
  stringsAsFactors = FALSE)
print(provenance, row.names = FALSE)

## 3. Data design and train/test split ------------------------------------
data_design <- data.frame(
  source = spec$data$source,
  samples = spec$data$expected_sample_count,
  markers = spec$data$expected_marker_count,
  chromosome = spec$data$chromosome,
  training_samples = spec$split$expected_train_count,
  test_samples = spec$split$expected_test_count,
  split_seed = spec$split$seed,
  alignment = paste(unlist(spec$data$alignment), collapse = "; "),
  scaling = spec$data$preprocessing$genotypes,
  summary_statistics = spec$data$preprocessing$summary_statistics,
  stringsAsFactors = FALSE)
print(data_design, row.names = FALSE)

## 4. Scenario and annotation design --------------------------------------
scenario_design <- do.call(rbind, lapply(names(spec$scenarios), function(id) {
  x <- spec$scenarios[[id]]
  data.frame(scenario = id,
    replicates = spec$supported_profiles[[profile]]$replicate_count,
    target_h2 = spec$controls$simulation$h2,
    expected_active = spec$controls$simulation$target_expected_nonnull,
    effect_distribution = x$effect_distribution,
    annotation_relationship = x$annotation_policy,
    stringsAsFactors = FALSE)
}))
print(scenario_design, row.names = FALSE)

preview_ids <- sprintf("marker_%05d", seq_len(spec$data$expected_marker_count))
annotation_design <- construct_annotation_design(preview_ids, spec)
annotation_truth <- construct_annotation_truth(annotation_design, spec)
annotation_summary <- annotation_design_summary(annotation_design,
  annotation_truth, spec)
print(annotation_summary, row.names = FALSE)
print(annotation_truth$informative_annotations)
print(annotation_truth$uninformative_annotations)

## 5. Coordinates, seeds, methods, priors and controls --------------------
final_coordinates <- benchmark_annotation_seeds(spec, profile, mode = "final")
qualification_coordinates <- benchmark_annotation_seeds(spec, profile,
  mode = "qualification")
print(utils::head(final_coordinates, 12L), row.names = FALSE)
print(qualification_coordinates, row.names = FALSE)

method_design <- do.call(rbind, lapply(names(spec$methods), function(id) {
  x <- spec$methods[[id]]
  data.frame(method = id, interface = x$interface,
  model = if (is.null(x$annotation_model)) x$native_method else
    paste(x$native_method, x$annotation_model, sep = "/"),
  annotation_aware = x$annotation_aware,
  chains = spec$qualification$nchains,
  mixture_variances = paste(spec$controls$priors$bayesr_mixture_var,
    collapse = "/"),
  prior_active_probability = spec$controls$priors$bayesr_active_probability,
  scheduler_or_ld = x$computational_policy,
  stringsAsFactors = FALSE)
}))
print(method_design, row.names = FALSE)

## 6. Qualification design and convergence thresholds --------------------
qualification_design <- list(
  entries = qualification_coordinates,
  maximum_history = spec$qualification$maximum_history,
  candidate_burnins = spec$qualification$candidate_burnins,
  candidate_retained = spec$qualification$candidate_retained,
  required_chains = spec$qualification$required_chains,
  thresholds = spec$qualification$thresholds,
  failure_rule = spec$qualification$failure_rule)
print(qualification_design)

estimands <- data.frame(
  family = c("annotation prior", "marker prioritisation", "prediction",
    "genetic recovery", "parameter recovery", "computation"),
  definition = c(
    "Same draw-wise probit-stick marker component prior for BED and CSR",
    "Marker PIP, causal rank, top-k recovery, precision and recall",
    "Phenotype correlation, NMSE and calibration on held-out samples",
    "Genetic-value correlation and effect recovery",
    "Heritability, variance and component-probability recovery",
    "Convergence-qualified runtime"),
  stringsAsFactors = FALSE)
print(estimands, row.names = FALSE)

## 7. Shared execution -----------------------------------------------------
# validate_only is the safe default. Qualification and final execution require
# an explicit SBLR_BENCH_MODE; final also requires a passing frozen decision.
results <- run_benchmark(spec = spec, output_dir = output_dir,
  profile = profile, resume = TRUE,
  validate_only = identical(mode, "validate_only"), mode = mode)

## 8. Qualification status and final summaries ----------------------------
fit_status <- results$status
simulation_truth <- results$truth
annotation_prior_summary <- results$annotation_prior_summary
marker_results <- results$marker_results
parameter_estimates <- results$estimates
annotation_metrics <- results$metrics
prediction_metrics <- results$prediction_metrics
convergence <- results$convergence
runtime <- results$runtime

print(table(fit_status$status, useNA = "ifany"))
if (identical(mode, "validate_only")) {
  message(nrow(results$qualification_coordinates), " qualification entries; ",
    nrow(results$coordinate_grid), " final coordinates; zero fits requested.")
} else if (any(fit_status$status %in% c("failed", "missing"))) {
  stop("Study 06 contains failed or missing required coordinates.")
}

for (x in list(annotation_prior_summary, marker_results, parameter_estimates,
    annotation_metrics, prediction_metrics, convergence, runtime)) {
  if (is.data.frame(x)) print(utils::head(x, 12L), row.names = FALSE)
}

## 9. Named plots ----------------------------------------------------------
annotation_prior_plot <- causal_pip_plot <- causal_rank_plot <-
  prediction_plot <- genetic_recovery_plot <- parameter_recovery_plot <-
  convergence_plot <- runtime_plot <- NULL

if (is.data.frame(annotation_metrics) && nrow(annotation_metrics)) {
  annotation_prior_plot <- plot_annotation_prior_recovery(annotation_metrics)
  causal_pip_plot <- plot_annotation_marker_recovery(annotation_metrics)
  causal_rank_plot <- causal_pip_plot
  print(annotation_prior_plot)
  print(causal_pip_plot)
}
if (is.data.frame(parameter_estimates) && nrow(parameter_estimates) &&
    all(c("scenario", "method", "parameter", "truth", "posterior_mean",
      "lower_95", "upper_95", "status") %in% names(parameter_estimates))) {
  parameter_recovery_plot <-
    plot_annotation_parameter_recovery(parameter_estimates)
  print(parameter_recovery_plot)
}
if (is.data.frame(prediction_metrics) && nrow(prediction_metrics)) {
  prediction_plot <- plot_prediction_metrics(prediction_metrics,
    metric_ids = "genetic_value_correlation")
  genetic_recovery_plot <- plot_effect_recovery(prediction_metrics)
  print(prediction_plot)
  print(genetic_recovery_plot)
}
if (is.data.frame(convergence) && nrow(convergence)) {
  convergence_plot_data <- transform(convergence,
    burnin_candidate = burnin, retained_draw_candidate = retained)
  convergence_plot <- plot_convergence_rhat(convergence_plot_data,
    threshold = spec$qualification$thresholds$rhat)
  print(convergence_plot)
}
if (is.data.frame(runtime) && nrow(runtime)) {
  runtime_plot <- plot_benchmark_runtime(runtime)
  print(runtime_plot)
}

## 10. Output inventory and extension points ------------------------------
print(benchmark_output_inventory(results), row.names = FALSE)
message("Study 06 ", mode, " workflow complete. Outputs: ", output_dir)

# Extension points: add prespecified annotation scenarios in spec.R, implement
# their deterministic construction in annotation-design.R, and add methods only
# when the pinned sblr interface exposes comparable retained-draw quantities.
