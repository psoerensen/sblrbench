# Study 02 exact prediction workflow. Run from the repository root or any
# descendant. Fitting and checkpoint mechanics remain in run_benchmark().

find_repository_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "studies", "02_prediction", "spec.R")))
      return(path)
    parent <- dirname(path)
    if (identical(parent, path))
      stop("Could not locate the sblrbench repository root.", call. = FALSE)
    path <- parent
  }
}

root <- find_repository_root()
old_working_directory <- setwd(root)
on.exit(setwd(old_working_directory), add = TRUE)
isolated_library <- file.path("results", "local",
  "current_benchmark_refresh", "rlib")
if (dir.exists(isolated_library))
  .libPaths(unique(c(normalizePath(isolated_library, winslash = "/"),
    .libPaths())))

# Setup and provenance -----------------------------------------------------
library(sblr)
library(sblrbench)
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("Study 02 plotting requires ggplot2.", call. = FALSE)

spec_path <- file.path("studies", "02_prediction", "spec.R")
spec <- read_benchmark_spec(spec_path)
validate_benchmark_spec(spec)
profile <- Sys.getenv("SBLR_BENCH_PROFILE", "benchmark")
profile_settings <- resolve_benchmark_profile(spec, profile)
output_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path("results", "local", "02_prediction"))
installed_sblr <- benchmark_package_provenance("sblr")
if (!identical(installed_sblr$sha, spec$packages$sblr$sha))
  stop("Installed sblr SHA does not match the Study 02 specification.",
    call. = FALSE)
print(list(study = spec$study, profile = profile,
  expected_sblr_sha = spec$packages$sblr$sha,
  installed_sblr_version = installed_sblr$version,
  installed_sblr_sha = installed_sblr$sha,
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)))

# Data design --------------------------------------------------------------
data_design <- benchmark_data_summary(spec, profile)
print(data_design)

# Scenario design ----------------------------------------------------------
scenario_design <- benchmark_scenario_table(spec, profile)
print(scenario_design)

# Coordinates and seeds ----------------------------------------------------
coordinate_seeds <- benchmark_coordinate_table(spec, profile)
coordinate_summary <- aggregate(replicate ~ scenario + method,
  coordinate_seeds, function(x) length(unique(x)))
names(coordinate_summary)[3L] <- "replicate_count"
print(coordinate_summary)
print(utils::head(coordinate_seeds, 8L))
cat("Coordinate count:", nrow(coordinate_seeds), "\n")

# Methods, priors, and controls -------------------------------------------
method_design <- benchmark_method_table(spec, profile)
print(method_design)

# Execution ----------------------------------------------------------------
results <- run_benchmark(spec = spec, output_dir = output_dir,
  profile = profile, resume = TRUE, validate_only = FALSE)

# Status and validation ----------------------------------------------------
fit_status <- results$status
expected_key <- with(coordinate_seeds,
  paste(scenario, replicate, method, sep = "::"))
observed_key <- with(fit_status,
  paste(scenario, replicate, method, sep = "::"))
if (nrow(fit_status) != nrow(coordinate_seeds) ||
    anyDuplicated(observed_key) || !setequal(expected_key, observed_key))
  stop("Study 02 fit status does not cover the expected coordinates.",
    call. = FALSE)
if (any(fit_status$status != "ok"))
  stop("Study 02 contains failed required fits; inspect fit_status.",
    call. = FALSE)
fit_status_counts <- as.data.frame(table(status = fit_status$status),
  stringsAsFactors = FALSE)
checkpoint_summary <- aggregate(reused ~ method, fit_status, sum)
runtime_summary <- benchmark_runtime_summary(results$runtime)
simulation_oracle_status <- results$oracle
if (is.null(simulation_oracle_status) ||
    any(simulation_oracle_status$status != "passed"))
  stop("Study 02 simulation-oracle validation is incomplete.", call. = FALSE)
print(fit_status_counts)
print(checkpoint_summary)
print(runtime_summary)
print(simulation_oracle_status)

# Result tables ------------------------------------------------------------
simulation_truth <- results$truth
prediction_metrics <- results$metrics
effect_recovery <- prediction_metrics[prediction_metrics$metric %in%
  c("effect_rmse", "genetic_value_correlation", "genetic_value_rmse"),
  , drop = FALSE]
runtime <- results$runtime
convergence <- results$convergence
marker_summaries <- results$marker_results
prediction_metric_summary_table <- prediction_metric_summary(
  prediction_metrics)
print(prediction_metric_summary_table)
print(utils::head(effect_recovery, 12L))
print(runtime_summary)
if (!is.null(convergence)) print(utils::head(convergence, 12L))
if (!is.null(marker_summaries)) print(utils::head(marker_summaries, 8L))
cat("Truth rows:", nrow(simulation_truth), "\n")

# Standard figures ---------------------------------------------------------
prediction_plot <- plot_prediction_metrics(prediction_metrics)
calibration_plot <- plot_prediction_calibration(prediction_metrics)
effect_recovery_plot <- plot_effect_recovery(prediction_metrics)
runtime_plot <- plot_benchmark_runtime(runtime)
print(prediction_plot)
print(calibration_plot)
print(effect_recovery_plot)
print(runtime_plot)

figure_dir <- results$paths$figures
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(figure_dir, "prediction-performance.png"),
  prediction_plot, width = 10, height = 6, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "prediction-calibration.png"),
  calibration_plot, width = 10, height = 6, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "effect-recovery.png"),
  effect_recovery_plot, width = 10, height = 6, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "runtime.png"), runtime_plot,
  width = 9, height = 5, dpi = 160)

# Output inventory and extension points -----------------------------------
output_inventory <- benchmark_output_inventory(results)
print(output_inventory)
cat("Tables:", results$paths$tables, "\n",
  "Figures:", results$paths$figures, "\n",
  "Manifest:", results$paths$manifest, "\n",
  "Session information:", results$paths$session_info, "\n",
  "Checkpoints:", results$paths$checkpoints, "\n", sep = "")
cat("Study 02 prediction workflow complete for profile `", profile,
  "` (", profile_settings$replicate_count, " replicates).\n", sep = "")

# Extend the audited workflow through spec.R and the task-specific shared
# functions; do not place new simulations, dispatch, or checkpoint logic here.
