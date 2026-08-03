# Study 03 exact parameter-estimation workflow. Low-level preparation,
# simulation, fitting, checkpointing, and extraction remain shared.

find_repository_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "studies",
          "03_parameter_estimation", "spec.R"))) return(path)
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
  stop("Study 03 plotting requires ggplot2.", call. = FALSE)

spec_path <- file.path("studies", "03_parameter_estimation", "spec.R")
spec <- read_benchmark_spec(spec_path)
validate_benchmark_spec(spec)
profile <- Sys.getenv("SBLR_BENCH_PROFILE", "benchmark")
profile_settings <- resolve_benchmark_profile(spec, profile)
output_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path("results", "local", "03_parameter_estimation"))
installed_sblr <- benchmark_package_provenance("sblr")
if (!identical(installed_sblr$sha, spec$packages$sblr$sha))
  stop("Installed sblr SHA does not match the Study 03 specification.",
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

# Estimands ----------------------------------------------------------------
estimand_design <- benchmark_estimand_table(spec)
print(estimand_design)

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
  stop("Study 03 fit status does not cover the expected coordinates.",
    call. = FALSE)
if (any(fit_status$status != "ok"))
  stop("Study 03 contains failed required fits; inspect fit_status.",
    call. = FALSE)
fit_status_counts <- as.data.frame(table(status = fit_status$status),
  stringsAsFactors = FALSE)
checkpoint_summary <- aggregate(reused ~ method, fit_status, sum)
simulation_oracle_status <- results$oracle
if (is.null(simulation_oracle_status) ||
    any(simulation_oracle_status$status != "passed"))
  stop("Study 03 simulation-oracle validation is incomplete.", call. = FALSE)
convergence <- results$convergence
convergence_status <- if (is.null(convergence)) data.frame(
  status = "unavailable", count = 0L) else data.frame(
    true_traces_available = sum(convergence$true_traces_available),
    compact_summaries_available = sum(convergence$compact_summary_available),
    final_states_available = sum(convergence$final_states_available),
    coordinate_count = nrow(convergence))
print(fit_status_counts)
print(checkpoint_summary)
print(simulation_oracle_status)
print(convergence_status)

# Result tables ------------------------------------------------------------
simulation_truth <- results$truth
parameter_estimates <- results$estimates
parameter_metrics <- results$metrics
variance_components <- parameter_estimates[
  parameter_estimates$estimand_id %in% c("effect_variance",
    "total_marker_effect_variance", "genetic_variance",
    "residual_variance", "heritability"), , drop = FALSE]
component_probabilities <- parameter_estimates[
  parameter_estimates$estimand_id == "causal_proportion", , drop = FALSE]
effect_recovery <- NULL # Study 03 defines scalar parameter recovery, not effect metrics.
runtime <- results$runtime
marker_summaries <- results$marker_results
print(parameter_metrics)
print(utils::head(parameter_estimates, 12L))
print(utils::head(variance_components, 12L))
print(utils::head(component_probabilities, 12L))
print(benchmark_runtime_summary(runtime))
if (!is.null(convergence)) print(utils::head(convergence, 12L))
if (!is.null(marker_summaries)) print(utils::head(marker_summaries, 8L))
cat("Simulation-truth rows:", nrow(simulation_truth), "\n")

# Standard figures ---------------------------------------------------------
heritability_recovery_plot <- plot_parameter_recovery(parameter_estimates,
  "heritability")
genetic_variance_plot <- plot_parameter_recovery(parameter_estimates,
  "genetic_variance")
residual_variance_plot <- plot_parameter_recovery(parameter_estimates,
  "residual_variance")
parameter_bias_plot <- plot_parameter_bias(parameter_estimates)
component_probability_plot <- plot_component_probabilities(
  component_probabilities)
runtime_plot <- plot_benchmark_runtime(runtime)
print(heritability_recovery_plot)
print(genetic_variance_plot)
print(residual_variance_plot)
print(parameter_bias_plot)
print(component_probability_plot)
print(runtime_plot)

figure_dir <- results$paths$figures
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(figure_dir, "heritability-recovery.png"),
  heritability_recovery_plot, width = 9, height = 5, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "genetic-variance-recovery.png"),
  genetic_variance_plot, width = 9, height = 5, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "residual-variance-recovery.png"),
  residual_variance_plot, width = 9, height = 5, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "parameter-bias.png"),
  parameter_bias_plot, width = 11, height = 9, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "active-probability.png"),
  component_probability_plot, width = 9, height = 5, dpi = 160)
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
cat("Study 03 parameter-estimation workflow complete for profile `", profile,
  "` (", profile_settings$replicate_count, " replicates).\n", sep = "")

# Extend the audited workflow through spec.R and the task-specific shared
# functions; do not place new estimands, dispatch, or extraction logic here.
