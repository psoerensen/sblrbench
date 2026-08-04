# Study 06 exact retained-low-rank LD-operator workflow.
# Run from the repository root. Workshop output is structural only and is not
# suitable for method-performance or operator-equivalence claims.

# Setup and provenance -------------------------------------------------------
root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..",
  winslash = "/", mustWork = TRUE)
library(sblr)
library(sblrbench)
spec <- read_benchmark_spec(file.path(root, "studies/06_ld_operator/spec.R"))
source(file.path(root, "studies/06_ld_operator/operator-design.R"))
options(sblrbench.ld_operator_runner = study06_reference_operator_runner)
profile <- Sys.getenv("SBLRBENCH_PROFILE", "benchmark")
output_dir <- Sys.getenv("SBLRBENCH_OUTPUT_DIR",
  file.path(root, "results/local/06_ld_operator"))
installed <- benchmark_package_provenance("sblr")
print(data.frame(
  study = spec$study, profile = profile,
  expected_sblr_sha = spec$packages$sblr$sha,
  installed_sblr_version = installed$version,
  installed_sblr_sha = installed$sha,
  output_dir = output_dir,
  main_capsule = spec$frozen_capsules$main,
  supplemental_capsule = spec$frozen_capsules$supplemental,
  stringsAsFactors = FALSE
))

# Data and simulation design ------------------------------------------------
data_design <- benchmark_data_summary(spec, profile)
scenario_design <- benchmark_scenario_table(spec, profile)
coordinates <- benchmark_coordinate_table(spec, profile)
cat("Study 06 coordinates:", nrow(coordinates), "\n")
print(data_design)
print(scenario_design)
print(utils::head(coordinates, 12L))

# Operator design -----------------------------------------------------------
operator_design <- study06_operator_table(spec)
block_design <- data.frame(
  policy = spec$operators$block$policy,
  selected_block_size = spec$operators$block$size,
  sensitivity_sizes = paste(spec$operators$block$sensitivity_sizes,
    collapse = ";"),
  eigen_tolerance = spec$operators$eigen$tolerance,
  operator_contract = spec$operators$contract,
  stringsAsFactors = FALSE
)
tolerance_design <- data.frame(
  tolerance = names(spec$operators$equivalence_tolerances),
  value = unlist(spec$operators$equivalence_tolerances),
  row.names = NULL
)
print(operator_design)
print(block_design)
print(tolerance_design)

# Methods, priors, and convergence controls --------------------------------
methods <- benchmark_method_table(spec, profile)
convergence_design <- data.frame(
  maximum_nit = spec$controls$convergence$maximum_nit,
  candidate_burnin = paste(spec$controls$convergence$candidate_burnin,
    collapse = ";"),
  candidate_retained = paste(spec$controls$convergence$candidate_retained,
    collapse = ";"),
  rhat_threshold = spec$controls$convergence$thresholds$rhat,
  bulk_ess_threshold = spec$controls$convergence$thresholds$ess_bulk,
  tail_ess_threshold = spec$controls$convergence$thresholds$ess_tail,
  relative_mcse_threshold = spec$controls$convergence$thresholds$relative_mcse,
  stringsAsFactors = FALSE
)
print(methods)
print(convergence_design)

# Execution -----------------------------------------------------------------
# The shared runner owns validation, scientific-identity checkpoints,
# extraction, metrics, and predictable local outputs. Operator construction
# remains in operator-design.R. Existing frozen capsules are never rewritten.
results <- run_benchmark(
  spec = spec,
  output_dir = output_dir,
  profile = profile,
  resume = TRUE,
  validate_only = FALSE
)

# Status and result tables --------------------------------------------------
fit_status <- results$status
operator_summary <- results$operator_summary
operator_comparisons <- results$operator_comparisons
eigenvalue_summary <- results$eigenvalue_summary
convergence <- results$convergence
recovery_metrics <- results$recovery_metrics
runtime <- results$runtime

print(table(fit_status$status, useNA = "ifany"))
expected <- nrow(benchmark_coordinates(spec, profile))
if (nrow(fit_status) != expected || any(fit_status$status != "ok"))
  stop("Study 06 required coordinates are missing or failed.", call. = FALSE)
if (nrow(operator_summary)) print(utils::head(operator_summary, 12L))
if (nrow(operator_comparisons)) print(operator_comparisons)
if (nrow(eigenvalue_summary)) print(utils::head(eigenvalue_summary, 12L))
if (nrow(convergence)) print(table(convergence$status, useNA = "ifany"))
if (nrow(recovery_metrics)) print(utils::head(recovery_metrics, 12L))
if (nrow(runtime)) print(benchmark_runtime_summary(runtime))

# Named figures -------------------------------------------------------------
operator_error_plot <- plot_operator_errors(operator_comparisons)
retained_rank_plot <- plot_operator_retained_rank(eigenvalue_summary)
eigenvalue_plot <- plot_operator_spectrum(eigenvalue_summary)
recovery_plot <- plot_operator_recovery(recovery_metrics)
runtime_plot <- plot_benchmark_runtime(runtime)

print(operator_error_plot)
print(retained_rank_plot)
print(eigenvalue_plot)
print(recovery_plot)
print(runtime_plot)

# Supplemental evidence ----------------------------------------------------
supplemental_design <- study06_supplemental_design(spec)
supplemental_report <- file.path(root,
  "studies/06_ld_operator/sbayesr_ld_robustness/report.qmd")
print(supplemental_design)
cat("Supplemental SBayesR report:", supplemental_report, "\n")
cat("Supplemental capsule:", spec$frozen_capsules$supplemental, "\n")

# Output inventory ----------------------------------------------------------
print(benchmark_output_inventory(results))
cat("Study 06 LD-operator workflow complete.\n")
