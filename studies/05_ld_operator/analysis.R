# Study 05 exact integrated LD-operator workflow.
# Run from the repository root. Workshop output is structural only and is not
# suitable for method-performance or operator-equivalence claims.

# Setup and provenance -------------------------------------------------------
root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..",
  winslash = "/", mustWork = TRUE)
library(sblr)
library(sblrbench)
spec <- read_benchmark_spec(file.path(root, "studies/05_ld_operator/spec.R"))
source(file.path(root, "studies/05_ld_operator/operator-design.R"))
options(sblrbench.ld_operator_runner = study05_reference_operator_runner)
profile <- Sys.getenv("SBLRBENCH_PROFILE", "benchmark")
output_dir <- Sys.getenv("SBLRBENCH_OUTPUT_DIR",
  file.path(root, "results/local/05_ld_operator"))
installed <- benchmark_package_provenance("sblr")
print(data.frame(
  study = spec$study, profile = profile,
  expected_sblr_sha = spec$packages$sblr$sha,
  installed_sblr_version = installed$version,
  installed_sblr_sha = installed$sha,
  output_dir = output_dir,
  capsule = spec$frozen_capsule,
  stringsAsFactors = FALSE
))

# Data and simulation design ------------------------------------------------
data_design <- benchmark_data_summary(spec, profile)
scenario_design <- benchmark_scenario_table(spec, profile)
coordinates <- benchmark_coordinate_table(spec, profile)
cat("Study 05 coordinates:", nrow(coordinates), "\n")
print(data_design)
print(scenario_design)
print(utils::head(coordinates, 12L))

# Operator design -----------------------------------------------------------
operator_design <- study05_operator_table(spec)
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
  stop("Study 05 required coordinates are missing or failed.", call. = FALSE)
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

# SBayesR sensitivity to LD approximation ----------------------------------
sbayesr_design <- study05_sbayesr_design(spec)
print(sbayesr_design)
sbayesr_fit_status <- results$sbayesr_fit_status
scheduler_comparison <- results$sbayesr_scheduler
sbayesr_variance <- results$sbayesr_variance
conditional_audit <- results$sbayesr_conditionals
quadratic_audit <- results$sbayesr_quadratics
sbayesr_recovery <- results$sbayesr_recovery
sbayesr_eigenvalues <- results$sbayesr_eigenvalues
print(sbayesr_fit_status)
if (any(sbayesr_fit_status$status != "ok"))
  stop("Study 05 SBayesR evidence contains a failed fit.", call. = FALSE)

# Scheduler comparison ------------------------------------------------------
scheduler_plot_data <- scheduler_comparison[
  scheduler_comparison$quantity %in% c("heritability", "vgs", "ves"), ]
scheduler_comparison_plot <- ggplot2::ggplot(scheduler_plot_data,
  ggplot2::aes(label, mean, colour = variant)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~quantity, scales = "free_y") +
  theme_sblrbench() + ggplot2::labs(x = NULL, y = "Posterior mean") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "none")

# Exact, sparse, and block-eigen posterior summaries ------------------------
variance_plot_data <- sbayesr_variance[
  sbayesr_variance$quantity %in% c("heritability", "vgs", "ves"), ]
sbayesr_variance_plot <- ggplot2::ggplot(variance_plot_data,
  ggplot2::aes(label, mean, colour = variant)) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~quantity, scales = "free_y") +
  theme_sblrbench() + ggplot2::labs(x = NULL, y = "Posterior mean") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "none")

# Corrected-score, quadratic, spectral, and recovery audits -----------------
conditional_plot_data <- conditional_audit[
  conditional_audit$effect_state != "all_zero", ]
conditional_score_plot <- ggplot2::ggplot(conditional_plot_data,
  ggplot2::aes(effect_state, max_score_sparse_difference,
    colour = parameter_set)) + ggplot2::geom_point(size = 2) +
  theme_sblrbench() + ggplot2::labs(x = NULL,
    y = "Maximum hard-sparse corrected-score difference") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
quadratic_error_plot <- ggplot2::ggplot(quadratic_audit,
  ggplot2::aes(effect_state, abs(quadratic_error), colour = operator)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08)) +
  ggplot2::scale_y_log10() + theme_sblrbench() +
  ggplot2::labs(x = NULL, y = "Absolute quadratic error (log scale)") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
sbayesr_recovery_plot <- ggplot2::ggplot(sbayesr_recovery,
  ggplot2::aes(label, effect_correlation, colour = variant)) +
  ggplot2::geom_point(size = 2) + theme_sblrbench() +
  ggplot2::labs(x = NULL, y = "Effect correlation") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "none")

print(scheduler_comparison_plot)
print(sbayesr_variance_plot)
print(conditional_score_plot)
print(quadratic_error_plot)
print(sbayesr_recovery_plot)
cat("Integrated capsule:", spec$frozen_capsule, "\n")

# Output inventory ----------------------------------------------------------
print(benchmark_output_inventory(results))
cat("Study 05 LD-operator workflow complete.\n")
