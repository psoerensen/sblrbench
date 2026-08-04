# Study 01 exact fine-mapping workflow
# Run from the repository root. Set SBLRBENCH_PROFILE=workshop for a one-
# replicate structural run; workshop output is not evidence for performance.

# Setup and provenance -------------------------------------------------------
root <- normalizePath(if (file.exists("DESCRIPTION")) "." else "../..",
  winslash = "/", mustWork = TRUE)
library(sblr)
library(sblrbench)
spec <- read_benchmark_spec(file.path(root,"studies/01_finemapping/spec.R"))
profile <- Sys.getenv("SBLRBENCH_PROFILE","benchmark")
output_dir <- Sys.getenv("SBLRBENCH_OUTPUT_DIR",
  file.path(root,"results/local/01_finemapping"))
installed <- benchmark_package_provenance("sblr")
print(data.frame(study=spec$study,profile=profile,
  expected_sblr_sha=spec$packages$sblr$sha,
  installed_sblr_version=installed$version,installed_sblr_sha=installed$sha,
  output_dir=output_dir,frozen_capsule=spec$frozen_capsule))

# Data, scenario, and causal design -----------------------------------------
data_design <- benchmark_data_summary(spec,profile)
scenario_design <- benchmark_scenario_table(spec,profile)
causal_design <- data.frame(
  causal_markers=spec$controls$simulation$n_causal,
  minimum_distance_bp=spec$causal_design$min_distance_bp,
  minimum_maf=spec$causal_design$min_maf,
  maximum_maf=spec$causal_design$max_maf,
  selection_rule=spec$causal_design$selection)
print(data_design); print(scenario_design); print(causal_design)

# Locus and credible-set design ---------------------------------------------
locus_design <- data.frame(
  algorithm=spec$locus_design$algorithm,
  target_probability=spec$locus_design$credible_set_target,
  minimum_r_squared=spec$locus_design$min_r2,
  marker_pip_cutoff=spec$locus_design$pip_cutoff,
  locus_pip_cutoff=spec$locus_design$locus_pip_cutoff,
  maximum_locus_distance_bp=spec$locus_design$max_locus_distance)
print(locus_design)

# Coordinates, seeds, methods, priors, and controls -------------------------
coordinates <- benchmark_coordinate_table(spec,profile)
cat("Coordinates:",nrow(coordinates),"\n")
print(utils::head(coordinates,12L))
methods <- benchmark_method_table(spec,profile)
print(methods)

# Estimands and metrics ------------------------------------------------------
estimands <- spec$estimands
metrics <- data.frame(metric=spec$metrics,stringsAsFactors=FALSE)
print(estimands); print(metrics)

# Execution -----------------------------------------------------------------
results <- run_benchmark(spec=spec,output_dir=output_dir,profile=profile,
  resume=TRUE,validate_only=FALSE)

# Status and validation -----------------------------------------------------
fit_status <- results$status
print(table(fit_status$status,useNA="ifany"))
expected <- nrow(benchmark_coordinates(spec,profile))
if(nrow(fit_status)!=expected || any(fit_status$status!="ok"))
  stop("Study 01 has missing or failed required coordinates.",call.=FALSE)
if(!is.null(results$oracle) && any(results$oracle$status!="passed"))
  stop("Study 01 simulation oracle failed.",call.=FALSE)
cat("Reused checkpoints:",sum(fit_status$reused),"of",expected,"\n")

# Result tables -------------------------------------------------------------
simulation_truth <- results$truth
loci <- results$loci
marker_results <- results$marker_results
credible_sets <- results$credible_sets
finemapping_metrics <- results$metrics
convergence <- results$convergence
runtime <- results$runtime
print(utils::head(loci,10L))
print(aggregate(value ~ method + metric,finemapping_metrics,mean))
print(benchmark_runtime_summary(runtime))

# Standard named figures ----------------------------------------------------
causal_pip_plot <- plot_causal_marker_pip(marker_results)
causal_rank_plot <- plot_causal_marker_rank(marker_results)
credible_set_size_plot <- plot_credible_set_size(credible_sets)
credible_set_coverage_plot <- plot_credible_set_coverage(credible_sets)
pip_calibration_plot <- plot_pip_calibration(marker_results)
runtime_plot <- plot_benchmark_runtime(runtime)
print(causal_pip_plot); print(causal_rank_plot)
print(credible_set_size_plot); print(credible_set_coverage_plot)
print(pip_calibration_plot); print(runtime_plot)

# Output inventory and extension points ------------------------------------
output_inventory <- benchmark_output_inventory(results)
print(output_inventory)
cat("Study 01 complete. Extend scenarios, loci, methods, or metrics in spec.R;",
  "keep fitting and checkpoint mechanics in the shared framework.\n")
