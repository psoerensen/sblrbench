# Study 04 exact convergence workflow. Data preparation, simulation, fitting,
# checkpointing, true-trace extraction, and diagnostic calculations remain in
# the shared framework.

find_repository_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "studies", "04_convergence", "spec.R")))
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
  stop("Study 04 plotting requires ggplot2.", call. = FALSE)

spec <- read_benchmark_spec(file.path("studies", "04_convergence", "spec.R"))
validate_benchmark_spec(spec)
matched_spec <- sblrbench:::benchmark_matched_spec(spec)
profile <- Sys.getenv("SBLR_BENCH_PROFILE", "benchmark")
profile_settings <- resolve_benchmark_profile(spec, profile)
output_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path("results", "local", "04_convergence"))
installed_sblr <- benchmark_package_provenance("sblr")
if (!identical(installed_sblr$sha, spec$packages$sblr$sha))
  stop("Installed sblr SHA does not match the Study 04 specification.",
    call. = FALSE)
print(list(study = spec$study, profile = profile,
  expected_sblr_version = spec$packages$sblr$version,
  expected_sblr_sha = spec$packages$sblr$sha,
  installed_sblr_version = installed_sblr$version,
  installed_sblr_sha = installed_sblr$sha,
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
  frozen_capsules = spec$frozen_capsules))

# Matched Study 03 design --------------------------------------------------
data_design <- benchmark_data_summary(matched_spec,
  if (profile == "benchmark") "benchmark" else "workshop")
scenario_design <- benchmark_scenario_table(matched_spec,
  if (profile == "benchmark") "benchmark" else "workshop")
scenario_design <- scenario_design[scenario_design$scenario %in%
  spec$matched_grid$scenario, , drop = FALSE]
coordinate_seeds <- benchmark_convergence_seeds(spec, profile)
coordinate_summary <- aggregate(replicate ~ stage + scenario + method,
  coordinate_seeds, function(x) length(unique(x)))
names(coordinate_summary)[4L] <- "replicate_count"
print(data_design)
print(scenario_design)
print(spec$matched_grid)
print(coordinate_summary)
print(utils::head(coordinate_seeds, 12L))
cat("Coordinate count:", nrow(coordinate_seeds), "\n")
cat("Seed rule: Study 03 simulation seeds; Study 04 fit base",
  spec$seeds$fit_base, "with architecture, replicate, method, and",
  spec$seeds$chain_stride, "chain strides.\n")

# Methods, priors, and controls -------------------------------------------
selection_controls <- spec$controls$selection
method_design <- do.call(rbind, lapply(seq_len(nrow(spec$matched_grid)),
  function(i) {
    method_id <- spec$matched_grid$method[[i]]
    method <- matched_spec$methods[[method_id]]
    data.frame(scenario = spec$matched_grid$scenario[[i]], method = method_id,
      label = method$label, interface = method$interface,
      model = method$native_method, prior_class = method$prior_class,
      selection_nit = selection_controls$nit,
      selection_nburn = selection_controls$nburn,
      nthin = selection_controls$nthin,
      nchains = selection_controls$nchains,
      ncores = selection_controls$ncores,
      h2_prior = matched_spec$controls$priors$h2,
      inclusion_probability = if (method$prior_class == "BayesC")
        matched_spec$controls$priors$bayesc_inclusion_probability else
        matched_spec$controls$priors$bayesr_active_probability,
      mixture_multipliers = if (method$prior_class == "BayesR")
        paste(matched_spec$controls$priors$bayesr_mixture_var,
          collapse = "; ") else NA_character_,
      update_flags = "convergence=core; keep_chains=TRUE; keep_traces=TRUE",
      stringsAsFactors = FALSE)
  }))
print(method_design)

# Diagnostic design -------------------------------------------------------
diagnostic_design <- benchmark_convergence_design(spec)
convergence_quantities <- diagnostic_design$quantities
burnin_candidates <- diagnostic_design$candidates$burnin
retained_draw_candidates <- diagnostic_design$candidates$retained
convergence_thresholds <- diagnostic_design$thresholds
recommendation_rules <- diagnostic_design$recommendation_rules
print(convergence_quantities)
print(burnin_candidates)
print(retained_draw_candidates)
print(convergence_thresholds)
print(recommendation_rules)

# Execution ----------------------------------------------------------------
results <- run_benchmark(spec = spec, output_dir = output_dir,
  profile = profile, resume = TRUE, validate_only = FALSE)

# Status and validation ----------------------------------------------------
fit_status <- results$status
expected_key <- with(coordinate_seeds,
  paste(stage, scenario, replicate, method, sep = "::"))
observed_key <- with(fit_status,
  paste(stage, scenario, replicate, method, sep = "::"))
if (nrow(fit_status) != nrow(coordinate_seeds) ||
    anyDuplicated(observed_key) || !setequal(expected_key, observed_key))
  stop("Study 04 fit status does not cover the expected coordinates.",
    call. = FALSE)
if (any(fit_status$status != "ok"))
  stop("Study 04 has failed or missing required coordinates; inspect fit_status.",
    call. = FALSE)
fit_status_counts <- as.data.frame(table(stage = fit_status$stage,
  status = fit_status$status), stringsAsFactors = FALSE)
checkpoint_summary <- aggregate(reused ~ stage + method, fit_status, sum)
print(fit_status_counts)
print(checkpoint_summary)
if (is.null(results$convergence) ||
    any(results$convergence$status == "indeterminate"))
  stop("Required true-chain convergence quantities are unavailable or indeterminate.",
    call. = FALSE)

# Result tables ------------------------------------------------------------
coordinate_grid <- results$coordinate_grid
convergence_summary <- results$convergence
candidate_summary <- results$candidate_summary
recommendations <- results$recommendations
validation_replicates <- results$validation_replicates
validation_summary <- results$validation_summary
burnin_stability <- results$burnin_stability
runtime <- results$runtime
print(utils::head(convergence_summary, 16L))
print(candidate_summary)
print(recommendations)
if (!is.null(validation_replicates)) print(validation_replicates)
if (!is.null(validation_summary)) print(validation_summary)
print(benchmark_runtime_summary(runtime))

# Named figures ------------------------------------------------------------
selection_convergence <- convergence_summary[
  convergence_summary$stage == "selection", , drop = FALSE]
rhat_plot <- plot_convergence_rhat(selection_convergence,
  spec$diagnostics$thresholds$rhat)
ess_plot <- plot_convergence_ess(selection_convergence,
  spec$diagnostics$thresholds$ess_bulk)
mcse_plot <- plot_convergence_mcse(selection_convergence,
  spec$diagnostics$thresholds$relative_mcse)
relative_mcse_plot <- mcse_plot
stability_plot <- plot_convergence_stability(burnin_stability,
  spec$diagnostics$thresholds$standardized_mean_shift)
runtime_plot <- plot_benchmark_runtime(runtime)
print(rhat_plot)
print(ess_plot)
print(mcse_plot)
print(relative_mcse_plot)
print(stability_plot)
print(runtime_plot)

figure_dir <- results$paths$figures
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(figure_dir, "rhat.png"), rhat_plot,
  width = 10, height = 6, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "ess.png"), ess_plot,
  width = 10, height = 8, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "relative-mcse.png"),
  relative_mcse_plot, width = 10, height = 6, dpi = 160)
ggplot2::ggsave(file.path(figure_dir, "burnin-stability.png"), stability_plot,
  width = 10, height = 6, dpi = 160)
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
cat("Study 04 convergence workflow complete for profile `", profile,
  "` with ", nrow(coordinate_grid), " matched coordinates.\n", sep = "")

# Extend the audited workflow through spec.R and the focused shared functions;
# do not place trace extraction, diagnostics, thresholds, or fit dispatch here.
