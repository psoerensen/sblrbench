# Worked parameter-estimation example
# This compact example uses the same public-data and model contracts as the
# benchmark. Short single-chain settings demonstrate mechanics, not convergence.
library(qgg)
library(sblr)
library(sblrbench)

source(file.path("studies", "01_finemapping", "setup_example_data.R"))
for (f in c("estimands.R", "simulation.R", "methods.R", "metrics.R", "pilot.R"))
  source(file.path("studies", "03_parameter_estimation", f))
config <- source(file.path("studies", "03_parameter_estimation", "config.R"), local = TRUE)$value
config$sample_limit <- 600L
config$profiles$development$nit <- 100L
config$profiles$development$nburn <- 50L
paths <- .study03_paths()
files <- .study01_example_files(paths$data_dir, config$example_data)
Glist <- .study01_load_glist(paths, files)
markers <- .study01_run_qc(Glist, config)
markers$marker_ids <- markers$marker_ids[seq_len(min(600L, length(markers$marker_ids)))]
ids <- .study01_selected_ids(Glist, config$sample_limit)
working <- .study01_set_rsids_ld(Glist, config$chr, markers$marker_ids)
Z <- .study01_extract_genotypes(working, config$chr, ids, markers$marker_ids)
spec <- .study03_replicate_specs(config)[[1L]]
simulation <- .study03_simulate(spec, Z, config)
stopifnot(check_oracle_genetic_values(simulation)$ok)
truth <- .study03_truth(simulation, config)

# Fit BED BayesC, then extract aligned trace draws and nonlinear quantities.
method <- .study03_method_specs(config)[[1L]]
fit <- .study03_fit(method, simulation, stats = NULL, working, config)
stopifnot(fit$status == "ok")
registry <- .study03_estimand_registry()
draws <- .study03_extract_draws(fit$native_fit, method$id, registry, ncol(Z))
posterior <- .study03_summarise_draws(draws)
recovery <- .study03_recovery(posterior, truth, method, registry,
  config$relative_error_tolerance)
print(recovery[, c("estimand_id", "truth", "posterior_mean", "bias",
  "lower_95", "upper_95", "covered_95", "interval_width_95")])
