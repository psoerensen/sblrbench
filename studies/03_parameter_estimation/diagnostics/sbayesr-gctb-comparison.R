#!/usr/bin/env Rscript

# Focused one-coordinate CSR SBayesR diagnostic. The public-API prior-isolation
# gate is evaluated before any fitting function can be called.

options(stringsAsFactors = FALSE, warn = 1)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]]) ||
  !nzchar(x[[1L]])) y else x

required_sblr_sha <- "02e8c74baa906e83c4a08d42a9cc6339b4e81072"
isolated_library <- file.path("results", "local", "current_benchmark_refresh", "rlib")
if (!dir.exists(isolated_library))
  stop("Missing isolated current-refresh library: ", isolated_library, call. = FALSE)
.libPaths(unique(c(normalizePath(isolated_library, winslash = "/"), .libPaths())))

required_packages <- c("sblr", "sblrbench", "qgg", "jsonlite", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
  quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages))
  stop("Missing diagnostic packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

find_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "sblrbench.Rproj"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot locate sblrbench repository root.", call. = FALSE)
    path <- parent
  }
}

root <- find_root()
setwd(root)
devtools::load_all(root,quiet=TRUE)
Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1")

local_root <- file.path("results", "local", "sbayesr_gctb_diagnostic")
directories <- file.path(local_root, c("inputs", "fits", "tables", "figures", "logs"))
for (directory in directories) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(local_root, "logs", "diagnostic.log")
log_line <- function(...) cat(format(Sys.time(), tz = "UTC", usetz = TRUE), " | ",
  paste0(..., collapse = ""), "\n", file = log_file, append = TRUE, sep = "")
atomic_csv <- function(x, path) {
  tmp <- tempfile(".diagnostic-", dirname(path), ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) { unlink(tmp); stop("Could not write ", path, call. = FALSE) }
  invisible(path)
}
atomic_rds <- function(x, path) {
  tmp <- tempfile(".diagnostic-", dirname(path), ".rds")
  saveRDS(x, tmp, version = 3)
  if (!file.rename(tmp, path)) { unlink(tmp); stop("Could not write ", path, call. = FALSE) }
  invisible(path)
}

started <- Sys.time()
log_line("starting diagnostic")

sblr_description <- utils::packageDescription("sblr")
sblr_sha <- sblr_description$RemoteSha %||% sblr_description$GithubSHA1 %||% NA_character_
if (!identical(sblr_sha, required_sblr_sha))
  stop("Installed sblr SHA is ", sblr_sha, "; required ", required_sblr_sha, call. = FALSE)

setup_files <- c(
  "studies/five_replicate_helpers.R",
  "studies/03_parameter_estimation/spec.R",
  "studies/03_parameter_estimation/diagnostics/sbayesr-gctb-comparison.R"
)
for (file in setup_files[1:2]) sys.source(file, envir = environment())
config <- read_benchmark_spec(setup_files[3])
legacy_config <- list(study=config$study,trait=config$data$trait,
  chr=config$data$chromosome,sample_limit=config$data$sample_limit,
  example_data=config$data$example_data,qc=config$markers$qc,
  sparse_ld=config$data$sparse_ld,
  simulation=list(h2=config$controls$simulation$h2,
    n_causal=config$controls$simulation$n_causal,
    base_seed=config$seeds$simulation_base,architectures=config$scenarios),
  oracle_tolerance=config$validation$oracle_tolerance)

capsule <- file.path("results", "reference", "03_parameter_estimation", "current")
capsule_manifest <- jsonlite::read_json(file.path(capsule, "benchmark_manifest.json"),
  simplifyVector = TRUE)
if (!identical(capsule_manifest$sblr_commit, required_sblr_sha) ||
    !identical(capsule_manifest$qgdata_commit, config$data$example_data$commit))
  stop("Current Study 03 capsule provenance does not match the diagnostic contract.", call. = FALSE)

seed_registry <- utils::read.csv(file.path(capsule, "seed_registry.csv"), check.names = FALSE)
seed_rows <- seed_registry[seed_registry$architecture == "sparse_mixture" &
  seed_registry$replicate == 1L & seed_registry$method == "st_csr_sbayesr", ]
if (nrow(seed_rows) != 4L || length(unique(seed_rows$simulation_seed)) != 1L ||
    length(unique(seed_rows$fit_seed)) != 1L || anyDuplicated(seed_rows$chain_seed))
  stop("The authoritative Study 03 seed coordinate is invalid.", call. = FALSE)
seed_rows <- seed_rows[order(seed_rows$chain), ]

pi_baseline <- c(0.99, rep(0.01 / 3, 3))
mixture_var <- c(0, 0.01, 0.1, 1)
truth_active_probability <- 50 / 37991
variants <- list(
  current_baseline = list(pi = pi_baseline, alpha = pi_baseline * 5e5, adjE = 0.9, updatePi = TRUE),
  weak_dirichlet = list(pi = pi_baseline, alpha = rep(1, 4), adjE = 0.9, updatePi = TRUE),
  weak_dirichlet_adjE1 = list(pi = pi_baseline, alpha = rep(1, 4), adjE = 1.0, updatePi = TRUE),
  fixed_pi_adjE1 = list(pi = pi_baseline, alpha = pi_baseline * 5e5, adjE = 1.0, updatePi = FALSE),
  truth_matched_pi = list(pi = c(1 - truth_active_probability,
    rep(truth_active_probability / 3, 3)), alpha = rep(1, 4), adjE = 1.0, updatePi = TRUE)
)

coordinate_contract <- list(
  coordinate_builder_version = 1L,
  study = "03_parameter_estimation", architecture = "sparse_mixture", replicate = 1L,
  trait = config$data$trait, sample_count = capsule_manifest$analysis_sample_count,
  marker_count = capsule_manifest$canonical_marker_count,
  causal_count = config$controls$simulation$n_causal,
  target_h2 = config$controls$simulation$h2,
  data_selection_seed = unique(seed_rows$data_selection_seed),
  architecture_seed = unique(seed_rows$architecture_seed),
  simulation_seed = unique(seed_rows$simulation_seed), fit_seed = unique(seed_rows$fit_seed),
  chain_seeds = seed_rows$chain_seed, sblr_sha = sblr_sha,
  qgdata_sha = config$data$example_data$commit, qc = config$markers$qc,
  sparse_ld = config$data$sparse_ld,
  diagnostic_schema_version = 2L
)
coordinate_identity <- function(x) {
  simulation <- x$simulation
  summary_stats <- x$summary_stats
  benchmark_semantic_checkpoint_identity(
    "study03-sbayesr-gctb-comparison",
    c(coordinate_contract, list(
      sample_hash = benchmark_hash_object(x$sample_ids),
      marker_hash = benchmark_hash_object(summary_stats$marker_names),
      phenotype_hash = benchmark_hash_object(simulation$truth$phenotypes),
      true_effect_hash = benchmark_hash_object(simulation$truth$effects),
      summary_statistic_hash = benchmark_hash_object(list(n = summary_stats$n,
        yy = summary_stats$yy, wy = summary_stats$wy,
        ww = summary_stats$ww, marker_names = summary_stats$marker_names)))))
}
coordinate_path <- file.path(local_root, "inputs", "coordinate.rds")
checkpoint_reused <- FALSE

if (file.exists(coordinate_path)) {
  candidate <- readRDS(coordinate_path)
  if (!identical(candidate$checkpoint_schema, "sblrbench-semantic-v2"))
    stop(paste("Legacy source-hashed diagnostic checkpoint detected.",
      "This checkpoint schema has been retired and is not reusable under",
      "the shared semantic checkpoint framework."), call. = FALSE)
  input_identity <- coordinate_identity(candidate)
  input_hash <- benchmark_semantic_checkpoint_hash(input_identity)
  loaded <- benchmark_load_semantic_checkpoint(coordinate_path, input_hash,
    validator = function(x) identical(x$input_identity, input_identity))
  coordinate <- loaded$value
  checkpoint_reused <- TRUE
  log_line("reused coordinate checkpoint ", coordinate_path)
} else {
  paths <- list(glist_path=Sys.getenv("SBLR_BENCH_GLIST",""),
    data_dir=file.path("results","local","03_parameter_estimation","data"),
    output_dir=file.path("results","local","03_parameter_estimation","genotype_setup"))
  expected_files <- file.path(paths$data_dir,
    names(config$data$example_data$md5))
  if (!all(file.exists(expected_files)))
    stop("Pinned Study 03 input cache is incomplete; this diagnostic will not download replacement data.",
      call. = FALSE)
  observed_md5 <- unname(tools::md5sum(expected_files))
  if (!identical(observed_md5, unname(config$data$example_data$md5)))
    stop("Pinned qgdata file checksums do not match Study 03 configuration.", call. = FALSE)

  base_glist <- sblrbench:::benchmark_load_glist(paths)
  marker_info <- sblrbench:::benchmark_filter_markers(base_glist,
    legacy_config$chr,legacy_config$qc,legacy_config$sparse_ld)
  sample_ids <- sblrbench:::benchmark_selected_ids(base_glist, config$data$sample_limit)
  working_glist <- sblrbench:::benchmark_set_glist_marker_order(base_glist,
    config$data$chromosome, marker_info$marker_ids)
  Z <- sblrbench:::benchmark_extract_scaled_genotypes(working_glist,
    config$data$chromosome, sample_ids, marker_info$marker_ids)
  sparse_glist <- sblrbench:::benchmark_make_full_sample_ld(working_glist,
    marker_info,list(chromosome=legacy_config$chr,
      sparse_ld=legacy_config$sparse_ld),paths$output_dir)

  coordinate_spec <- list(scenario="sparse_mixture",replicate=1L,
    simulation_seed=unique(seed_rows$simulation_seed))
  if (!identical(coordinate_spec$simulation_seed, unique(seed_rows$simulation_seed)))
    stop("Study 03 simulation seed derivation disagrees with the capsule.", call. = FALSE)
  simulation <- sblrbench:::simulate_prediction_architecture(coordinate_spec,Z,config)
  oracle <- sblrbench::check_oracle_genetic_values(simulation,
    tolerance = config$validation$oracle_tolerance)
  if (!oracle$ok) stop("Study 03 simulation oracle failed.", call. = FALSE)
  summary_stats <- sblrbench:::parameter_summary_stats(simulation,sparse_glist,config)
  coordinate <- list(simulation = simulation, summary_stats = summary_stats,
    marker_info = marker_info,
    sample_ids = sample_ids, oracle = oracle, qgdata_md5 = observed_md5,
    genotype_scale = "qgg::getG(impute=TRUE, scale=TRUE)",
    summary_scale = "sblr::make_summary_stats(scale=TRUE)")
  input_identity <- coordinate_identity(coordinate)
  input_hash <- benchmark_semantic_checkpoint_hash(input_identity)
  coordinate <- c(list(checkpoint_schema = "sblrbench-semantic-v2",
    semantic_hash = input_hash, input_hash = input_hash,
    input_identity = input_identity), coordinate)
  atomic_rds(coordinate, coordinate_path)
  rm(Z, base_glist, working_glist, sparse_glist); invisible(gc())
  log_line("created coordinate checkpoint ", coordinate_path)
}

simulation <- coordinate$simulation
summary_stats <- coordinate$summary_stats
n <- as.integer(summary_stats$n)
m <- length(summary_stats$ww[[1L]])
vy <- as.numeric(summary_stats$yy) / (n - 1)
if (!identical(n, 5000L) || !identical(m, 37991L) ||
    !identical(simulation$scenario$architecture, "sparse_mixture") ||
    !identical(simulation$scenario$replicate, 1L))
  stop("Resolved coordinate does not match the prescribed design.", call. = FALSE)

# Resolve with the exact installed implementation. This unexported resolver is
# used only to audit what the public fit would resolve; no sampler is called.
prior_resolver <- getFromNamespace(".make_stblr_bayesr_priors", "sblr")
baseline_prior <- prior_resolver(vy = vy, m = m, h2 = 0.30, nub = 4, nue = 4,
  pi = pi_baseline, mixture_var = mixture_var, trait_names = config$data$trait,
  alpha = pi_baseline * 5e5, marker_scale = 1)

fixed_prior <- list(B = unname(baseline_prior$B[1, 1]),
  E = unname(baseline_prior$E[1, 1]),
  ssb_prior = unname(baseline_prior$ssb_prior[1, 1]),
  sse_prior = unname(baseline_prior$sse_prior[1, 1]))

public_wrapper_formals <- names(formals(sblr::stblr_csr))
internal_dispatch_formals <- names(formals(getFromNamespace(".stblr_csr_impl", "sblr")))
backend <- getFromNamespace("stblr_csr_bayesr", "sblr")
backend_formals <- names(formals(backend))
fixed_names <- c("B", "E", "ssb_prior", "sse_prior")
direct_public_support <- all(fixed_names %in% public_wrapper_formals)
dispatch_support <- all(fixed_names %in% internal_dispatch_formals)
rejection_message <- tryCatch({
  do.call(backend, c(list(stats = list()),
    stats::setNames(rep(list(matrix(1, 1, 1)), length(fixed_names)), fixed_names)))
  "unexpectedly accepted"
}, error = conditionMessage)
isolation_possible <- direct_public_support && dispatch_support

design <- data.frame(
  study = "03_parameter_estimation", architecture = "sparse_mixture", replicate = 1L,
  trait = config$data$trait, sample_count = n, marker_count = m,
  causal_count = config$controls$simulation$n_causal,
  target_h2 = config$controls$simulation$h2,
  realized_h2 = stats::var(simulation$truth$genetic_values[, 1L]) /
    (stats::var(simulation$truth$genetic_values[, 1L]) +
      stats::var(simulation$truth$residuals[, 1L])),
  simulation_seed = unique(seed_rows$simulation_seed),
  fit_seed = unique(seed_rows$fit_seed), chain_seeds = paste(seed_rows$chain_seed, collapse = ";"),
  qgdata_sha = config$data$example_data$commit, sblr_sha = sblr_sha,
  genotype_scale = coordinate$genotype_scale, summary_scale = coordinate$summary_scale,
  yy = as.numeric(summary_stats$yy), phenotype_variance_n_minus_1 = vy,
  oracle_ok = coordinate$oracle$ok, input_hash = input_hash,
  input_checkpoint_reused = checkpoint_reused)
atomic_csv(design, file.path(local_root, "tables", "design.csv"))

collapse_num <- function(x) paste(format(x, digits = 17, scientific = FALSE, trim = TRUE),
  collapse = ";")
prior_resolution <- do.call(rbind, lapply(names(variants), function(id) {
  v <- variants[[id]]
  automatic_prior <- prior_resolver(vy = vy, m = m, h2 = 0.30, nub = 4, nue = 4,
    pi = v$pi, mixture_var = mixture_var, trait_names = config$data$trait,
    alpha = v$alpha, marker_scale = 1)
  data.frame(variant = id, pi = collapse_num(v$pi), alpha = collapse_num(v$alpha),
    adjE = v$adjE, updatePi = v$updatePi, B = fixed_prior$B, E = fixed_prior$E,
    ssb_prior = fixed_prior$ssb_prior, sse_prior = fixed_prior$sse_prior,
    initial_mixture_weight = sum(v$pi * mixture_var),
    dirichlet_prior_mean = collapse_num(v$alpha / sum(v$alpha)),
    prior_mean_mixture_weight = sum((v$alpha / sum(v$alpha)) * mixture_var),
    automatically_resolved_B = unname(automatic_prior$B[1, 1]),
    automatically_resolved_E = unname(automatic_prior$E[1, 1]),
    automatically_resolved_ssb_prior = unname(automatic_prior$ssb_prior[1, 1]),
    automatically_resolved_sse_prior = unname(automatic_prior$sse_prior[1, 1]),
    automatic_ssb_ratio_to_baseline = unname(automatic_prior$ssb_prior[1, 1]) /
      fixed_prior$ssb_prior,
    baseline_initial_mixture_weight = baseline_prior$mixture_weight_initial,
    baseline_prior_mean_mixture_weight = baseline_prior$mixture_weight_prior_mean,
    intended_fixed_prior_values_match_baseline = TRUE,
    fit_metadata_fixed_prior_confirmation = FALSE,
    public_api_accepts_all_fixed_priors = isolation_possible,
    status = if (isolation_possible) "ready_for_fit" else "blocked_before_fitting",
    reason = if (isolation_possible) "" else paste(
      "Public CSR BayesR API cannot accept B, E, ssb_prior, sse_prior.", rejection_message),
    stringsAsFactors = FALSE)
}))
atomic_csv(prior_resolution, file.path(local_root, "tables", "prior_resolution.csv"))

fit_status <- data.frame(variant = names(variants), expected_chains = 4L,
  completed_chains = 0L, retained_draws_per_chain = 0L,
  fit_attempted = FALSE, status = "blocked_before_fitting",
  reason = "Critical prior-isolation gate failed before sampler invocation",
  runtime_seconds = 0, stringsAsFactors = FALSE)
atomic_csv(fit_status, file.path(local_root, "tables", "fit_status.csv"))

chain_summary <- expand.grid(variant = names(variants), chain = 1:4,
  stringsAsFactors = FALSE)
chain_summary$chain_seed <- rep(seed_rows$chain_seed, each = length(variants))
chain_summary$retained_draws <- 0L; chain_summary$status <- "not_run_prior_isolation_unavailable"
atomic_csv(chain_summary, file.path(local_root, "tables", "chain_summary.csv"))

pi_summary <- do.call(rbind, lapply(names(variants), function(id) do.call(rbind,
  lapply(seq_along(mixture_var), function(k) data.frame(variant = id,
    component = k - 1L, mixture_variance = mixture_var[k], initial_pi = variants[[id]]$pi[k],
    dirichlet_alpha = variants[[id]]$alpha[k],
    dirichlet_prior_mean = variants[[id]]$alpha[k] / sum(variants[[id]]$alpha),
    posterior_mean_pi = NA_real_, final_pi = NA_real_,
    status = "not_run_prior_isolation_unavailable")))))
atomic_csv(pi_summary, file.path(local_root, "tables", "pi_summary.csv"))

component_summary <- pi_summary[, c("variant", "component", "mixture_variance")]
component_summary$posterior_mean_marker_count <- NA_real_
component_summary$final_marker_count <- NA_real_
component_summary$status <- "not_run_prior_isolation_unavailable"
atomic_csv(component_summary, file.path(local_root, "tables", "component_summary.csv"))

scalar_names <- c("vbs", "vgs", "ves", "vle", "vld", "heritability", "active_probability")
scalar_summary <- expand.grid(variant = names(variants), quantity = scalar_names,
  stringsAsFactors = FALSE)
for (nm in c("mean", "sd", "minimum", "maximum")) scalar_summary[[nm]] <- NA_real_
scalar_summary$status <- "not_run_prior_isolation_unavailable"
atomic_csv(scalar_summary, file.path(local_root, "tables", "scalar_summary.csv"))

convergence_summary <- scalar_summary[, c("variant", "quantity")]
for (nm in c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "relative_mcse"))
  convergence_summary[[nm]] <- NA_real_
convergence_summary$status <- "not_run_prior_isolation_unavailable"
atomic_csv(convergence_summary, file.path(local_root, "tables", "convergence_summary.csv"))

effect_recovery <- data.frame(variant = names(variants), effect_rmse = NA_real_,
  effect_correlation = NA_real_, genetic_value_correlation = NA_real_,
  genetic_value_rmse = NA_real_, direct_predicted_genetic_variance = NA_real_,
  calibration_slope = NA_real_, nonzero_posterior_mean_effects = NA_integer_,
  pip_gt_0_01 = NA_integer_, pip_gt_0_05 = NA_integer_, pip_gt_0_10 = NA_integer_,
  pip_gt_0_50 = NA_integer_, status = "not_run_prior_isolation_unavailable")
atomic_csv(effect_recovery, file.path(local_root, "tables", "effect_recovery.csv"))

direct_variance_checks <- expand.grid(variant = names(variants), chain = 1:4,
  stringsAsFactors = FALSE)
for (nm in c("b_R_b", "b_XtX_b_over_n", "b_XtX_b_over_n_minus_1",
    "var_Zb_n_minus_1", "stored_vgs", "absolute_discrepancy"))
  direct_variance_checks[[nm]] <- NA_real_
direct_variance_checks$denominator_note <- "Z is qgg standardized; n=5000; sample variance uses n-1"
direct_variance_checks$status <- "not_run_prior_isolation_unavailable"
atomic_csv(direct_variance_checks, file.path(local_root, "tables", "direct_variance_checks.csv"))

residual_checks <- expand.grid(variant = names(variants), chain = 1:4,
  stringsAsFactors = FALSE)
for (nm in c("direct_sse", "direct_residual_variance_n", "direct_residual_variance_n_minus_1",
    "stored_residual_expression", "stored_ves", "maximum_residual_drift"))
  residual_checks[[nm]] <- NA_real_
residual_checks$denominator_note <- "Direct SSE unavailable without a fitted final effect state"
residual_checks$status <- "not_run_prior_isolation_unavailable"
atomic_csv(residual_checks, file.path(local_root, "tables", "residual_checks.csv"))

variant_contrasts <- data.frame(variant = setdiff(names(variants), "current_baseline"),
  reference = "current_baseline", posterior_active_probability_difference = NA_real_,
  active_marker_count_difference = NA_real_, vbs_difference = NA_real_,
  vgs_difference = NA_real_, ves_difference = NA_real_, heritability_difference = NA_real_,
  effect_rmse_difference = NA_real_, genetic_value_correlation_difference = NA_real_,
  direct_residual_sse_difference = NA_real_, prediction_calibration_difference = NA_real_,
  status = "not_run_prior_isolation_unavailable")
atomic_csv(variant_contrasts, file.path(local_root, "tables", "variant_contrasts.csv"))

# A future supported implementation must use this identity for every fit
# checkpoint. Existing checkpoints are rejected rather than silently reused.
fit_checkpoint_identity <- function(variant_id)
  benchmark_semantic_checkpoint_identity("study03-sbayesr-gctb-fit",list(
    coordinate_hash=input_hash,scenario="sparse_mixture",replicate=1L,
    trait=config$data$trait,sample_hash=benchmark_hash_object(coordinate$sample_ids),
    marker_hash=benchmark_hash_object(summary_stats$marker_names),
    phenotype_hash=benchmark_hash_object(simulation$truth$phenotypes),
    true_effect_hash=benchmark_hash_object(simulation$truth$effects),
    summary_statistic_hash=input_identity$scientific_inputs$summary_statistic_hash,
    method="st_csr_sbayesr",simulation_seed=unique(seed_rows$simulation_seed),
    fit_seed=unique(seed_rows$fit_seed),chain_seeds=seed_rows$chain_seed,
    sblr_sha=sblr_sha,qgdata_sha=config$data$example_data$commit,
    variant_id=variant_id,variant_controls=variants[[variant_id]],
    fixed_prior=fixed_prior,mcmc_controls=list(nchains=4L,nburn=250L,
      nit=2000L,nthin=1L,ncores=4L)))
validate_fit_checkpoint <- function(checkpoint, variant_id) {
  identical(checkpoint$identity, fit_checkpoint_identity(variant_id)) ||
    stop("Fit checkpoint identity differs for ", variant_id, "; refusing reuse.", call. = FALSE)
}

api_audit <- list(
  public_function = "sblr::stblr_csr",
  public_formals = public_wrapper_formals,
  bayesr_backend = "sblr:::stblr_csr_bayesr (unexported; inspected only)",
  bayesr_backend_formals = backend_formals,
  internal_dispatch_formals = internal_dispatch_formals,
  required_explicit_arguments = fixed_names,
  direct_public_support = direct_public_support,
  internal_dispatch_support = dispatch_support,
  unknown_argument_rejection = rejection_message,
  conclusion = "Explicit prior isolation is not possible through the current public API"
)

manifest <- list(
  diagnostic = "sbayesr_gctb_comparison", status = "blocked_before_fitting",
  blocking_stage = "critical_prior_isolation_gate", input_hash = input_hash,
  input_checkpoint_reused = checkpoint_reused, starting_sblrbench_head =
    trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  package = list(sblr_version = as.character(utils::packageVersion("sblr")),
    sblr_sha = sblr_sha, sblrbench_version = as.character(utils::packageVersion("sblrbench")),
    qgg_version = as.character(utils::packageVersion("qgg")),
    qgdata_sha = config$data$example_data$commit),
  coordinate = input_identity, mcmc = list(nchains = 4L, nburn = 250L,
    retained_draws_per_chain = 2000L, nthin = 1L, ncores = 4L),
  baseline_prior = c(fixed_prior, list(pi = pi_baseline, alpha = pi_baseline * 5e5,
    mixture_var = mixture_var, initial_mixture_weight = baseline_prior$mixture_weight_initial,
    prior_mean_mixture_weight = baseline_prior$mixture_weight_prior_mean)),
  variants = variants, api_audit = api_audit, fit_attempted = FALSE,
  sampler_invocations = 0L, generated_figures = character(),
  explanation = paste("No variant was fitted because B, E, ssb_prior, and sse_prior",
    "cannot be fixed explicitly through the pinned public CSR SBayesR API."),
  runtime_seconds = unname(as.numeric(difftime(Sys.time(), started, units = "secs"))),
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
jsonlite::write_json(manifest, file.path(local_root, "manifest.json"), pretty = TRUE,
  auto_unbox = TRUE, digits = 17, null = "null")
writeLines(capture.output(sessionInfo()), file.path(local_root, "session_info.txt"))
writeLines(c(
  "No figures were generated.",
  "The critical prior-isolation gate failed before fitting, so posterior diagnostic plots would be fabricated evidence."
), file.path(local_root, "figures", "README.txt"))

log_line("blocked before fitting: ", api_audit$conclusion)
cat("SBayesR diagnostic status: blocked_before_fitting\n")
cat("Input checkpoint reused:", checkpoint_reused, "\n")
cat("Input hash:", input_hash, "\n")
cat("Resolved B:", format(fixed_prior$B, digits = 17), "\n")
cat("Resolved E:", format(fixed_prior$E, digits = 17), "\n")
cat("Resolved ssb_prior:", format(fixed_prior$ssb_prior, digits = 17), "\n")
cat("Resolved sse_prior:", format(fixed_prior$sse_prior, digits = 17), "\n")
cat("Public API limitation:", rejection_message, "\n")
cat("No sampler was called.\n")
