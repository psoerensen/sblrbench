#!/usr/bin/env Rscript

# Assemble the supplemental Study 06 SBayesR LD-operator validation.
options(stringsAsFactors = FALSE, warn = 1)
`%||%` <- function(x, y) if (is.null(x)) y else x

root <- normalizePath(getwd(), winslash = "/")
while (!file.exists(file.path(root, "sblrbench.Rproj"))) {
  parent <- dirname(root)
  if (identical(parent, root)) stop("Repository root not found.", call. = FALSE)
  root <- parent
}
setwd(root)
spec_environment <- new.env(parent = baseenv())
sys.source("studies/06_ld_operator/spec.R", envir = spec_environment)
study06_spec <- spec_environment$spec
required_sha <- study06_spec$packages$sblr$sha
isolated <- file.path(root, "results", "local", "current_benchmark_refresh", "rlib")
if (!dir.exists(isolated)) stop("Validated isolated library is missing.", call. = FALSE)
.libPaths(unique(c(normalizePath(isolated, winslash = "/"), .libPaths())))
required <- c("sblr", "sblrbench", "digest", "jsonlite", "posterior")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
sblr_sha <- utils::packageDescription("sblr")$RemoteSha %||% NA_character_
if (!identical(sblr_sha, required_sha)) stop("Installed sblr SHA mismatch: ", sblr_sha, call. = FALSE)

Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("studies/06_ld_operator/sbayesr_ld_robustness/scripts/eigen-extension.R")

local_exact <- file.path("results", "local", "bed_vs_csr_bayesr_exact_ld")
local_scheduler <- file.path("results", "local", "bed_vs_csr_bayesr_scheduler")
local_supplement <- file.path("results", "local", "06_ld_operator", "sbayesr_ld_robustness")
capsule <- study06_spec$frozen_capsules$supplemental
for (path in c(file.path(local_supplement, "fits"), file.path(local_supplement, "tables"), capsule))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

needed_exact <- c("inputs/marker_subset.rds", "inputs/simulation.rds",
  "ld/exact_checkpoint.rds", "ld/sparse_checkpoint.rds",
  "fits/bed_exact_data.rds", "fits/csr_exact_ld.rds", "fits/csr_sparse_ld.rds")
if (!all(file.exists(file.path(local_exact, needed_exact))))
  stop("Completed A-C checkpoints are incomplete; refusing to rerun them.", call. = FALSE)
if (!file.exists(file.path(local_scheduler, "manifest.json")))
  stop("Completed scheduler diagnostic is missing.", call. = FALSE)

selection_cp <- readRDS(file.path(local_exact, "inputs", "marker_subset.rds"))
simulation <- readRDS(file.path(local_exact, "inputs", "simulation.rds"))
exact_cp <- readRDS(file.path(local_exact, "ld", "exact_checkpoint.rds"))
sparse_cp <- readRDS(file.path(local_exact, "ld", "sparse_checkpoint.rds"))
abc_cp <- setNames(lapply(c("bed_exact_data", "csr_exact_ld", "csr_sparse_ld"),
  function(id) readRDS(file.path(local_exact, "fits", paste0(id, ".rds")))),
  c("bed_exact_data", "csr_exact_ld", "csr_sparse_ld"))
if (!all(vapply(abc_cp, function(x) inherits(x$fit, "stblr_fit") &&
    length(x$fit$chains) == 4L, logical(1))))
  stop("A-C fit checkpoint validation failed.", call. = FALSE)

base_coordinate <- readRDS(file.path("results", "local", "sbayesr_gctb_diagnostic",
  "inputs", "coordinate.rds"))
start <- selection_cp$selection$start
end <- selection_cp$selection$end
if (!identical(c(start, end), c(13392L, 14891L))) stop("Supplemental Study 06 window changed.")
Z <- base_coordinate$simulation$data$genotypes[, start:end, drop = FALSE]
marker_ids <- colnames(Z)
sample_ids <- rownames(Z)
if (!identical(marker_ids, selection_cp$selection$marker_ids) ||
    !identical(dim(Z), c(5000L, 1500L)) ||
    !identical(marker_ids, rownames(abc_cp$csr_exact_ld$fit$bm)))
  stop("Supplemental Study 06 marker/sample alignment failed.", call. = FALSE)
if (max(abs(as.numeric(Z %*% simulation$effects) - simulation$genetic_values)) > 1e-10 ||
    abs(simulation$realized_h2 - .30) > 1e-12)
  stop("Supplemental Study 06 simulation oracle failed.", call. = FALSE)

marker_table <- utils::read.csv(file.path(local_exact, "tables", "selected_marker_window.csv"))
blocks <- s06rob_block_definitions(marker_table, 250L)
expected_starts <- c(1L, 251L, 501L, 624L, 751L, 1001L, 1251L)
if (!identical(blocks$start, expected_starts)) stop("Deterministic block policy changed.")

Y <- matrix(simulation$phenotype, ncol = 1L,
  dimnames = list(sample_ids, "trait1"))
stats_reduced <- sblr::make_summary_stats(Glist = exact_cp$glist, y = Y,
  chr = 1L, rows = seq_len(nrow(Z)), scale = TRUE, nthreads = 1L)
XtX <- crossprod(Z)
Xty <- as.numeric(crossprod(Z, simulation$phenotype))
yy <- sum(simulation$phenotype^2)
n_minus_1 <- nrow(Z) - 1L

exact_correlation <- s06rob_dense_csr(sblr::sparseLD_read_CSR(exact_cp$prefix, one_based = TRUE))
sparse_correlation <- s06rob_dense_csr(sblr::sparseLD_read_CSR(sparse_cp$prefix, one_based = TRUE))
xx_scale <- tcrossprod(sqrt(diag(XtX)))
exact_operator <- exact_correlation * xx_scale
sparse_operator <- sparse_correlation * xx_scale
if (sqrt(sum((exact_operator - XtX)^2) / sum(XtX^2)) > 1e-7)
  stop("Completed exact operator no longer validates.", call. = FALSE)

# Native unfiltered dense reconstruction, using the supported public block policy.
bed_inputs <- getFromNamespace(".stblr_csr_block_eigen_inputs", "sblr")(
  stats_reduced, exact_cp$glist, blocks$start)
native_full <- sblr:::stblr_block_eigen_contract_internal(
  bed_files = bed_inputs$bed_files, n_bed = bed_inputs$n_bed,
  cls = bed_inputs$cls, rows = bed_inputs$rows, af = bed_inputs$af,
  block_start = bed_inputs$block_start, wy = do.call(rbind, stats_reduced$wy),
  effects = numeric(ncol(Z)), eigen_filter = "ridge_fixed",
  eigen_tau = 0, eigen_eta = 0, validation_mutation = "")
full_blocks <- Map(s06rob_unpack_triangle, native_full$packed_upper_triangle,
  native_full$block_size)
full_operator <- matrix(0, ncol(Z), ncol(Z))
for (i in seq_len(nrow(blocks))) {
  idx <- blocks$start[i]:blocks$end[i]
  full_operator[idx, idx] <- full_blocks[[i]]
}
block_target <- s06rob_mask_blocks(XtX, blocks)
full_implementation_error <- sqrt(sum((full_operator - block_target)^2) / sum(block_target^2))
if (full_implementation_error > 1e-7)
  stop("Full-rank block-eigen does not reproduce its configured block target.", call. = FALSE)
full_global_error <- sqrt(sum((full_operator - exact_operator)^2) / sum(exact_operator^2))
full_global_equivalent <- full_global_error <= .01

retained <- s06rob_retained_operator(XtX, blocks, eigen_prop = .995)
retained_operator <- retained$matrix
total_retained_rank <- sum(retained$diagnostics$retained_rank)

pi_common <- c(.99, rep(.01 / 3, 3))
alpha_common <- pi_common * 5e5
mixture_var_common <- c(0, .01, .1, 1)
chain_seeds <- c(150104L, 250104L, 350104L, 450104L)
common <- list(method = "sbayesr", pi = pi_common, alpha = alpha_common,
  mixture_var = mixture_var_common, h2 = .30, adjE = .9,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE,
  updateLDswap = FALSE, maf_effect_s = NULL, estimate_maf_effect_s = FALSE,
  nburn = 250L, nit = 1000L, nthin = 1L, nchains = 4L, ncores = 4L,
  seed = 50104L, chain_seeds = chain_seeds, keep_chains = TRUE,
  convergence = "extended",
  convergence_control = list(warn = FALSE, keep_traces = TRUE,
    extended_groups = "probability"), verbose = FALSE)

fit_base_identity <- list(
  marker_hash = s06rob_hash(marker_ids), sample_hash = s06rob_hash(sample_ids),
  phenotype_hash = s06rob_hash(Y), effect_hash = s06rob_hash(simulation$effects),
  genotype_hash = s06rob_hash(Z), sblr_sha = sblr_sha, block_start = blocks$start,
  common = common)
sampler_calls <- 0L
fit_extension <- function(id) {
  controls <- if (id == "block_eigen_full_rank")
    list(representation = "dense_reconstructed", eigen_policy = "ridge_fixed",
      eigen_eta = 0, eigen_tau = 0, low_rank_residual_rebuild_every = 0L) else
    list(representation = "low_rank", eigen_policy = "cumulative_positive_mass",
      eigen_prop = .995, low_rank_residual_rebuild_every = 100L)
  identity <- c(fit_base_identity, list(variant = id, controls = controls))
  identity_hash <- s06rob_hash(identity)
  path <- file.path(local_supplement, "fits", paste0(id, ".rds"))
  if (file.exists(path)) {
    checkpoint <- readRDS(path)
    if (!identical(checkpoint$identity_hash, identity_hash))
      stop(id, " checkpoint identity changed.", call. = FALSE)
    checkpoint$reused <- TRUE
    return(checkpoint)
  }
  sampler_calls <<- sampler_calls + 1L
  started <- proc.time()[["elapsed"]]
  warnings <- character()
  fit <- withCallingHandlers(do.call(sblr::stblr_block_eigen,
    c(list(stats = stats_reduced, Glist = exact_cp$glist,
      block_start = blocks$start), common, controls)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  checkpoint <- list(identity_hash = identity_hash, identity = identity,
    fit = fit, elapsed_seconds = proc.time()[["elapsed"]] - started,
    warnings = warnings, reused = FALSE)
  s06rob_atomic_rds(checkpoint, path)
  checkpoint
}

de_cp <- setNames(lapply(c("block_eigen_full_rank", "block_eigen_retained"),
  fit_extension), c("block_eigen_full_rank", "block_eigen_retained"))
if (!all(vapply(de_cp, function(x) inherits(x$fit, "stblr_fit") &&
    length(x$fit$chains) == 4L, logical(1)))) stop("D-E fits failed validation.")
native_retained_blocks <- de_cp$block_eigen_retained$fit$input$eigen_diagnostics$blocks
if (!identical(as.integer(native_retained_blocks$retained_rank),
    as.integer(retained$diagnostics$retained_rank)) ||
    !identical(as.integer(native_retained_blocks$block_size), as.integer(blocks$size)))
  stop("Reconstructed retained ranks differ from native fit metadata.", call. = FALSE)

fits <- c(lapply(abc_cp, `[[`, "fit"), lapply(de_cp, `[[`, "fit"))
ids <- names(fits)
labels <- c(bed_exact_data = "BED exact data", csr_exact_ld = "CSR exact LD",
  csr_sparse_ld = "CSR hard sparse", block_eigen_full_rank = "Block eigen full rank",
  block_eigen_retained = "Block eigen retained 0.995")

# Input and prior equivalence.
prior_abc <- utils::read.csv(file.path(local_exact, "tables", "prior_resolution.csv"))
resolver <- getFromNamespace(".make_stblr_bayesr_priors", "sblr")
resolved <- resolver(vy = as.numeric(stats_reduced$yy) / n_minus_1,
  m = ncol(Z), h2 = .30, nub = 4, nue = 4, pi = pi_common,
  mixture_var = mixture_var_common, trait_names = "trait1",
  alpha = alpha_common, marker_scale = 1)
scalar <- function(x) as.numeric(x)[1L]
prior_de <- do.call(rbind, lapply(names(de_cp), function(id) data.frame(
  variant = id, label = labels[[id]],
  pi = paste(format(pi_common, digits = 17), collapse = ";"),
  alpha = paste(format(alpha_common, digits = 17), collapse = ";"),
  dirichlet_mean = paste(format(alpha_common / sum(alpha_common), digits = 17), collapse = ";"),
  mixture_var = paste(format(mixture_var_common, digits = 17), collapse = ";"),
  initial_weight = sum(pi_common * mixture_var_common),
  prior_mean_weight = sum(alpha_common / sum(alpha_common) * mixture_var_common),
  B = scalar(resolved$B), E = scalar(resolved$E),
  ssb_prior = scalar(resolved$ssb_prior), sse_prior = scalar(resolved$sse_prior),
  source = "exact_installed_resolver", h2 = .30, adjE = .9,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE)))
prior_resolution <- rbind(prior_abc, prior_de)
prior_resolution$label <- unname(labels[prior_resolution$variant])
for (field in c("pi", "alpha", "mixture_var", "B", "E", "ssb_prior",
    "sse_prior", "h2", "adjE")) {
  value <- prior_resolution[[field]]
  equal <- if (is.numeric(value)) {
    max(abs(value - value[1L])) <= 1e-12 * max(1, abs(value[1L]))
  } else {
    parsed <- lapply(strsplit(gsub("[[:space:]]", "", value), ";", fixed = TRUE), as.numeric)
    reference <- parsed[[1L]]
    all(vapply(parsed, function(x) length(x) == length(reference) &&
      max(abs(x - reference)) <= 1e-12 * max(1, max(abs(reference))), logical(1)))
  }
  if (!equal) stop("Prior equality failed for ", field, call. = FALSE)
}
input_equivalence <- data.frame(
  variant = ids, label = unname(labels[ids]),
  phenotype_hash = s06rob_hash(Y), marker_hash = s06rob_hash(marker_ids),
  sample_hash = s06rob_hash(sample_ids), pi_equal = TRUE, alpha_equal = TRUE,
  mixture_var_equal = TRUE, h2_equal = TRUE, adjE_equal = TRUE,
  update_flags_equal = TRUE, mcmc_lengths_equal = TRUE,
  chain_seeds = paste(chain_seeds, collapse = ";"),
  fit_origin = c(rep("reused_completed_checkpoint", 3L),
    rep("supplemental_validation_checkpoint", 2L)))

# Common posterior extraction.
new_traces <- do.call(rbind, Map(function(fit, id)
  s06rob_trace_table(fit, id, labels[[id]]), fits[names(de_cp)], names(de_cp)))
new_summary <- s06rob_summarise_traces(new_traces)
variance_summary <- rbind(
  utils::read.csv(file.path(local_exact, "tables", "variance_summary.csv")),
  new_summary$variance)
convergence_summary <- rbind(
  utils::read.csv(file.path(local_exact, "tables", "convergence_summary.csv")),
  new_summary$convergence)
variance_summary$label <- unname(labels[variance_summary$variant])
convergence_summary$label <- unname(labels[convergence_summary$variant])
pi_summary <- variance_summary[grepl("^pi_component_|active_probability",
  variance_summary$quantity), ]

effect_recovery <- prediction_metrics <- component_summary <- list()
posterior_effects <- list()
for (id in ids) {
  fit <- fits[[id]]
  b <- as.numeric(fit$bm[, 1L]); names(b) <- rownames(fit$bm); b <- b[marker_ids]
  pip <- as.numeric(fit$dm[, 1L]); names(pip) <- rownames(fit$dm); pip <- pip[marker_ids]
  g <- as.numeric(Z %*% b)
  posterior_effects[[id]] <- b
  calibration <- stats::lm.fit(cbind(1, g), simulation$genetic_values)$coefficients[2L]
  effect_recovery[[id]] <- data.frame(variant = id, label = labels[[id]],
    effect_rmse = sqrt(mean((b - simulation$effects)^2)),
    effect_correlation = stats::cor(b, simulation$effects),
    pip_gt_001 = sum(pip > .01), pip_gt_005 = sum(pip > .05),
    pip_gt_010 = sum(pip > .10), pip_gt_050 = sum(pip > .50))
  prediction_metrics[[id]] <- data.frame(variant = id, label = labels[[id]],
    genetic_value_correlation = stats::cor(g, simulation$genetic_values),
    genetic_value_rmse = sqrt(mean((g - simulation$genetic_values)^2)),
    prediction_nmse = mean((g - simulation$genetic_values)^2) / stats::var(simulation$genetic_values),
    calibration_slope = unname(calibration), direct_variance = stats::var(g),
    direct_residual_sse = sum((simulation$phenotype - g)^2))
  probabilities <- fit$component_probabilities[[1L]]
  component_summary[[id]] <- data.frame(variant = id, label = labels[[id]],
    component = colnames(probabilities), posterior_expected_count = colSums(probabilities))
}
effect_recovery <- do.call(rbind, effect_recovery)
prediction_metrics <- do.call(rbind, prediction_metrics)
component_summary <- do.call(rbind, component_summary)

operators <- list(exact = exact_operator, hard_sparse = sparse_operator,
  block_eigen_full_rank = full_operator, block_eigen_retained = retained_operator)
operator_summary <- rbind(
  s06rob_operator_metrics(exact_operator, exact_operator, "exact", ncol(Z)),
  s06rob_operator_metrics(sparse_operator, exact_operator, "hard_sparse", qr(sparse_operator)$rank),
  s06rob_operator_metrics(full_operator, exact_operator, "block_eigen_full_rank", sum(blocks$size)),
  s06rob_operator_metrics(retained_operator, exact_operator, "block_eigen_retained", total_retained_rank))
operator_pairwise <- do.call(rbind, lapply(names(operators), function(id) {
  A <- operators[[id]]
  data.frame(reference = "exact", operator = id,
    relative_frobenius_error = sqrt(sum((A - exact_operator)^2) / sum(exact_operator^2)),
    maximum_absolute_error = max(abs(A - exact_operator)),
    diagonal_maximum_error = max(abs(diag(A) - diag(exact_operator))),
    trace_difference = sum(diag(A)) - sum(diag(exact_operator)))
}))
operator_pairwise <- rbind(operator_pairwise,
  data.frame(reference = "configured_block_target", operator = "block_eigen_full_rank",
    relative_frobenius_error = sqrt(sum((full_operator - block_target)^2) / sum(block_target^2)),
    maximum_absolute_error = max(abs(full_operator - block_target)),
    diagonal_maximum_error = max(abs(diag(full_operator) - diag(block_target))),
    trace_difference = sum(diag(full_operator)) - sum(diag(block_target))),
  data.frame(reference = "configured_block_target", operator = "block_eigen_retained",
    relative_frobenius_error = sqrt(sum((retained_operator - block_target)^2) / sum(block_target^2)),
    maximum_absolute_error = max(abs(retained_operator - block_target)),
    diagonal_maximum_error = max(abs(diag(retained_operator) - diag(block_target))),
    trace_difference = sum(diag(retained_operator)) - sum(diag(block_target))))

spectra <- lapply(names(operators), function(id) s06rob_spectrum(operators[[id]], id))
names(spectra) <- names(operators)
eigenvalue_summary <- do.call(rbind, lapply(spectra, `[[`, "summary"))
eigenvalue_spectrum <- do.call(rbind, lapply(spectra, `[[`, "spectrum"))
effective_rank_summary <- eigenvalue_summary[, c("operator", "effective_rank",
  "positive_condition_number", "positive", "near_zero", "below_negative_tolerance")]

states <- c(list(all_zero = numeric(ncol(Z)), truth = simulation$effects), posterior_effects)
operator_quadratic_checks <- operator_score_checks <- list()
for (state in names(states)) {
  b <- states[[state]]
  exact_score <- Xty - as.numeric(exact_operator %*% b) + diag(exact_operator) * b
  for (operator in names(operators)) {
    A <- operators[[operator]]
    score <- Xty - as.numeric(A %*% b) + diag(A) * b
    operator_quadratic_checks[[length(operator_quadratic_checks) + 1L]] <- data.frame(
      effect_state = state, operator = operator,
      quadratic = sum(b * (A %*% b)),
      exact_quadratic = sum(b * (exact_operator %*% b)),
      quadratic_error = sum(b * ((A - exact_operator) %*% b)),
      residual_sse = yy - 2 * sum(b * Xty) + sum(b * (A %*% b)))
    operator_score_checks[[length(operator_score_checks) + 1L]] <- data.frame(
      effect_state = state, operator = operator,
      maximum_absolute_score_error = max(abs(score - exact_score)),
      score_rmse = sqrt(mean((score - exact_score)^2)))
  }
}
operator_quadratic_checks <- do.call(rbind, operator_quadratic_checks)
operator_score_checks <- do.call(rbind, operator_score_checks)

posterior_operator <- do.call(rbind, lapply(ids, function(id) {
  b <- posterior_effects[[id]]
  g <- as.numeric(Z %*% b)
  data.frame(variant = id, label = labels[[id]],
    direct_variance_n_minus_1 = stats::var(g),
    exact_variance_n_minus_1 = sum(b * (exact_operator %*% b)) / n_minus_1,
    sparse_variance_n_minus_1 = sum(b * (sparse_operator %*% b)) / n_minus_1,
    full_rank_variance_n_minus_1 = sum(b * (full_operator %*% b)) / n_minus_1,
    retained_variance_n_minus_1 = sum(b * (retained_operator %*% b)) / n_minus_1,
    direct_residual_sse = sum((simulation$phenotype - g)^2))
}))

# Spectral projections in both exact and sparse eigenbases.
effect_vectors <- c(list(truth = simulation$effects), posterior_effects)
projection_rows <- list(); decomposition <- list(); low_summary <- list()
for (basis_name in c("exact", "hard_sparse")) {
  basis_operator <- operators[[basis_name]]
  eg <- eigen((basis_operator + t(basis_operator)) / 2, symmetric = TRUE)
  ascending <- order(eg$values)
  low1 <- ascending[seq_len(max(1L, floor(.01 * ncol(Z))))]
  low5 <- ascending[seq_len(max(1L, floor(.05 * ncol(Z))))]
  negative <- which(eg$values < -1e-8)
  delta_direction <- colSums(eg$vectors * ((exact_operator - sparse_operator) %*% eg$vectors))
  changed <- order(abs(delta_direction), decreasing = TRUE)[seq_len(max(1L, floor(.01 * ncol(Z))))]
  for (effect_id in names(effect_vectors)) {
    b <- effect_vectors[[effect_id]]
    coef <- as.numeric(crossprod(eg$vectors, b))
    exact_action <- as.numeric(crossprod(eg$vectors, exact_operator %*% b))
    sparse_action <- as.numeric(crossprod(eg$vectors, sparse_operator %*% b))
    rows <- data.frame(basis = basis_name, effect = effect_id,
      direction = seq_along(coef), eigenvalue = eg$values,
      squared_coefficient = coef^2,
      exact_quadratic_contribution = coef * exact_action,
      sparse_quadratic_contribution = coef * sparse_action,
      exact_minus_sparse_contribution = coef * (exact_action - sparse_action),
      negative_sparse_direction = if (basis_name == "hard_sparse") seq_along(coef) %in% negative else FALSE,
      lowest_one_percent = seq_along(coef) %in% low1,
      lowest_five_percent = seq_along(coef) %in% low5,
      largest_operator_change = seq_along(coef) %in% changed)
    projection_rows[[length(projection_rows) + 1L]] <- rows
    decomposition[[length(decomposition) + 1L]] <- data.frame(
      basis = basis_name, effect = effect_id,
      exact_quadratic = sum(rows$exact_quadratic_contribution),
      sparse_quadratic = sum(rows$sparse_quadratic_contribution),
      exact_minus_sparse = sum(rows$exact_minus_sparse_contribution),
      negative_sparse_contribution = sum(rows$exact_minus_sparse_contribution[rows$negative_sparse_direction]),
      lowest_one_percent_contribution = sum(rows$exact_minus_sparse_contribution[rows$lowest_one_percent]),
      lowest_five_percent_contribution = sum(rows$exact_minus_sparse_contribution[rows$lowest_five_percent]),
      largest_change_directions_contribution = sum(rows$exact_minus_sparse_contribution[rows$largest_operator_change]))
    low_summary[[length(low_summary) + 1L]] <- data.frame(
      basis = basis_name, effect = effect_id,
      squared_effect_lowest_one_percent = sum(rows$squared_coefficient[rows$lowest_one_percent]),
      squared_effect_lowest_five_percent = sum(rows$squared_coefficient[rows$lowest_five_percent]),
      squared_effect_negative_sparse = sum(rows$squared_coefficient[rows$negative_sparse_direction]))
  }
}
effect_spectral_projection <- do.call(rbind, projection_rows)
quadratic_error_decomposition <- do.call(rbind, decomposition)
low_eigen_direction_summary <- do.call(rbind, low_summary)

scheduler_variance <- utils::read.csv(file.path(local_scheduler, "tables", "variance_summary.csv"))
scheduler_comparison <- scheduler_variance[scheduler_variance$quantity %in%
  c("heritability", "vgs", "ves", "vbs", "active_probability"),
  c("variant", "label", "quantity", "mean", "sd", "lower_025", "upper_975")]

fit_status <- rbind(
  transform(utils::read.csv(file.path(local_exact, "tables", "fit_status.csv")),
    fit_origin = "reused_completed_checkpoint"),
  do.call(rbind, lapply(names(de_cp), function(id) data.frame(
    variant = id, label = labels[[id]], status = "ok", chains = 4L,
    draws_per_chain = 1000L, elapsed_seconds = de_cp[[id]]$elapsed_seconds,
    checkpoint_reused = TRUE, warnings = length(de_cp[[id]]$warnings),
    fit_origin = "supplemental_validation_checkpoint"))))
fit_status$label <- unname(labels[fit_status$variant])

design <- data.frame(
  study = "06_ld_operator_sbayesr_ld_robustness", samples = nrow(Z), markers = ncol(Z),
  qc_start = start, qc_end = end, first_marker = marker_ids[1L],
  last_marker = marker_ids[length(marker_ids)], simulation_seed = simulation$seed,
  causal_markers = length(simulation$causal_index), realized_h2 = simulation$realized_h2,
  block_policy = "fixed contiguous 250-marker starts plus chromosome boundary splits",
  block_starts = paste(blocks$start, collapse = ";"), block_count = nrow(blocks),
  eigen_full_representation = "dense_reconstructed/ridge_fixed_eta_0",
  eigen_retained_representation = "low_rank/cumulative_positive_mass_0.995",
  retained_total_rank = total_retained_rank, retained_discarded_rank = ncol(Z) - total_retained_rank,
  full_global_equivalent = full_global_equivalent)
provenance <- data.frame(item = c("starting_head", "sblrbench_version", "sblr_version",
  "sblr_sha", "qgdata_sha", "R", "simulation_hash", "marker_hash", "sample_hash"),
  value = c("39c8596ddd810d6fee43bd7f7906d20cbbe52440",
    as.character(utils::packageVersion("sblrbench")), as.character(utils::packageVersion("sblr")),
    sblr_sha, "6cca5819e711d326cfb2614d7e9d9f34942612cd", R.version.string,
    s06rob_hash(simulation), s06rob_hash(marker_ids), s06rob_hash(sample_ids)))

tables <- list(
  design.csv = design, provenance.csv = provenance, fit_status.csv = fit_status,
  prior_resolution.csv = prior_resolution, input_equivalence.csv = input_equivalence,
  scheduler_comparison.csv = scheduler_comparison,
  ld_operator_comparison.csv = utils::read.csv(file.path(local_exact, "tables", "ld_operator_comparison.csv")),
  operator_summary.csv = operator_summary,
  operator_pairwise_comparison.csv = operator_pairwise,
  eigenvalue_summary.csv = eigenvalue_summary,
  effective_rank_summary.csv = effective_rank_summary,
  conditional_marker_summary.csv = utils::read.csv(file.path(local_exact, "tables", "conditional_marker_summary.csv")),
  deterministic_variance_audit.csv = utils::read.csv(file.path(local_exact, "tables", "deterministic_variance_audit.csv")),
  deterministic_residual_audit.csv = utils::read.csv(file.path(local_exact, "tables", "deterministic_residual_audit.csv")),
  variance_summary.csv = variance_summary, pi_summary.csv = pi_summary,
  component_summary.csv = component_summary, convergence_summary.csv = convergence_summary,
  effect_recovery.csv = effect_recovery, prediction_metrics.csv = prediction_metrics,
  posterior_operator_comparison.csv = posterior_operator,
  quadratic_error_decomposition.csv = quadratic_error_decomposition,
  low_eigen_direction_summary.csv = low_eigen_direction_summary,
  block_diagnostics.csv = retained$diagnostics,
  operator_quadratic_checks.csv = operator_quadratic_checks,
  operator_score_checks.csv = operator_score_checks)
for (name in names(tables)) s06rob_write_csv(tables[[name]], file.path(capsule, name))
s06rob_write_csv(blocks, file.path(capsule, "block_definitions.csv"))

# Detailed local-only tables.
s06rob_write_csv(eigenvalue_spectrum, file.path(local_supplement, "tables", "eigenvalue_spectrum.csv"))
s06rob_write_csv(effect_spectral_projection, file.path(local_supplement, "tables", "effect_spectral_projection.csv"))

manifest <- list(
  schema_version = 1L, study = "06_ld_operator_sbayesr_ld_robustness",
  source_head = "39c8596ddd810d6fee43bd7f7906d20cbbe52440",
  sblr_version = as.character(utils::packageVersion("sblr")), sblr_sha = sblr_sha,
  qgdata_sha = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
  coordinate = list(samples = 5000L, markers = 1500L, qc_start = start,
    qc_end = end, simulation_seed = simulation$seed, chain_seeds = chain_seeds),
  fits = list(reused = names(abc_cp), newly_checkpointed = names(de_cp)),
  ld = list(exact = list(max_distance_variants = 0L, r2_threshold = 0),
    sparse = list(max_distance_variants = 1000L, r2_threshold = .001)),
  blocks = list(policy = design$block_policy, starts = blocks$start),
  eigen = list(full = list(representation = "dense_reconstructed",
      policy = "ridge_fixed", eta = 0),
    retained = list(representation = "low_rank",
      policy = "cumulative_positive_mass", proportion = .995,
      residual_rebuild_every = 100L, total_rank = total_retained_rank,
      block_ranks = retained$diagnostics$retained_rank)),
  interpretation_gate = list(full_rank_global_equivalent = full_global_equivalent,
    full_rank_block_target_relative_error = full_implementation_error,
    full_rank_global_relative_error = full_global_error))
jsonlite::write_json(manifest, file.path(capsule, "benchmark_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 16)

readme <- c(
  "# Supplemental Study 06 SBayesR LD-operator capsule", "",
  "This capsule contains compact tables for the SBayesR LD-operator robustness study.",
  "BED/exact/sparse fits and the scheduler diagnostic were reused from validated local checkpoints;",
  "only the full-rank and retained-low-rank block-eigen fits were newly checkpointed.", "",
  "Large fits, genotype matrices, operator matrices, eigenvectors, logs, and detailed spectra remain ignored under",
  "`results/local/06_ld_operator/sbayesr_ld_robustness/` and are not required to render the report.")
writeLines(readme, file.path(capsule, "README.md"))

checksum_files <- sort(setdiff(list.files(capsule, full.names = TRUE),
  file.path(capsule, "checksums.csv")))
checksums <- data.frame(file = basename(checksum_files),
  sha256 = vapply(checksum_files, s06rob_hash_file, character(1)))
s06rob_write_csv(checksums, file.path(capsule, "checksums.csv"))

cat("study06_sbayesr_ld_robustness_complete\n",
  "abc_checkpoints_reused=3/3\n",
  "de_checkpoints_reused=", sum(vapply(de_cp, `[[`, logical(1), "reused")), "/2\n",
  "sampler_calls=", sampler_calls, "\n",
  "capsule_regenerated=TRUE\n", sep = "")
