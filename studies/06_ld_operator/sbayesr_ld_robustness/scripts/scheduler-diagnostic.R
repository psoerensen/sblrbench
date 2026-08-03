#!/usr/bin/env Rscript

# One-coordinate comparison of scheduled BED BayesR, full-sweep BED BayesR,
# and CSR SBayesR. Generated fits and evidence remain under results/local.

options(stringsAsFactors = FALSE, warn = 1)
`%||%` <- function(x, y) if (is.null(x)) y else x

required_sblr_sha <- "02e8c74baa906e83c4a08d42a9cc6339b4e81072"
isolated_library <- file.path("results", "local", "current_benchmark_refresh", "rlib")
if (!dir.exists(isolated_library)) stop("Missing isolated library: ", isolated_library, call. = FALSE)
.libPaths(unique(c(normalizePath(isolated_library, winslash = "/"), .libPaths())))
required <- c("sblr", "qgg", "digest", "jsonlite", "posterior", "ggplot2")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)

find_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) && file.exists(file.path(path, "sblrbench.Rproj"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot locate repository root.", call. = FALSE)
    path <- parent
  }
}
root <- find_root()
setwd(root)
devtools::load_all(root,quiet=TRUE)
Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

local_root <- file.path("results", "local", "bed_vs_csr_bayesr_scheduler")
dirs <- file.path(local_root, c("inputs", "fits", "tables", "figures", "logs"))
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(local_root, "logs", "diagnostic.log")
log_line <- function(...) cat(format(Sys.time(), tz = "UTC", usetz = TRUE), " | ",
  paste0(..., collapse = ""), "\n", file = log_path, append = TRUE, sep = "")
atomic_rds <- function(x, path) {
  tmp <- tempfile(".checkpoint-", dirname(path), ".rds")
  saveRDS(x, tmp, version = 3)
  if (!file.rename(tmp, path)) { unlink(tmp); stop("Cannot write ", path, call. = FALSE) }
  invisible(path)
}
write_csv <- function(x, name) {
  path <- file.path(local_root, "tables", name)
  tmp <- tempfile(".table-", dirname(path), ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) { unlink(tmp); stop("Cannot write ", path, call. = FALSE) }
  invisible(path)
}
hash_object <- benchmark_hash_object
hash_file <- function(x) unname(benchmark_file_sha256(x))
collapse_num <- function(x) paste(format(as.numeric(x), digits = 17, scientific = FALSE, trim = TRUE), collapse = ";")

started <- Sys.time()
sampler_calls <- 0L
setup_files <- c("studies/01_finemapping/setup_example_data.R",
  "studies/five_replicate_helpers.R", "studies/03_parameter_estimation/spec.R")
for (f in setup_files[1:2]) sys.source(f, envir = environment())
config <- read_benchmark_spec(setup_files[3])
legacy_config <- list(chr=config$data$chromosome,trait=config$data$trait,
  sample_limit=config$data$sample_limit,example_data=config$data$example_data,
  qc=config$markers$qc,sparse_ld=config$data$sparse_ld,
  simulation=list(h2=config$controls$simulation$h2,
    n_causal=config$controls$simulation$n_causal,
    base_seed=config$seeds$simulation_base,architectures=config$scenarios),
  oracle_tolerance=config$validation$oracle_tolerance)

sblr_desc <- utils::packageDescription("sblr")
sblr_sha <- sblr_desc$RemoteSha %||% sblr_desc$GithubSHA1 %||% NA_character_
if (!identical(sblr_sha, required_sblr_sha))
  stop("Installed sblr SHA is ", sblr_sha, "; required ", required_sblr_sha, call. = FALSE)
capsule <- file.path("results", "reference", "03_parameter_estimation", "current")
manifest03 <- jsonlite::read_json(file.path(capsule, "benchmark_manifest.json"), simplifyVector = TRUE)
if (!identical(manifest03$sblr_commit, sblr_sha) || !identical(manifest03$qgdata_commit, config$data$example_data$commit))
  stop("Frozen Study 03 provenance differs from the diagnostic contract.", call. = FALSE)

seed_registry <- utils::read.csv(file.path(capsule, "seed_registry.csv"), check.names = FALSE)
seed_rows <- seed_registry[seed_registry$architecture == "sparse_mixture" &
  seed_registry$replicate == 1L & seed_registry$method == "st_csr_sbayesr", ]
seed_rows <- seed_rows[order(seed_rows$chain), ]
chain_seeds <- c(140104L, 240104L, 340104L, 440104L)
if (nrow(seed_rows) != 4L || !identical(as.integer(seed_rows$simulation_seed), rep(7002L, 4L)) ||
    !identical(as.integer(seed_rows$fit_seed), rep(40104L, 4L)) ||
    !identical(as.integer(seed_rows$chain_seed), chain_seeds))
  stop("Frozen seed registry does not match the prescribed coordinate.", call. = FALSE)

# Reuse and validate the prior diagnostic's exact coordinate checkpoint.
old_coordinate_path <- file.path("results", "local", "sbayesr_gctb_diagnostic", "inputs", "coordinate.rds")
if (!file.exists(old_coordinate_path)) stop("Missing prior coordinate checkpoint: ", old_coordinate_path, call. = FALSE)
coordinate <- readRDS(old_coordinate_path)
if(!identical(coordinate$checkpoint_schema,"sblrbench-semantic-v2"))
  stop(paste("Legacy source-hashed diagnostic checkpoint detected.",
    "This checkpoint schema has been retired and is not reusable under",
    "the shared semantic checkpoint framework."),call.=FALSE)
simulation <- coordinate$simulation
stats <- coordinate$summary_stats
Z <- simulation$data$genotypes
if (!identical(dim(Z), c(5000L, 37991L)) || !identical(stats$n, 5000L) ||
    !identical(length(stats$ww[[1L]]), 37991L) ||
    !identical(simulation$scenario$architecture, "sparse_mixture") ||
    !identical(simulation$scenario$replicate, 1L) ||
    !identical(simulation$provenance$seed, 7002L) ||
    !identical(rownames(Z), coordinate$sample_ids) ||
    !identical(colnames(Z), stats$marker_names) ||
    !identical(colnames(Z), simulation$data$marker_ids) ||
    sum(simulation$truth$effects[, 1L] != 0) != 50L)
  stop("Prior coordinate checkpoint does not match the exact design.", call. = FALSE)
oracle <- sblrbench::check_oracle_genetic_values(simulation, tolerance = config$validation$oracle_tolerance)
realized_h2 <- stats::var(simulation$truth$genetic_values[, 1L]) /
  (stats::var(simulation$truth$genetic_values[, 1L]) + stats::var(simulation$truth$residuals[, 1L]))
if (!oracle$ok || abs(realized_h2 - 0.30) > 1e-12) stop("Coordinate oracle or heritability validation failed.", call. = FALSE)

paths <- list(glist_path=Sys.getenv("SBLR_BENCH_GLIST",""),
  data_dir=file.path("results","local","03_parameter_estimation","data"),
  output_dir=file.path("results","local","03_parameter_estimation","genotype_setup"))
expected_data <- file.path(paths$data_dir, names(config$data$example_data$md5))
if (!all(file.exists(expected_data)) || !identical(unname(tools::md5sum(expected_data)), unname(config$data$example_data$md5)))
  stop("Pinned qgdata cache is incomplete or has changed.", call. = FALSE)
base_glist <- .study01_load_glist(paths)
marker_info <- .study01_run_qc(base_glist, legacy_config)
sample_ids <- .study01_selected_ids(base_glist, config$data$sample_limit)
working_glist <- .study01_set_rsids_ld(base_glist, config$data$chromosome, marker_info$marker_ids)
sparse_cache <- file.path(paths$output_dir, "ld_sparse_bed_glist.rds")
if (!file.exists(sparse_cache)) stop("Existing Study 03 sparse-LD cache is missing; refusing LD construction.", call. = FALSE)
Glist <- readRDS(sparse_cache)
ld_prefix <- Glist$sparseLD$prefix
ld_files <- sort(Sys.glob(paste0(ld_prefix, "*")))
if (!length(ld_files) || !identical(marker_info$marker_ids, simulation$data$marker_ids) ||
    !identical(sample_ids, simulation$data$sample_ids)) stop("Current genotype/LD alignment differs from the checkpoint.", call. = FALSE)
recomputed_stats <- sblrbench:::parameter_summary_stats(simulation,Glist,config)
if (!identical(recomputed_stats$marker_names, stats$marker_names) ||
    !isTRUE(all.equal(recomputed_stats$yy, stats$yy, tolerance = 0)) ||
    !isTRUE(all.equal(recomputed_stats$wy, stats$wy, tolerance = 0)) ||
    !isTRUE(all.equal(recomputed_stats$ww, stats$ww, tolerance = 0)))
  stop("Recomputed summary statistics differ from the prior coordinate checkpoint.", call. = FALSE)

phenotype_hash <- hash_object(simulation$truth$phenotypes)
effect_hash <- hash_object(simulation$truth$effects)
sample_hash <- hash_object(simulation$data$sample_ids)
marker_hash <- hash_object(simulation$data$marker_ids)
ld_hashes <- vapply(ld_files, hash_file, character(1))
pi_common <- c(0.99, rep(0.01 / 3, 3))
mixture_var_common <- c(0, 0.01, 0.1, 1)
alpha_common <- pi_common * 5e5
common <- list(pi = pi_common, alpha = alpha_common, mixture_var = mixture_var_common,
  h2 = 0.30, adjE = 0.9, updateB = TRUE, updateE = TRUE, updatePi = TRUE,
  nburn = 250L, nit = 2000L, nthin = 1L, nchains = 4L, ncores = 4L,
  seed = 40104L, chain_seeds = chain_seeds, nub = 4, nue = 4)
schedulers <- list(
  bed_scheduled_current = list(full_sweep_every = 10L, null_skip_base = 50L,
    null_skip_max = 200L, candidate_threshold = 1e-3, candidate_lifetime = 20L,
    skip_nulls_burnin_only = FALSE),
  bed_full_sweep = list(full_sweep_every = 0L, null_skip_base = 1L,
    null_skip_max = 1L, candidate_threshold = 0, candidate_lifetime = 0L,
    skip_nulls_burnin_only = FALSE),
  csr_current = list(full_sweep_every = NA_integer_, null_skip_base = NA_integer_,
    null_skip_max = NA_integer_, candidate_threshold = NA_real_, candidate_lifetime = NA_integer_,
    skip_nulls_burnin_only = NA))
labels <- c(bed_scheduled_current = "BED scheduled", bed_full_sweep = "BED full sweep", csr_current = "CSR current")

input_identity <- list(diagnostic_schema_version = 2L,
  prior_coordinate_hash = coordinate$semantic_hash,
  simulation_seed = 7002L, phenotype_hash = phenotype_hash, effect_hash = effect_hash,
  sample_hash = sample_hash, marker_hash = marker_hash, qgdata_sha = config$data$example_data$commit,
  sblr_sha = sblr_sha, ld_prefix = ld_prefix,
  ld_hashes = ld_hashes, common = common, schedulers = schedulers)
input_identity <- benchmark_semantic_checkpoint_identity(
  "study06-sbayesr-scheduler",input_identity)
input_hash <- benchmark_semantic_checkpoint_hash(input_identity)
input_checkpoint <- file.path(local_root, "inputs", "coordinate_validation.rds")
input_reused <- FALSE
if (file.exists(input_checkpoint)) {
  old <- benchmark_load_semantic_checkpoint(input_checkpoint,input_hash)$value
  input_reused <- TRUE
} else atomic_rds(list(checkpoint_schema="sblrbench-semantic-v2",
  semantic_hash=input_hash,input_hash = input_hash, identity = input_identity,
  validation = list(oracle = oracle, realized_h2 = realized_h2, summary_stats_exact = TRUE)), input_checkpoint)

fit_one <- function(id) {
  identity <- benchmark_semantic_checkpoint_identity(
    "study06-sbayesr-scheduler-fit",
    list(coordinate_hash = input_hash, scenario = "sparse_mixture", replicate = 1L,
      trait = config$data$trait, sample_hash = sample_hash, marker_hash = marker_hash,
      phenotype_hash = phenotype_hash, true_effect_hash = effect_hash,
      ld_prefix = ld_prefix,
      ld_file_hashes = ld_hashes, variant = id, method_controls = common,
      scheduler_controls = schedulers[[id]], simulation_seed = 7002L,
      fit_seed = 40104L, chain_seeds = chain_seeds, sblr_sha = sblr_sha,
      qgdata_sha = config$data$example_data$commit,
      backend = if (grepl("^bed", id)) "stblr_bed/bayesr" else "stblr_csr/sbayesr"))
  identity_hash <- benchmark_semantic_checkpoint_hash(identity)
  path <- file.path(local_root, "fits", paste0(id, ".rds"))
  if (file.exists(path)) {
    checkpoint <- benchmark_load_semantic_checkpoint(path, identity_hash)$value
    checkpoint$checkpoint_reused <- TRUE
    log_line("reused fit checkpoint ", id)
    return(checkpoint)
  }
  sampler_calls <<- sampler_calls + 1L
  log_line("starting fit ", id)
  args <- c(common, list(keep_chains = TRUE, convergence = "extended",
    convergence_control = list(warn = FALSE, keep_traces = TRUE, extended_groups = "probability"),
    verbose = FALSE))
  t0 <- proc.time()[["elapsed"]]
  warnings <- messages <- character()
  fit <- withCallingHandlers({
    if (grepl("^bed", id)) {
      do.call(sblr::stblr_bed, c(list(y = simulation$truth$phenotypes, Glist = Glist,
        rows = seq_len(nrow(Z)), method = "bayesr"), args, schedulers[[id]]))
    } else {
      do.call(sblr::stblr_csr, c(list(stats = stats, Glist = Glist, method = "sbayesr",
        updateLDswap = FALSE, estimate_maf_effect_s = FALSE, maf_effect_s = NULL), args))
    }
  }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
     message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") })
  elapsed <- unname(proc.time()[["elapsed"]] - t0)
  checkpoint <- list(checkpoint_schema = "sblrbench-semantic-v2", semantic_hash = identity_hash,
    identity_hash = identity_hash, identity = identity,
    fit = fit, elapsed_seconds = elapsed, warnings = warnings, messages = messages,
    checkpoint_reused = FALSE, completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  atomic_rds(checkpoint, path)
  log_line("completed fit ", id, " in ", elapsed, " seconds")
  checkpoint
}

fit_ids <- names(schedulers)
checkpoints <- setNames(lapply(fit_ids, fit_one), fit_ids)
fits <- lapply(checkpoints, `[[`, "fit")
if (any(vapply(fits, is.null, logical(1)))) stop("A primary fit is missing.", call. = FALSE)
if (!identical(as.integer(fits$bed_full_sweep$input$full_sweep_every), 0L))
  stop("BED full-sweep backend did not record full_sweep_every = 0.", call. = FALSE)

# Verify exact common prior resolution and scientifically relevant controls.
get_scalar <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)[1L]
prior_resolver <- getFromNamespace(".make_stblr_bayesr_priors", "sblr")
resolved_common <- prior_resolver(vy = as.numeric(stats$yy) / (as.numeric(stats$n) - 1),
  m = length(stats$ww[[1L]]), h2 = common$h2, nub = common$nub, nue = common$nue,
  pi = common$pi, mixture_var = common$mixture_var, trait_names = "trait1",
  alpha = common$alpha, marker_scale = 1)
prior_resolution <- do.call(rbind, lapply(fit_ids, function(id) {
  z <- fits[[id]]$input
  metadata_complete <- all(vapply(c("B", "E", "ssb_prior", "sse_prior"),
    function(nm) !is.null(z[[nm]]), logical(1)))
  data.frame(variant = id, label = labels[[id]], pi = collapse_num(z$pi), alpha = collapse_num(z$alpha),
    dirichlet_prior_mean = collapse_num(z$alpha / sum(z$alpha)), mixture_var = collapse_num(z$mixture_var),
    initial_mixture_weight = sum(z$pi * z$mixture_var),
    prior_mean_mixture_weight = sum((z$alpha / sum(z$alpha)) * z$mixture_var),
    B = if (metadata_complete) get_scalar(z$B) else get_scalar(resolved_common$B),
    E = if (metadata_complete) get_scalar(z$E) else get_scalar(resolved_common$E),
    ssb_prior = if (metadata_complete) get_scalar(z$ssb_prior) else get_scalar(resolved_common$ssb_prior),
    sse_prior = if (metadata_complete) get_scalar(z$sse_prior) else get_scalar(resolved_common$sse_prior),
    resolved_value_source = if (metadata_complete) "fit_metadata" else "exact_installed_resolver_from_recorded_fit_inputs",
    fit_metadata_has_fixed_prior_quantities = metadata_complete,
    h2 = z$h2, adjE = z$adjE, updateB = z$updateB, updateE = z$updateE, updatePi = z$updatePi,
    stringsAsFactors = FALSE)
}))
equal_fields <- c("pi", "alpha", "mixture_var", "B", "E", "ssb_prior", "sse_prior", "h2", "adjE")
prior_equal <- all(vapply(equal_fields, function(nm) {
  x <- prior_resolution[[nm]]
  if (is.numeric(x)) max(abs(x - x[1L])) <= 1e-12 * max(1, abs(x[1L])) else length(unique(x)) == 1L
}, logical(1)))
strict_prior_metadata_confirmation <- all(prior_resolution$fit_metadata_has_fixed_prior_quantities)
if (!prior_equal) stop("Automatically resolved common prior quantities differ; interpretation is blocked.", call. = FALSE)
write_csv(prior_resolution, "prior_resolution.csv")

scheduler_design <- do.call(rbind, lapply(fit_ids, function(id) {
  s <- schedulers[[id]]; bed <- grepl("^bed", id); total <- common$nburn + common$nit
  data.frame(variant = id, label = labels[[id]], scheduled = identical(id, "bed_scheduled_current"),
    full_sweep_every = s$full_sweep_every, null_skip_base = s$null_skip_base,
    null_skip_max = s$null_skip_max, candidate_threshold = s$candidate_threshold,
    candidate_lifetime = s$candidate_lifetime, skip_nulls_burnin_only = s$skip_nulls_burnin_only,
    nominal_total_iterations = total,
    nominal_mandatory_full_sweeps = if (!bed || identical(s$full_sweep_every, 0L)) total else floor(total / s$full_sweep_every),
    all_markers_updated_each_iteration = !identical(id, "bed_scheduled_current"),
    pi_updated_each_iteration = TRUE, actual_marker_update_count_available = FALSE, stringsAsFactors = FALSE)
}))
write_csv(scheduler_design, "scheduler_design.csv")

# Extract retained chain traces from the public convergence bundle.
extract_traces <- function(fit, id) {
  b <- fit$convergence_traces
  if (is.null(b$values) || !identical(dim(b$values)[1:2], c(common$nit, common$nchains)))
    stop("Missing or mis-sized convergence trace bundle for ", id, call. = FALSE)
  q <- b$quantities
  rows <- vector("list", dim(b$values)[3L])
  for (j in seq_len(dim(b$values)[3L])) {
    group <- as.character(q$group[j])
    component <- if ("component_name" %in% names(q)) as.character(q$component_name[j]) else NA_character_
    quantity <- if (identical(group, "component_pi")) paste0("pi_", component) else group
    rows[[j]] <- data.frame(variant = id, label = labels[[id]], iteration = rep(seq_len(common$nit), common$nchains),
      chain = rep(seq_len(common$nchains), each = common$nit), quantity = quantity,
      component = component, value = as.vector(b$values[, , j]), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  wide <- function(nm) matrix(out$value[out$quantity == nm], nrow = common$nit, ncol = common$nchains)
  vg <- wide("vgs"); ve <- wide("ves")
  out <- rbind(out, data.frame(variant = id, label = labels[[id]],
    iteration = rep(seq_len(common$nit), common$nchains), chain = rep(seq_len(common$nchains), each = common$nit),
    quantity = "heritability", component = NA_character_, value = as.vector(vg / (vg + ve))))
  null_name <- unique(out$quantity[grepl("^pi_component_0$", out$quantity)])
  if (length(null_name) == 1L) {
    p0 <- wide(null_name)
    out <- rbind(out, data.frame(variant = id, label = labels[[id]],
      iteration = rep(seq_len(common$nit), common$nchains), chain = rep(seq_len(common$nchains), each = common$nit),
      quantity = "active_probability", component = NA_character_, value = as.vector(1 - p0)))
  }
  out
}
traces <- do.call(rbind, Map(extract_traces, fits, fit_ids))
if (any(!is.finite(traces$value))) stop("Non-finite retained scalar trace.", call. = FALSE)

diagnose <- function(x) {
  mat <- matrix(x$value, nrow = common$nit, ncol = common$nchains)
  sd_all <- stats::sd(as.numeric(mat)); mcse <- posterior::mcse_mean(mat)
  data.frame(variant = x$variant[1], label = x$label[1], quantity = x$quantity[1],
    rhat = posterior::rhat(mat), bulk_ess = posterior::ess_bulk(mat), tail_ess = posterior::ess_tail(mat),
    mcse_mean = mcse, relative_mcse = if (is.finite(sd_all) && sd_all > 0) mcse / sd_all else NA_real_)
}
convergence_summary <- do.call(rbind, lapply(split(traces,
  interaction(traces$variant, traces$quantity, drop = TRUE)), diagnose))
write_csv(convergence_summary, "convergence_summary.csv")

scalar_summary <- do.call(rbind, lapply(split(traces,
  interaction(traces$variant, traces$quantity, drop = TRUE)), function(x) data.frame(
    variant = x$variant[1], label = x$label[1], quantity = x$quantity[1], draws = nrow(x),
    mean = mean(x$value), sd = stats::sd(x$value), median = stats::median(x$value),
    lower_025 = unname(stats::quantile(x$value, .025)), upper_975 = unname(stats::quantile(x$value, .975)),
    minimum = min(x$value), maximum = max(x$value))))
write_csv(scalar_summary, "variance_summary.csv")

chain_summary <- do.call(rbind, lapply(split(traces,
  interaction(traces$variant, traces$quantity, traces$chain, drop = TRUE)), function(x) data.frame(
    variant = x$variant[1], label = x$label[1], chain = x$chain[1], chain_seed = chain_seeds[x$chain[1]],
    quantity = x$quantity[1], draws = nrow(x), mean = mean(x$value), sd = stats::sd(x$value),
    final = tail(x$value, 1L), minimum = min(x$value), maximum = max(x$value))))
write_csv(chain_summary, "chain_summary.csv")

pi_summary <- scalar_summary[grepl("^pi_component_|^active_probability$", scalar_summary$quantity), ]
pi_summary$summary_type <- "retained_draw_posterior"
write_csv(pi_summary, "pi_summary.csv")

component_summary <- do.call(rbind, lapply(fit_ids, function(id) {
  fit <- fits[[id]]
  cp_list <- fit$component_probabilities %||% fit$comp_prob
  if (is.null(cp_list) || !is.list(cp_list) || is.null(dim(cp_list[[1L]])))
    stop("Aggregate component probabilities are unavailable for ", id, call. = FALSE)
  cp <- cp_list[[1L]]
  aggregate <- data.frame(variant = id, label = labels[[id]], chain = NA_integer_,
    component = colnames(cp), summary_type = "aggregate_marker_posterior_expected_count",
    count = colSums(cp), stringsAsFactors = FALSE)
  chain_rows <- do.call(rbind, lapply(seq_along(fit$chains), function(k) {
    ch <- fit$chains[[k]]; out <- list()
    if (!is.null(ch$component_probabilities)) out[[length(out) + 1L]] <- data.frame(
      variant = id, label = labels[[id]], chain = k, component = colnames(ch$component_probabilities),
      summary_type = "chain_marker_posterior_expected_count", count = colSums(ch$component_probabilities))
    state <- ch$component %||% ch$state
    if (!is.null(state)) out[[length(out) + 1L]] <- data.frame(variant = id, label = labels[[id]], chain = k,
      component = paste0("component_", 0:3), summary_type = "chain_final_sampled_membership_count",
      count = tabulate(as.integer(state) + 1L, nbins = 4L))
    do.call(rbind, out)
  }))
  rbind(aggregate, chain_rows)
}))
write_csv(component_summary, "component_summary.csv")

# Posterior-mean marker and prediction recovery using current definitions.
effect_recovery <- prediction_metrics <- list()
predictions <- list()
for (id in fit_ids) {
  b <- as.numeric(fits[[id]]$bm[, 1L]); names(b) <- rownames(fits[[id]]$bm)
  b <- b[simulation$data$marker_ids]
  truth_b <- simulation$truth$effects[, 1L]
  pip <- as.numeric(fits[[id]]$dm[, 1L]); names(pip) <- rownames(fits[[id]]$dm); pip <- pip[simulation$data$marker_ids]
  if (any(!is.finite(b)) || any(!is.finite(pip)) || any(pip < 0 | pip > 1)) stop("Invalid marker summaries for ", id, call. = FALSE)
  ghat <- as.numeric(Z %*% b); predictions[[id]] <- ghat
  gtrue <- simulation$truth$genetic_values[, 1L]; y <- simulation$truth$phenotypes[, 1L]
  calibration <- stats::lm.fit(cbind(1, ghat), gtrue)$coefficients
  effect_recovery[[id]] <- data.frame(variant = id, label = labels[[id]],
    effect_rmse = sqrt(mean((b - truth_b)^2)), effect_correlation = stats::cor(b, truth_b),
    nonzero_posterior_mean_effects = sum(b != 0), pip_gt_001 = sum(pip > .01),
    pip_gt_005 = sum(pip > .05), pip_gt_010 = sum(pip > .10), pip_gt_050 = sum(pip > .50))
  prediction_metrics[[id]] <- data.frame(variant = id, label = labels[[id]],
    genetic_value_correlation = stats::cor(ghat, gtrue), genetic_value_rmse = sqrt(mean((ghat - gtrue)^2)),
    direct_predicted_genetic_variance = stats::var(ghat), phenotype_prediction_correlation = stats::cor(ghat, y),
    prediction_nmse = mean((ghat - gtrue)^2) / stats::var(gtrue), calibration_intercept = calibration[1L],
    calibration_slope = calibration[2L])
}
effect_recovery <- do.call(rbind, effect_recovery); prediction_metrics <- do.call(rbind, prediction_metrics)
write_csv(effect_recovery, "effect_recovery.csv"); write_csv(prediction_metrics, "prediction_metrics.csv")

# Direct identities. Public chain objects expose chain posterior means and final
# component states, but not chain-final effect vectors; those are not invented.
csr <- sblr::sparseLD_read_CSR(ld_prefix, one_based = TRUE)
csr_matvec <- getFromNamespace(".stblr_csr_matvec", "sblr")
direct_var <- direct_resid <- list()
for (id in fit_ids) {
  fit <- fits[[id]]
  vectors <- list(aggregate_final_effect = as.numeric(fit$b[, 1L]), posterior_mean_effect = as.numeric(fit$bm[, 1L]))
  for (k in seq_along(fit$chains)) if (!is.null(fit$chains[[k]]$bm))
    vectors[[paste0("chain_", k, "_posterior_mean_effect")]] <- as.numeric(fit$chains[[k]]$bm)
  for (nm in names(vectors)) {
    b <- vectors[[nm]]; g <- as.numeric(Z %*% b); r <- simulation$truth$phenotypes[, 1L] - g
    ch <- if (grepl("^chain_", nm)) as.integer(sub("^chain_([0-9]+).*", "\\1", nm)) else NA_integer_
    stored_vg <- if (is.na(ch)) mean(tail(subset(traces, variant == id & quantity == "vgs")$value, common$nchains)) else
      tail(subset(traces, variant == id & quantity == "vgs" & chain == ch)$value, 1L)
    stored_ve <- if (is.na(ch)) mean(tail(subset(traces, variant == id & quantity == "ves")$value, common$nchains)) else
      tail(subset(traces, variant == id & quantity == "ves" & chain == ch)$value, 1L)
    op <- sum(b * csr_matvec(csr, b))
    direct_var[[length(direct_var) + 1L]] <- data.frame(variant = id, label = labels[[id]], effect_state = nm,
      chain = ch, denominator_n = nrow(Z), denominator_n_minus_1 = nrow(Z) - 1L,
      genotype_scale = "qgg standardized", var_Xb_n_minus_1 = stats::var(g), sum_Xb_sq_over_n = sum(g^2) / nrow(Z),
      b_XtX_b_over_n = sum(g^2) / nrow(Z), csr_b_R_b = op, stored_vgs_comparator = stored_vg,
      stored_comparison_is_draw_identity = FALSE)
    direct_resid[[length(direct_resid) + 1L]] <- data.frame(variant = id, label = labels[[id]], effect_state = nm,
      chain = ch, residual_sse = sum(r^2), residual_variance_n_minus_1 = stats::var(r),
      residual_mean_square_n = sum(r^2) / nrow(Z), stored_ves_comparator = stored_ve,
      stored_comparison_is_draw_identity = FALSE)
  }
}
direct_var <- do.call(rbind, direct_var); direct_resid <- do.call(rbind, direct_resid)
write_csv(direct_var, "direct_variance_checks.csv"); write_csv(direct_resid, "direct_residual_checks.csv")

# Primary contrasts.
contrast_pairs <- list(c("bed_full_sweep", "bed_scheduled_current"),
  c("csr_current", "bed_full_sweep"), c("csr_current", "bed_scheduled_current"))
metric_values <- merge(
  reshape(scalar_summary[scalar_summary$quantity %in% c("heritability", "vgs", "ves", "vbs", "active_probability") ,
    c("variant", "quantity", "mean")], idvar = "variant", timevar = "quantity", direction = "wide"),
  merge(effect_recovery, prediction_metrics, by = c("variant", "label")), by = "variant")
names(metric_values) <- sub("^mean\\.", "", names(metric_values))
component_means <- reshape(pi_summary[grepl("^pi_component_", pi_summary$quantity), c("variant", "quantity", "mean")],
  idvar = "variant", timevar = "quantity", direction = "wide")
names(component_means) <- sub("^mean\\.", "", names(component_means))
metric_values <- merge(metric_values, component_means, by = "variant")
active_counts <- aggregate(count ~ variant, component_summary[component_summary$summary_type == "aggregate_marker_posterior_expected_count" & component_summary$component != "component_0", ], sum)
names(active_counts)[2L] <- "active_marker_count"
metric_values <- merge(metric_values, active_counts, by = "variant")
contrast_metrics <- setdiff(names(metric_values), c("variant", "label"))
variant_contrasts <- do.call(rbind, lapply(contrast_pairs, function(p) do.call(rbind, lapply(contrast_metrics, function(nm) {
  data.frame(contrast = paste(p[1], "-", p[2]), focal = p[1], reference = p[2], metric = nm,
    focal_value = metric_values[metric_values$variant == p[1], nm],
    reference_value = metric_values[metric_values$variant == p[2], nm],
    difference = metric_values[metric_values$variant == p[1], nm] - metric_values[metric_values$variant == p[2], nm])
}))))
write_csv(variant_contrasts, "variant_contrasts.csv")

# Frozen benchmark anchors (not posterior samples from this diagnostic).
pe <- utils::read.csv(file.path(capsule, "parameter_estimates.csv"), check.names = FALSE)
anchor_methods <- c("st_bed_bayesr", "st_csr_sbayesr", "st_bed_bayesc", "st_csr_sbayesc")
frozen <- pe[pe$architecture == "sparse_mixture" & pe$replicate == 1L & pe$method %in% anchor_methods,
  c("architecture", "replicate", "method", "estimand_id", "truth", "posterior_mean", "posterior_sd", "lower_95", "upper_95", "status")]
frozen$evidence_source <- "frozen Study 03 current capsule"
write_csv(frozen, "frozen_benchmark_anchors.csv")

fit_status <- do.call(rbind, lapply(fit_ids, function(id) data.frame(variant = id, label = labels[[id]],
  status = "ok", completed_chains = length(fits[[id]]$chains), retained_draws_per_chain = common$nit,
  elapsed_seconds = checkpoints[[id]]$elapsed_seconds, checkpoint_reused = checkpoints[[id]]$checkpoint_reused,
  warning_count = length(checkpoints[[id]]$warnings), message_count = length(checkpoints[[id]]$messages),
  sampler_called_this_invocation = !checkpoints[[id]]$checkpoint_reused, stringsAsFactors = FALSE)))
write_csv(fit_status, "fit_status.csv")
design <- data.frame(study = "03_parameter_estimation", architecture = "sparse_mixture", replicate = 1L,
  trait = "trait1", samples = 5000L, markers = 37991L, causal_markers = 50L, target_h2 = .30,
  realized_h2 = realized_h2, simulation_seed = 7002L, method_seed = 40104L,
  chain_seeds = paste(chain_seeds, collapse = ";"), nburn = 250L, nit = 2000L, nthin = 1L,
  nchains = 4L, ncores = 4L, input_hash = input_hash, input_checkpoint_reused = input_reused,
  sampler_calls_this_invocation = sampler_calls, stringsAsFactors = FALSE)
write_csv(design, "design.csv")
provenance <- data.frame(item = c("sblrbench_starting_head", "sblr_version", "sblr_sha", "qgdata_sha", "R_version",
  "platform", "prior_coordinate_hash", "phenotype_hash", "effect_hash", "sample_order_hash", "marker_order_hash", "ld_prefix"),
  value = c("39c8596ddd810d6fee43bd7f7906d20cbbe52440", as.character(utils::packageVersion("sblr")), sblr_sha,
    config$data$example_data$commit, R.version.string, R.version$platform, coordinate$semantic_hash, phenotype_hash, effect_hash,
    sample_hash, marker_hash, normalizePath(ld_prefix, winslash = "/", mustWork = FALSE)))
write_csv(provenance, "provenance.csv")

# Figures.
theme_diag <- function() ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "bottom")
save_plot <- function(p, name, width = 10, height = 6) ggplot2::ggsave(file.path(local_root, "figures", name), p,
  width = width, height = height, dpi = 140)
trace_plot <- function(quantity, file, ylab) {
  z <- traces[traces$quantity == quantity, ]
  p <- ggplot2::ggplot(z, ggplot2::aes(iteration, value, colour = factor(chain))) + ggplot2::geom_line(alpha = .75) +
    ggplot2::facet_wrap(~label, scales = "free_y") + ggplot2::labs(x = "Retained draw", y = ylab, colour = "Chain") + theme_diag()
  save_plot(p, file)
}
trace_plot("heritability", "heritability_traces.png", "Heritability")
trace_plot("vgs", "vgs_traces.png", "vgs")
trace_plot("ves", "ves_traces.png", "ves")
trace_plot("vbs", "vbs_traces.png", "vbs")
trace_plot("active_probability", "active_probability_traces.png", "Total non-null pi")
zpi <- traces[grepl("^pi_component_", traces$quantity), ]
save_plot(ggplot2::ggplot(zpi, ggplot2::aes(iteration, value, colour = factor(chain))) + ggplot2::geom_line(alpha = .7) +
  ggplot2::facet_grid(quantity ~ label, scales = "free_y") + ggplot2::labs(x = "Retained draw", y = "Component pi", colour = "Chain") + theme_diag(),
  "component_pi_traces.png", 11, 9)
for (spec in list(c("heritability", "posterior_heritability.png", "Heritability"),
                  c("vgs|ves", "posterior_vgs_ves.png", "Variance"))) {
  z <- traces[grepl(paste0("^(", spec[1], ")$"), traces$quantity), ]
  save_plot(ggplot2::ggplot(z, ggplot2::aes(label, value, colour = label)) + ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(data = aggregate(value ~ label + quantity + chain, z, mean), position = ggplot2::position_jitter(width = .08), size = 1.5) +
    ggplot2::facet_wrap(~quantity, scales = "free_y") + ggplot2::labs(x = NULL, y = spec[3], colour = NULL) + theme_diag(), spec[2])
}
save_plot(ggplot2::ggplot(component_summary[component_summary$summary_type == "chain_marker_posterior_expected_count", ],
  ggplot2::aes(component, count, colour = label, group = interaction(label, chain))) + ggplot2::geom_point() + ggplot2::geom_line() +
  ggplot2::facet_wrap(~label) + ggplot2::labs(x = NULL, y = "Posterior expected marker count", colour = NULL) + theme_diag(), "component_counts.png")
erlong <- reshape(effect_recovery[, c("label", "effect_rmse", "effect_correlation")], varying = c("effect_rmse", "effect_correlation"),
  v.names = "value", timevar = "metric", times = c("Effect RMSE", "Effect correlation"), direction = "long")
save_plot(ggplot2::ggplot(erlong, ggplot2::aes(label, value, fill = label)) + ggplot2::geom_col() + ggplot2::facet_wrap(~metric, scales = "free_y") +
  ggplot2::labs(x = NULL, y = NULL, fill = NULL) + theme_diag(), "effect_recovery.png")
pred_df <- do.call(rbind, lapply(fit_ids, function(id) data.frame(label = labels[[id]], truth = simulation$truth$genetic_values[, 1L], predicted = predictions[[id]])))
save_plot(ggplot2::ggplot(pred_df, ggplot2::aes(truth, predicted)) + ggplot2::geom_point(alpha = .15, size = .5) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::facet_wrap(~label) + theme_diag(), "true_vs_predicted_genetic_values.png")
save_plot(ggplot2::ggplot(pred_df, ggplot2::aes(predicted, truth)) + ggplot2::geom_point(alpha = .15, size = .5) +
  ggplot2::geom_smooth(method = "lm", se = FALSE) + ggplot2::facet_wrap(~label) + theme_diag(), "prediction_calibration.png")
save_plot(ggplot2::ggplot(direct_var[direct_var$effect_state == "posterior_mean_effect", ],
  ggplot2::aes(stored_vgs_comparator, var_Xb_n_minus_1, colour = label)) + ggplot2::geom_point(size = 3) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::labs(x = "Stored final-draw vgs comparator", y = "var(X posterior-mean b)", colour = NULL) + theme_diag(), "direct_vs_stored_genetic_variance.png")
save_plot(ggplot2::ggplot(direct_resid[direct_resid$effect_state == "posterior_mean_effect", ],
  ggplot2::aes(stored_ves_comparator, residual_variance_n_minus_1, colour = label)) + ggplot2::geom_point(size = 3) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::labs(x = "Stored final-draw ves comparator", y = "var(y - X posterior-mean b)", colour = NULL) + theme_diag(), "direct_vs_stored_residual.png")

writeLines(capture.output(utils::sessionInfo()), file.path(local_root, "session_info.txt"))
manifest <- list(schema_version = 1L, status = "complete", input_hash = input_hash,
  input_checkpoint_reused = input_reused, sampler_calls_this_invocation = sampler_calls,
  fit_checkpoint_reused = vapply(checkpoints, `[[`, logical(1), "checkpoint_reused"),
  prior_resolution_equal = prior_equal, strict_prior_metadata_confirmation = strict_prior_metadata_confirmation,
  runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  variants = fit_status, provenance = as.list(stats::setNames(provenance$value, provenance$item)))
jsonlite::write_json(manifest, file.path(local_root, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

# Stable tracked developer report with key numeric tables.
md_table <- function(x, digits = 5) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) format(round(v, digits), trim = TRUE, scientific = FALSE) else as.character(v))
  x[is.na(x)] <- "NA"
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(x)), collapse = "|"), "|"),
    apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
central <- scalar_summary[scalar_summary$quantity %in% c("vbs", "vgs", "ves", "heritability", "active_probability"),
  c("label", "quantity", "mean", "sd", "lower_025", "upper_975")]
conv_key <- convergence_summary[convergence_summary$quantity %in% c("vbs", "vgs", "ves", "heritability", "active_probability"), ]
anchor_key <- frozen[frozen$estimand_id %in% c("genetic_variance", "residual_variance", "heritability"),
  c("method", "estimand_id", "posterior_mean", "posterior_sd")]
component_key <- component_summary[component_summary$summary_type == "aggregate_marker_posterior_expected_count",
  c("label", "component", "count")]
worst_rhat <- max(convergence_summary$rhat, na.rm = TRUE)
minimum_bulk_ess <- min(convergence_summary$bulk_ess, na.rm = TRUE)
maximum_relative_mcse <- max(convergence_summary$relative_mcse, na.rm = TRUE)
report <- c("# BED versus CSR BayesR scheduler diagnostic", "", "## Scientific question", "",
  "This one-replicate developer diagnostic asks whether adaptive packed-BED marker scheduling explains the apparent BED BayesR advantage over CSR SBayesR. It is diagnostic rather than definitive.", "",
  "## Code-review motivation", "", "Scheduled BED performs mandatory full sweeps periodically and adaptively revisits active, candidate, and due markers between them. Full-sweep BED and CSR update every marker every iteration.", "",
  "## Exact coordinate and provenance", "", paste0("Study 03 `sparse_mixture`, replicate 1, trait 1; n=5,000; m=37,991; 50 causal markers; h2=0.30; simulation seed 7,002. Installed `sblr` ", sblr_desc$Version, " at `", sblr_sha, "`; qgdata `", config$data$example_data$commit, "`."), "",
  "## Variants and scheduler", "", md_table(scheduler_design[, c("label", "full_sweep_every", "null_skip_base", "null_skip_max", "candidate_threshold", "candidate_lifetime", "nominal_mandatory_full_sweeps")], 4), "",
  "Actual marker-update counts are unavailable from the public fit object. Mandatory full-sweep counts do not include active, candidate, or due-marker updates between scheduled sweeps.", "",
  "## Prior equality", "", paste0("Resolved numerical prior equality: **", if (prior_equal) "passed" else "failed", "**. Strict all-route fit-metadata confirmation: **", if (strict_prior_metadata_confirmation) "passed" else "unavailable", "**. The CSR public fit metadata omits `B`, `E`, `ssb_prior`, and `sse_prior`; its row is recomputed with the exact installed resolver from the recorded common fit inputs."), "", md_table(prior_resolution[, c("label", "B", "E", "ssb_prior", "sse_prior", "initial_mixture_weight", "resolved_value_source")], 8), "",
  "## Convergence", "", paste0("All registered central and component-probability quantities are supported: maximum R-hat ", format(round(worst_rhat, 5), nsmall = 5), ", minimum bulk ESS ", format(round(minimum_bulk_ess, 1), nsmall = 1), ", and maximum relative MCSE ", format(round(maximum_relative_mcse, 4), nsmall = 4), "."), "", md_table(conv_key[, c("label", "quantity", "rhat", "bulk_ess", "tail_ess", "mcse_mean", "relative_mcse")], 5), "",
  "## Variance and mixture results", "", md_table(central, 5), "", md_table(pi_summary[grepl("^pi_component_", pi_summary$quantity), c("label", "quantity", "mean", "sd")], 6), "", md_table(component_key, 3), "",
  "## Effect and prediction recovery", "", md_table(effect_recovery, 6), "", md_table(prediction_metrics, 6), "",
  "## Direct variance and residual checks", "", "The public chain records retain chain posterior-mean effects and final component states, but not chain-final effect vectors. Direct genotype and CSR-operator identities are therefore reported for aggregate-final and posterior-mean vectors; stored final-draw scalar comparisons are labelled as non-identities and must not be interpreted as posterior-mean equality.", "", md_table(direct_var[direct_var$effect_state == "posterior_mean_effect", c("label", "var_Xb_n_minus_1", "sum_Xb_sq_over_n", "csr_b_R_b", "stored_vgs_comparator")], 6), "", md_table(direct_resid[direct_resid$effect_state == "posterior_mean_effect", c("label", "residual_sse", "residual_variance_n_minus_1", "stored_ves_comparator")], 6), "",
  "## Frozen Study 03 comparison", "", "These are external frozen-capsule summaries, not samples pooled with the new diagnostic fits.", "", md_table(anchor_key, 6), "",
  "## Strongest supported interpretation", "", "**Outcome B: the scheduler does not explain the primary BED-versus-CSR BayesR discrepancy.** Full-sweep BED remains close to scheduled BED: heritability differs by -0.00189, vgs by -0.00241, ves by +0.00340, active probability by +0.0000046, effect RMSE by +0.0000056, and genetic-value correlation by -0.00042. In contrast, CSR versus full-sweep BED differs by +0.0901 in heritability, +0.1295 in vgs, -0.1301 in ves, +0.00026 in effect RMSE, and -0.0283 in genetic-value correlation. CSR has more diffuse marker inclusion (355.1 posterior expected active markers versus 303.7 for full-sweep BED) despite a similar strongly anchored total active pi.", "", "CSR phenotype-prediction correlation is higher, but its truth-based genetic-value NMSE is worse and its calibration slope is 0.849 versus approximately 1.05 for both BED routes. Thus the phenotype correlation alone is not evidence of better genetic-value recovery.", "",
  "## Limitations", "", "- One simulation replicate cannot establish a universal method ranking.", "- Scheduled and full-sweep algorithms consume RNG differently; identical chain seeds do not imply identical random streams.", "- Sparse CSR still differs from individual-level BED after scheduling is matched.", "- Native actual marker-update counters are not exposed.", "- Variance of a posterior-mean effect vector is not the posterior mean of draw-level variance.", "",
  "## Recommended next step", "", "Run a deterministic BED-versus-CSR conditional-update audit with matched small inputs and effectively complete LD. The audit should compare one-marker conditional probabilities, conditional effect moments, residual updates, and operator quadratic forms. Do not change package defaults from this single replicate alone.", "",
  "## Reproduction", "", "```powershell", "Rscript studies/06_ld_operator/sbayesr_ld_robustness/scripts/scheduler-diagnostic.R", "```", "", "Generated fits, tables, figures, logs, and manifests are local and untracked under `results/local/bed_vs_csr_bayesr_scheduler/`.")
writeLines(report, file.path(out, "legacy_scheduler_report.md"))

cat("BayesR scheduler diagnostic complete\n")
cat("input_checkpoint_reused=", input_reused, "\n", sep = "")
cat("fit_checkpoints_reused=", sum(vapply(checkpoints, `[[`, logical(1), "checkpoint_reused")), "/3\n", sep = "")
cat("sampler_calls=", sampler_calls, "\n", sep = "")
