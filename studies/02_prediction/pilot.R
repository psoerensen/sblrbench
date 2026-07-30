.study02_paths <- function() {
  glist_path <- Sys.getenv("SBLR_BENCH_GLIST", "")
  data_dir <- Sys.getenv("SBLR_BENCH_DATA_DIR", file.path("results", "local", "02_prediction", "data"))
  list(glist_path = glist_path, data_dir = data_dir,
    output_dir = file.path("results", "local", "02_prediction", "genotype_setup"))
}

.study02_set_training_af <- function(Glist, chr, marker_ids, af) {
  idx <- match(marker_ids, Glist$rsids[[chr]])
  if (anyNA(idx) || length(af) != length(idx) || any(!is.finite(af))) stop("Training allele frequencies are not aligned.", call. = FALSE)
  out <- Glist
  out$af[[chr]][idx] <- unname(af)
  out$maf[[chr]][idx] <- pmin(unname(af), 1 - unname(af))
  out$rsidsLD[[chr]] <- marker_ids
  out
}

.study02_extract_raw <- function(Glist, chr, sample_ids, marker_ids) {
  x <- qgg::getG(Glist = Glist, chr = chr, ids = sample_ids, rsids = marker_ids,
    impute = FALSE, scale = FALSE)
  if (!identical(rownames(x), sample_ids) || !identical(colnames(x), marker_ids)) stop("Raw genotype extraction lost canonical order.", call. = FALSE)
  x
}

.study02_make_ld <- function(Glist, split, marker_ids, config, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, paste0("training_ld_", split$split_id))
  cache <- paste0(prefix, "_glist.rds")
  if (file.exists(cache)) {
    x <- readRDS(cache)
    ids <- x$rsids[[config$chr]][x$sparseLD$cls[[1L]]]
    if (identical(ids, marker_ids) && identical(x$sparseLD$rows, split$train_rows) &&
        length(Sys.glob(paste0(x$sparseLD$prefix, "*")))) return(x)
  }
  x <- do.call(sblr::make_sparse_ld, c(list(Glist = Glist, rows = split$train_rows,
    out_prefix = prefix, chr = config$chr), config$sparse_ld))
  ids <- x$rsids[[config$chr]][x$sparseLD$cls[[1L]]]
  if (!identical(ids, marker_ids) || !identical(x$sparseLD$rows, split$train_rows) ||
      !identical(as.integer(x$sparseLD$reference_n), length(split$train_rows))) stop("Training sparse-LD provenance is invalid.", call. = FALSE)
  saveRDS(x, cache)
  x
}

.study02_replicate_specs <- function(config, replicate_override, architecture_override) {
  n <- if (nzchar(replicate_override)) suppressWarnings(as.integer(replicate_override)) else config$replicate_counts[["development"]]
  if (length(n) != 1L || is.na(n) || !n %in% unname(config$replicate_counts)) stop("SBLR_BENCH_REPLICATES must be 1, 5, or 10.", call. = FALSE)
  architectures <- names(config$simulation$architectures)
  if (nzchar(architecture_override)) {
    if (!architecture_override %in% architectures) stop("Unknown SBLR_BENCH_ARCHITECTURE.", call. = FALSE)
    architectures <- architecture_override
  }
  out <- list()
  for (a in architectures) for (i in seq_len(n)) out[[length(out) + 1L]] <- list(
    architecture = a, replicate = i,
    simulation_seed = as.integer(config$simulation$base_seed +
      match(a, names(config$simulation$architectures)) * 1000L + i))
  out
}

.study02_simulate <- function(spec, Z_all, config) {
  architecture <- config$simulation$architectures[[spec$architecture]]
  marker_ids <- colnames(Z_all); sample_ids <- rownames(Z_all)
  if (config$simulation$n_causal > length(marker_ids)) stop("n_causal exceeds the marker count.", call. = FALSE)
  set.seed(spec$simulation_seed)
  causal_index <- sort(sample.int(length(marker_ids), config$simulation$n_causal, replace = FALSE))
  causal_ids <- marker_ids[causal_index]
  if (architecture$effect_distribution == "single_normal") {
    component <- rep("single_normal", length(causal_index))
    raw_effect <- stats::rnorm(length(causal_index))
  } else if (architecture$effect_distribution == "variance_mixture") {
    component_index <- sample.int(length(architecture$mixture_var), length(causal_index),
      replace = TRUE, prob = architecture$mixture_prob)
    component <- paste0("variance_", architecture$mixture_var[component_index])
    raw_effect <- stats::rnorm(length(causal_index), sd = sqrt(architecture$mixture_var[component_index]))
  } else stop("Unknown effect distribution.", call. = FALSE)
  effects <- matrix(0, nrow = length(marker_ids), ncol = 1L,
    dimnames = list(marker_ids, config$trait))
  effects[causal_index, 1L] <- raw_effect
  genetic_values <- Z_all %*% effects
  target_vg <- config$simulation$h2 / (1 - config$simulation$h2)
  effect_scale <- sqrt(target_vg / stats::var(genetic_values[, 1L]))
  effects[, 1L] <- effects[, 1L] * effect_scale
  genetic_values <- Z_all %*% effects
  residual <- stats::rnorm(nrow(Z_all)); residual <- residual - mean(residual)
  residual <- residual / stats::sd(residual)
  residuals <- matrix(residual, ncol = 1L, dimnames = list(sample_ids, config$trait))
  phenotypes <- genetic_values + residuals
  observed_h2 <- stats::var(genetic_values[, 1L]) /
    (stats::var(genetic_values[, 1L]) + stats::var(residuals[, 1L]))
  raw <- list(y = phenotypes, W = Z_all, B = effects, G = genetic_values,
    E = residuals, causal = list(shared = causal_ids,
      specific = stats::setNames(list(character()), config$trait), all = causal_ids),
    rsids = marker_ids, ids = sample_ids, h2_target = config$simulation$h2,
    h2_observed = observed_h2, shared_idx = causal_index,
    specific_idx = stats::setNames(list(integer()), config$trait), causal_rsids = causal_ids)
  sim <- sblrbench::as_sblrbench_simulation(raw, study = config$study,
    architecture = spec$architecture, replicate = as.integer(spec$replicate),
    seed = spec$simulation_seed)
  sim$extras$effect_components <- data.frame(marker = causal_ids,
    component = component, raw_effect = raw_effect,
    final_effect = effects[causal_ids, 1L], stringsAsFactors = FALSE)
  sim$extras$effect_distribution <- architecture$effect_distribution
  sim$extras$effect_scale <- effect_scale
  sblrbench::validate_sblrbench_simulation(sim)
  sim
}

.study02_summary_stats <- function(simulation, Glist, split, config) {
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  stats <- sblr::make_summary_stats(Glist = Glist, y = y, chr = config$chr,
    rows = split$train_rows, scale = TRUE, nthreads = 1L)
  if (!identical(stats$marker_names, simulation$data$marker_ids) ||
      !identical(stats$trait_names, config$trait) ||
      !identical(as.integer(stats$n), length(split$train_ids)) ||
      !isTRUE(all.equal(unname(stats$af[[1L]]),
        unname(Glist$sparseLD$af[[1L]]), tolerance = 0))) {
    stop("Training summary statistics failed sample, marker, trait, or frequency checks.", call. = FALSE)
  }
  stats
}

.study02_method_specs <- function(config) {
  active <- c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr")
  if (!identical(config$methods, active) || any(config$methods %in% config$multitrait$methods)) stop("Active Study 02 methods must be the four configured ST methods.", call. = FALSE)
  map <- list(st_bed_bayesc = c("stblr_bed", "bayesc"),
    st_bed_bayesr = c("stblr_bed", "bayesr"),
    st_csr_sbayesc = c("stblr_csr", "sbayesc"),
    st_csr_sbayesr = c("stblr_csr", "sbayesr"))
  lapply(seq_along(active), function(i) list(id = active[[i]],
    interface = map[[active[[i]]]][[1]], native_method = map[[active[[i]]]][[2]],
    method_index = i, representation = if (grepl("bed", active[[i]])) "BED" else "CSR",
    prior_class = if (grepl("bayesr", active[[i]])) "BayesR" else "BayesC"))
}

.study02_fit <- function(method, simulation, stats, Glist, split, config) {
  seed <- as.integer(config$mcmc$seed_offset +
    match(simulation$scenario$architecture, names(config$simulation$architectures)) * 10000L +
    simulation$scenario$replicate * 100L + method$method_index)
  common <- config$mcmc[c("nit", "nburn", "nthin", "nchains", "ncores", "convergence")]
  common$seed <- seed; common$verbose <- FALSE; common$h2 <- config$priors$h2
  if (method$prior_class == "BayesR") {
    active <- config$priors$bayesr_active_probability
    common$pi <- c(1 - active, rep(active / 3, 3L))
    common$mixture_var <- config$priors$bayesr_mixture_var
  } else common$pi_init <- config$priors$bayesc_inclusion_probability
  spec <- sblrbench::new_sblr_native_method(method$id, method$id,
    method$interface, method$native_method,
    capabilities = c("posterior_effects", "pip", "scalar_trait",
      if (method$representation == "BED") "individual_level" else "summary_statistics"),
    metadata = list(development_settings = TRUE))
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  input <- if (method$representation == "BED") list(y = y, Glist = Glist,
    rows = split$train_rows) else list(stats = stats, Glist = Glist)
  tryCatch({
    result <- sblrbench::run_sblrbench_method(spec, fit_inputs = input, controls = common)
    sblrbench::validate_sblrbench_result(result, simulation)
    list(status = "ok", reason = "", method = method, mcmc_seed = seed,
      controls = common, result = result)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e),
    method = method, mcmc_seed = seed, controls = common, result = NULL))
}

.study02_predict <- function(fit, simulation, test_simulation, Z_test) {
  if (fit$status != "ok") return(fit)
  effects <- align_traits(align_markers(fit$result$estimates$effects,
    simulation$data$marker_ids), simulation$data$trait_names)
  prediction <- Z_test %*% effects
  fit$result <- sblrbench::add_sblrbench_predictions(fit$result, prediction, test_simulation)
  fit
}

.study02_metrics <- function(fit, test_simulation) {
  if (fit$status != "ok") return(data.frame(study = test_simulation$scenario$study,
    scenario = test_simulation$scenario$architecture,
    replicate = test_simulation$scenario$replicate, method_id = fit$method$id,
    trait = test_simulation$data$trait_names, metric = "method_fit",
    value = NA_real_, status = "failed", reason = fit$reason, stringsAsFactors = FALSE))
  sblrbench::evaluate_metrics(test_simulation, fit$result,
    metrics = c("prediction_correlation", "prediction_mse", "prediction_nmse",
      "phenotype_prediction_correlation", "prediction_calibration", "effect_rmse"))
}

.study02_comparisons <- function() data.frame(
  comparison_id = c("bayesr_vs_bayesc_bed", "sbayesr_vs_sbayesc_csr",
    "csr_vs_bed_bayesc", "csr_vs_bed_bayesr"),
  focal_method = c("st_bed_bayesr", "st_csr_sbayesr",
    "st_csr_sbayesc", "st_csr_sbayesr"),
  comparison_method = c("st_bed_bayesc", "st_csr_sbayesc",
    "st_bed_bayesc", "st_bed_bayesr"), stringsAsFactors = FALSE)

.study02_paired <- function(metrics) sblrbench::paired_method_advantages(
  metrics, .study02_comparisons())

.study02_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE); path
}
