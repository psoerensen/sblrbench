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
    if (identical(ids, marker_ids) && identical(x$sparseLD$rows, split$train_rows) && length(Sys.glob(paste0(x$sparseLD$prefix, "*")))) return(x)
  }
  x <- do.call(sblr::make_sparse_ld, c(list(Glist = Glist, rows = split$train_rows,
    out_prefix = prefix, chr = config$chr), config$sparse_ld))
  ids <- x$rsids[[config$chr]][x$sparseLD$cls[[1L]]]
  if (!identical(ids, marker_ids)) stop("Training sparse-LD marker order is invalid.", call. = FALSE)
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
    simulation_seed = as.integer(config$simulation$base_seed + match(a, names(config$simulation$architectures)) * 1000L + i))
  out
}

.study02_simulate <- function(spec, Z_all, config) {
  a <- config$simulation$architectures[[spec$architecture]]
  # mtsim 0.1.2 requires an un-named unit diagonal in the correlation matrix.
  rg <- matrix(config$simulation$shared_effect_correlation,
               length(config$traits), length(config$traits)); diag(rg) <- 1
  raw <- sblr::mtsim(W = Z_all, standardize_W = FALSE, nt = length(config$traits),
    n_shared = a$n_shared, n_specific = a$n_specific, h2 = config$simulation$h2,
    rg = rg, re = config$simulation$residual_correlation, seed = spec$simulation_seed)
  colnames(raw$y) <- colnames(raw$B) <- colnames(raw$G) <- colnames(raw$E) <- config$traits
  sim <- sblrbench::as_sblrbench_simulation(raw, study = config$study,
    architecture = spec$architecture, replicate = as.integer(spec$replicate), seed = spec$simulation_seed)
  if (length(sim$truth$causal$shared) != a$n_shared ||
      any(lengths(sim$truth$causal$specific) != a$n_specific)) stop("mtsim causal counts do not match the architecture.", call. = FALSE)
  sim
}

.study02_summary_stats <- function(simulation, Glist, split, config) {
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  full <- sblr::make_summary_stats(Glist = Glist, y = y, chr = config$chr,
    rows = split$train_rows, scale = TRUE, nthreads = 1L)
  if (!identical(full$marker_names, simulation$data$marker_ids) ||
      !identical(full$trait_names, config$traits) || !identical(as.integer(full$n), length(split$train_ids))) stop("Training summary statistics failed alignment checks.", call. = FALSE)
  single <- lapply(config$traits, function(t) sblr::make_summary_stats(
    Glist = Glist, y = y[, t], chr = config$chr, rows = split$train_rows,
    scale = TRUE, nthreads = 1L))
  names(single) <- config$traits
  list(multi = full, single = single)
}

.study02_method_specs <- function(config) {
  map <- list(st_bed_bayesr = c("stblr_bed", "bayesr"), mt_bed_bayesr = c("mtblr_bed", "bayesr"),
    st_csr_sbayesr = c("stblr_csr", "sbayesr"), mt_csr_sbayesr = c("mtblr_csr", "sbayesr"))
  lapply(seq_along(config$methods), function(i) list(id = config$methods[[i]],
    interface = map[[config$methods[[i]]]][[1]], native_method = map[[config$methods[[i]]]][[2]], method_index = i,
    representation = if (grepl("bed", config$methods[[i]])) "BED" else "CSR",
    model_scope = if (startsWith(config$methods[[i]], "mt_")) "MT" else "ST"))
}

.study02_one_fit <- function(method, fit_inputs, controls) {
  spec <- sblrbench::new_sblr_native_method(method$id, method$id, method$interface,
    method$native_method, capabilities = c("posterior_effects", "pip", if (method$model_scope == "MT") "multi_trait" else "scalar_trait",
      if (method$representation == "BED") "individual_level" else "summary_statistics"),
    metadata = list(development_settings = TRUE))
  sblrbench::run_sblrbench_method(spec, fit_inputs = fit_inputs, controls = controls)
}

.study02_fit <- function(method, simulation, stats, Glist, split, config) {
  seed <- as.integer(config$mcmc$seed_offset + match(simulation$scenario$architecture, names(config$simulation$architectures)) * 10000L +
    simulation$scenario$replicate * 100L + method$method_index)
  common <- config$mcmc[c("nit", "nburn", "nthin", "nchains", "ncores", "convergence")]
  common$verbose <- FALSE
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  tryCatch({
    if (method$model_scope == "ST") {
      runs <- vector("list", length(config$traits)); names(runs) <- config$traits
      for (j in seq_along(config$traits)) {
        input <- if (method$representation == "BED") list(y = y[, j, drop = FALSE], Glist = Glist,
          rows = split$train_rows) else list(stats = stats$single[[j]], Glist = Glist)
        active <- config$bayesr$active_probability
        st_pi <- c(1 - active, rep(active / (length(config$bayesr$mixture_var) - 1L),
          length(config$bayesr$mixture_var) - 1L))
        controls <- c(common, list(h2 = config$bayesr$h2[[j]], pi = st_pi,
          mixture_var = config$bayesr$mixture_var, seed = seed + j))
        runs[[j]] <- .study02_one_fit(method, input, controls)
      }
      extracted <- runs
      effects <- do.call(cbind, lapply(extracted, function(x) x$estimates$effects[, 1L])); colnames(effects) <- config$traits
      pip <- do.call(cbind, lapply(extracted, function(x) x$estimates$pip[, 1L])); colnames(pip) <- config$traits
      result <- sblrbench::new_sblrbench_result(method$id, effects = effects, pip = pip,
        convergence = lapply(extracted, function(x) x$diagnostics$convergence),
        warnings = unlist(lapply(extracted, function(x) x$diagnostics$warnings)),
        elapsed_seconds = sum(vapply(extracted, function(x) x$computation$elapsed_seconds, numeric(1))),
        provenance = list(trait_seeds = seed + seq_along(config$traits)), native_fit = runs)
    } else {
      input <- if (method$representation == "BED") list(y = y, Glist = Glist,
        rows = split$train_rows, residual_covariance = "diagonal") else list(stats = stats$multi,
          Glist = Glist, sample_overlap = "not_modeled")
      controls <- c(common, list(h2 = config$bayesr$h2, pi = config$bayesr$active_probability,
        mixture_var = config$bayesr$mixture_var, seed = seed))
      result <- .study02_one_fit(method, input, controls)
    }
    sblrbench::validate_sblrbench_result(result, simulation)
    list(status = "ok", reason = "", method = method, mcmc_seed = seed,
         controls = common, result = result)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e), method = method,
    mcmc_seed = seed, controls = common, result = NULL))
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
    scenario = test_simulation$scenario$architecture, replicate = test_simulation$scenario$replicate,
    method_id = fit$method$id, trait = rep(test_simulation$data$trait_names, each = 1L), metric = "method_fit",
    value = NA_real_, status = "failed", reason = fit$reason, stringsAsFactors = FALSE))
  sblrbench::evaluate_metrics(test_simulation, fit$result, metrics = c("prediction_correlation",
    "prediction_mse", "prediction_nmse", "phenotype_prediction_correlation",
    "prediction_calibration", "effect_rmse"))
}

.study02_paired <- function(metrics) {
  map <- data.frame(method = c("st_bed_bayesr", "mt_bed_bayesr", "st_csr_sbayesr", "mt_csr_sbayesr"),
    representation = rep(c("BED", "CSR"), each = 2), scope = rep(c("ST", "MT"), 2), stringsAsFactors = FALSE)
  sblrbench::paired_mt_advantages(metrics, map)
}

.study02_write_csv <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); utils::write.csv(x, path, row.names = FALSE); path }
