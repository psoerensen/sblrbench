.study06_method_for <- function(architecture, configuration, config) {
  if (!architecture %in% config$architectures ||
      !configuration %in% config$configurations)
    stop("Invalid Study 06 method coordinates.", call. = FALSE)
  bayesr <- identical(architecture, "sparse_mixture")
  native <- if (configuration == "bed")
    if (bayesr) "bayesr" else "bayesc"
  else if (bayesr) "sbayesr" else "sbayesc"
  interface <- if (configuration == "bed") "stblr_bed" else
    if (grepl("^block_eigen", configuration)) "stblr_block_eigen" else
      "stblr_csr"
  data.frame(
    method_id = paste(architecture, configuration, sep = "__"),
    architecture = architecture, configuration = configuration,
    interface = interface, native_method = native,
    model_family = if (bayesr) "BayesR/SBayesR" else "BayesC/SBayesC",
    operator_family = switch(configuration, bed = "packed_bed",
      full_csr = "full_csr", block_csr = "runtime_matched_block_csr",
      "block_eigen"),
    filter_policy = switch(configuration,
      block_eigen_unfiltered = "ridge_fixed_eta_zero",
      block_eigen_hard = "hard_truncate",
      block_eigen_ridge_fixed = "ridge_fixed", "none"),
    stringsAsFactors = FALSE)
}

.study06_method_grid <- function(config) do.call(rbind,
  lapply(config$architectures, function(a) do.call(rbind,
    lapply(config$configurations, function(k)
      .study06_method_for(a, k, config)))))

.study06_baseline_recommendations <- function(config) {
  x <- utils::read.csv(config$recommendation_source,
    stringsAsFactors = FALSE)
  expected <- c("st_bed_bayesc", "st_bed_bayesr",
    "st_csr_sbayesc", "st_csr_sbayesr")
  x <- x[match(expected, x$method), , drop = FALSE]
  if (nrow(x) != 4L || anyNA(x$method) ||
      any(x$recommendation_status != "available") ||
      any(x$recommended_nchains != 4L) ||
      any(x$recommended_ncores != 4L) ||
      any(x$recommended_nthin != 1L) ||
      any(x$recommended_nit_argument !=
        x$recommended_post_burnin_draws))
    stop("Frozen Study 04 recommendation contract failed.",
      call. = FALSE)
  x
}

.study06_baseline_id <- function(architecture, configuration) {
  bed <- identical(configuration, "bed")
  if (architecture == "sparse_homogeneous")
    if (bed) "st_bed_bayesc" else "st_csr_sbayesc"
  else if (bed) "st_bed_bayesr" else "st_csr_sbayesr"
}

.study06_controls <- function(method, config,
                              phase = c("benchmark", "convergence", "pilot"),
                              block_recommendations = NULL) {
  phase <- match.arg(phase)
  if (phase == "convergence") {
    z <- config$convergence
    return(list(nit = z$maximum_nit, nburn = z$maximum_nburn,
      nthin = z$nthin, nchains = z$nchains, ncores = z$ncores))
  }
  if (phase == "pilot") {
    x <- .study06_baseline_recommendations(config)
    id <- .study06_baseline_id(method$architecture, "full_csr")
    z <- x[x$method == id, , drop = FALSE]
    return(list(nit = z$recommended_nit_argument,
      nburn = z$recommended_nburn, nthin = z$recommended_nthin,
      nchains = z$recommended_nchains, ncores = z$recommended_ncores))
  }
  if (method$configuration %in% c("bed", "full_csr")) {
    x <- .study06_baseline_recommendations(config)
    id <- .study06_baseline_id(method$architecture,
      method$configuration)
    z <- x[x$method == id, , drop = FALSE]
    return(list(nit = z$recommended_nit_argument,
      nburn = z$recommended_nburn, nthin = z$recommended_nthin,
      nchains = z$recommended_nchains,
      ncores = z$recommended_ncores))
  }
  if (is.null(block_recommendations))
    stop("Block-operator recommendations are required.", call. = FALSE)
  z <- block_recommendations[
    block_recommendations$architecture == method$architecture &
      block_recommendations$configuration == method$configuration, ,
    drop = FALSE]
  if (nrow(z) != 1L || z$recommendation_status != "available")
    stop("No supported block recommendation for ",
      method$method_id, call. = FALSE)
  list(nit = z$nit, nburn = z$nburn, nthin = z$nthin,
    nchains = z$nchains, ncores = z$ncores)
}

.study06_fit_seeds <- function(method, replicate, config) {
  fit <- .study06_seed(config, method$architecture, replicate,
    method$configuration)
  chain <- vapply(1:4, function(i)
    .study06_seed(config, method$architecture, replicate,
      method$configuration, i), integer(1))
  c(fit_seed = fit, stats::setNames(chain, paste0("chain_", 1:4)))
}

.study06_fit <- function(method, simulation, stats, Glist, split,
                         blocks, runtime_glist, config,
                         phase = c("benchmark", "convergence", "pilot"),
                         block_recommendations = NULL,
                         selected_hard_tau = config$selected_hard_tau,
                         selected_ridge_eta = config$selected_ridge_eta) {
  phase <- match.arg(phase)
  controls <- .study06_controls(method, config, phase,
    block_recommendations)
  seeds <- .study06_fit_seeds(method, simulation$replicate, config)
  controls <- c(controls, list(seed = seeds[["fit_seed"]],
    chain_seeds = unname(seeds[paste0("chain_", 1:4)]),
    keep_chains = TRUE, convergence = "core",
    convergence_control = list(warn = FALSE, keep_traces = TRUE),
    verbose = FALSE, h2 = config$simulation$h2))
  if (grepl("BayesR", method$model_family, fixed = TRUE)) {
    controls$mixture_var <- c(0, 0.01, 0.1, 1)
    controls$pi <- c(0.99, rep(0.01 / 3, 3))
  } else controls$pi_init <- 0.01
  y <- simulation$phenotype[split$train_ids, , drop = FALSE]
  input <- switch(method$configuration,
    bed = list(y = y, Glist = Glist, rows = split$train_rows),
    full_csr = list(stats = stats, Glist = Glist),
    block_csr = list(stats = stats, Glist = runtime_glist),
    block_eigen_unfiltered = list(stats = stats, Glist = Glist,
      block_start = blocks$start),
    block_eigen_hard = list(stats = stats, Glist = Glist,
      block_start = blocks$start),
    block_eigen_ridge_fixed = list(stats = stats, Glist = Glist,
      block_start = blocks$start))
  if (method$configuration == "block_eigen_unfiltered") {
    controls$eigen_filter <- "ridge_fixed"
    controls$eigen_eta <- 0
  } else if (method$configuration == "block_eigen_hard") {
    if (!is.finite(selected_hard_tau) || selected_hard_tau <= 0)
      stop("A nonzero supported hard threshold is required.", call. = FALSE)
    controls$eigen_filter <- "hard_truncate"
    controls$eigen_tau <- selected_hard_tau
  } else if (method$configuration == "block_eigen_ridge_fixed") {
    if (!is.finite(selected_ridge_eta) || selected_ridge_eta <= 0)
      stop("A positive frozen fixed-ridge eta is required.", call. = FALSE)
    controls$eigen_filter <- "ridge_fixed"
    controls$eigen_eta <- selected_ridge_eta
  }
  native <- sblrbench::new_sblr_native_method(method$method_id,
    unname(config$configuration_labels[method$configuration]),
    method$interface, method$native_method,
    capabilities = c("scalar_trait", "multichain",
      "posterior_effects", "pip"))
  started <- Sys.time()
  elapsed <- proc.time()[["elapsed"]]
  tryCatch({
    result <- sblrbench::run_sblrbench_method(native,
      fit_inputs = input, controls = controls)
    fit <- result$native_fit
    if (is.null(fit$convergence_traces$values) ||
        dim(fit$convergence_traces$values)[2L] != 4L ||
        is.null(fit$chains) || length(fit$chains) != 4L ||
        any(vapply(fit$chains, function(x)
          is.null(x$chain_index), logical(1))) ||
        length(unique(vapply(fit$chains, `[[`, integer(1),
          "chain_index"))) != 4L)
      stop("Four identifiable chains were not retained.", call. = FALSE)
    list(status = "ok", reason = "", method = method,
      architecture = simulation$architecture,
      replicate = simulation$replicate, result = result, fit = fit,
      controls = controls, seeds = seeds, started_at = started,
      finished_at = Sys.time(),
      runtime = proc.time()[["elapsed"]] - elapsed,
      warnings = paste(result$diagnostics$warnings, collapse = " | "))
  }, error = function(e) list(status = "failed",
    reason = conditionMessage(e), method = method,
    architecture = simulation$architecture,
    replicate = simulation$replicate, result = NULL, fit = NULL,
    controls = controls, seeds = seeds, started_at = started,
    finished_at = Sys.time(),
    runtime = proc.time()[["elapsed"]] - elapsed, warnings = ""))
}
