.study06v2_method_for <- function(architecture, configuration, config) {
  stopifnot(architecture %in% config$architectures,
    configuration %in% config$configurations)
  bayesr <- identical(architecture, "sparse_mixture")
  low_rank <- startsWith(configuration, "low_rank_")
  data.frame(
    method_id = paste(architecture, configuration, sep = "__"),
    architecture = architecture,
    configuration = configuration,
    interface = if (configuration == "bed") "stblr_bed" else
      if (low_rank) "stblr_block_eigen" else "stblr_csr",
    native_method = if (configuration == "bed") {
      if (bayesr) "bayesr" else "bayesc"
    } else if (bayesr) "sbayesr" else "sbayesc",
    model_family = if (bayesr) "BayesR/SBayesR" else "BayesC/SBayesC",
    operator_family = if (configuration == "bed") "packed_bed" else
      if (configuration == "full_csr") "full_csr" else
        if (configuration == "block_csr") "block_csr" else
          "retained_low_rank",
    representation = if (low_rank) "low_rank" else NA_character_,
    eigen_prop = if (configuration == "low_rank_full")
      config$eigen_prop_full else if (configuration == "low_rank_0999")
        config$eigen_props[["low_rank_0999"]] else if (
          configuration == "low_rank_0995")
            config$eigen_props[["low_rank_0995"]] else NA_real_,
    stringsAsFactors = FALSE)
}

.study06v2_method_grid <- function(config) do.call(rbind,
  lapply(config$architectures, function(a) do.call(rbind,
    lapply(config$configurations, function(k)
      .study06v2_method_for(a, k, config)))))

.study06v2_seed <- function(config, architecture, replicate,
                            configuration = NULL, chain = NULL) {
  a <- match(architecture, config$architectures)
  k <- if (is.null(configuration)) 0L else
    match(configuration, config$configurations)
  if (is.na(a) || is.na(k) || !replicate %in% seq_len(config$replicate_count))
    stop("Invalid Study 06 v2 seed coordinates.", call. = FALSE)
  value <- config$seeds$simulation_base +
    (a - 1L) * config$seeds$architecture_stride +
    as.integer(replicate) * config$seeds$replicate_stride
  if (!is.null(configuration)) value <- value +
    config$seeds$fit_base + k * config$seeds$configuration_stride
  if (!is.null(chain)) {
    if (!chain %in% 1:4) stop("Chain must be 1--4.", call. = FALSE)
    value <- value + as.integer(chain) * config$seeds$chain_stride
  }
  as.integer(value %% .Machine$integer.max)
}

.study06v2_baseline_controls <- function(method, config) {
  source("studies/06_ld_operator/methods.R", local = TRUE)
  x <- .study06_baseline_recommendations(config)
  id <- .study06_baseline_id(method$architecture,
    if (method$configuration == "bed") "bed" else "full_csr")
  z <- x[x$method == id, , drop = FALSE]
  list(nit = z$recommended_nit_argument, nburn = z$recommended_nburn,
    nthin = z$recommended_nthin, nchains = z$recommended_nchains,
    ncores = z$recommended_ncores)
}

.study06v2_controls <- function(method, config, phase,
                                recommendations = NULL) {
  if (phase == "operator-pilot")
    return(list(nit = 50L, nburn = 25L, nthin = 1L,
      nchains = 4L, ncores = 4L))
  if (phase == "convergence")
    return(list(nit = config$convergence$maximum_nit,
      nburn = 0L, nthin = 1L, nchains = 4L, ncores = 4L))
  if (method$configuration %in% c("bed", "full_csr"))
    return(.study06v2_baseline_controls(method, config))
  if (is.null(recommendations))
    stop("Validated v2 convergence recommendations are required.",
      call. = FALSE)
  z <- recommendations[recommendations$architecture == method$architecture &
    recommendations$configuration == method$configuration, , drop = FALSE]
  if (nrow(z) != 1L || !identical(z$recommendation_status, "available"))
    stop("No supported v2 recommendation for ", method$method_id,
      call. = FALSE)
  list(nit = z$nit, nburn = z$nburn, nthin = z$nthin,
    nchains = z$nchains, ncores = z$ncores)
}

.study06v2_input_hash <- function(simulation, stats, method, controls,
                                  config) digest::digest(list(
  architecture = simulation$architecture,
  replicate = simulation$replicate,
  simulation_seed = simulation$simulation_seed,
  marker_ids = simulation$marker_ids,
  effects = simulation$effects,
  phenotype = simulation$phenotype,
  stats_marker_names = stats$marker_names,
  stats_rows = stats$rows,
  stats_wy = stats$wy,
  method = method,
  controls = controls,
  block_size = config$block_size), algo = "sha256")

.study06v2_fit <- function(method, simulation, stats, Glist, split, blocks,
                           block_glist, config, phase,
                           recommendations = NULL) {
  controls <- .study06v2_controls(method, config, phase, recommendations)
  fit_seed <- .study06v2_seed(config, method$architecture,
    simulation$replicate, method$configuration)
  chain_seeds <- vapply(1:4, function(chain)
    .study06v2_seed(config, method$architecture, simulation$replicate,
      method$configuration, chain), integer(1))
  controls <- c(controls, list(seed = fit_seed,
    chain_seeds = unname(chain_seeds), keep_chains = TRUE,
    convergence = "core",
    convergence_control = list(warn = FALSE, keep_traces = TRUE),
    verbose = FALSE, h2 = config$simulation$h2))
  if (grepl("BayesR", method$model_family, fixed = TRUE)) {
    controls$mixture_var <- c(0, 0.01, 0.1, 1)
    controls$pi <- c(0.99, rep(0.01 / 3, 3))
  } else controls$pi_init <- 0.01
  low_rank <- identical(method$operator_family, "retained_low_rank")
  if (low_rank) {
    controls$representation <- "low_rank"
    controls$eigen_prop <- as.numeric(method$eigen_prop)
  }
  .study06v2_assert_fit_spec(method$configuration, controls, config)
  input <- if (method$configuration == "bed") {
    list(y = simulation$phenotype[split$train_ids, , drop = FALSE],
      Glist = Glist, rows = split$train_rows)
  } else if (method$configuration == "full_csr") {
    list(stats = stats, Glist = Glist)
  } else if (method$configuration == "block_csr") {
    list(stats = stats, Glist = block_glist)
  } else list(stats = stats, Glist = Glist,
    block_start = blocks$start)
  native <- sblrbench::new_sblr_native_method(method$method_id,
    unname(config$configuration_labels[method$configuration]),
    method$interface, method$native_method,
    capabilities = c("scalar_trait", "multichain",
      "posterior_effects", "pip"))
  started <- Sys.time(); elapsed <- proc.time()[["elapsed"]]
  warning_text <- character()
  result <- tryCatch(withCallingHandlers(
    sblrbench::run_sblrbench_method(native, fit_inputs = input,
      controls = controls), warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = identity)
  runtime <- proc.time()[["elapsed"]] - elapsed
  if (inherits(result, "error")) return(list(status = "failed",
    reason = conditionMessage(result), method = method,
    architecture = simulation$architecture,
    replicate = simulation$replicate, controls = controls,
    result = NULL, fit = NULL, started_at = started,
    finished_at = Sys.time(), runtime = runtime,
    warnings = paste(unique(warning_text), collapse = " | ")))
  fit <- result$native_fit
  valid_chains <- !is.null(fit$chains) && length(fit$chains) == 4L &&
    identical(sort(vapply(fit$chains, `[[`, integer(1), "chain_index")), 1:4)
  if (!valid_chains) return(list(status = "failed",
    reason = "Four identifiable chains were not retained.", method = method,
    architecture = simulation$architecture,
    replicate = simulation$replicate, controls = controls,
    result = result, fit = fit, started_at = started,
    finished_at = Sys.time(), runtime = runtime,
    warnings = paste(unique(warning_text), collapse = " | ")))
  if (low_rank && (!identical(fit$input$operator_contract,
      config$operator_contract) ||
      !identical(fit$input$operator_representation, "low_rank") ||
      !isTRUE(all.equal(fit$input$eigen_prop, method$eigen_prop))))
    return(list(status = "failed",
      reason = "Low-rank fit metadata failed the v2 contract.", method = method,
      architecture = simulation$architecture,
      replicate = simulation$replicate, controls = controls,
      result = result, fit = fit, started_at = started,
      finished_at = Sys.time(), runtime = runtime,
      warnings = paste(unique(warning_text), collapse = " | ")))
  list(status = "ok", reason = "", method = method,
    architecture = simulation$architecture,
    replicate = simulation$replicate, controls = controls,
    seeds = c(fit_seed = fit_seed, stats::setNames(chain_seeds,
      paste0("chain_", 1:4))), result = result, fit = fit,
    source_sha = config$required_sblr_sha,
    package_version = config$required_sblr_version,
    operator_contract = if (low_rank) config$operator_contract else NA_character_,
    representation = if (low_rank) "low_rank" else NA_character_,
    eigen_prop = if (low_rank) method$eigen_prop else NA_real_,
    started_at = started, finished_at = Sys.time(), runtime = runtime,
    warnings = paste(unique(c(warning_text,
      result$diagnostics$warnings)), collapse = " | "))
}
