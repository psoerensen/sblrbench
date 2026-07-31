.study03_method_specs <- function(config) {
  map <- list(st_bed_bayesc = c("stblr_bed", "bayesc"),
    st_bed_bayesr = c("stblr_bed", "bayesr"),
    st_csr_sbayesc = c("stblr_csr", "sbayesc"),
    st_csr_sbayesr = c("stblr_csr", "sbayesr"))
  if (!identical(config$methods, names(map))) stop("Study 03 requires exactly four ST methods.", call. = FALSE)
  lapply(seq_along(map), function(i) list(id = names(map)[i], interface = map[[i]][1],
    native_method = map[[i]][2], method_index = i,
    representation = if (grepl("bed", names(map)[i])) "BED" else "CSR",
    prior_class = if (grepl("bayesr", names(map)[i])) "BayesR" else "BayesC"))
}

.study03_seed_registry <- function(specs, methods, config) {
  do.call(rbind, lapply(specs, function(s) do.call(rbind, lapply(methods, function(m) {
    fit_seed <- as.integer(config$seeds$method_offset +
      match(s$architecture, names(config$simulation$architectures)) * 10000L +
      s$replicate * 100L + m$method_index)
    nchains <- if (identical(config$profile, "five_replicate_development")) 4L else 1L
    data.frame(architecture = s$architecture, replicate = s$replicate,
      data_selection_seed = config$seeds$data_selection,
      architecture_seed = config$simulation$base_seed +
        match(s$architecture, names(config$simulation$architectures)) * 1000L,
      simulation_seed = s$simulation_seed,
      method = m$id, fit_seed = fit_seed, chain = seq_len(nchains),
      chain_seed = if (nchains == 1L) fit_seed else .five_replicate_chain_seeds(fit_seed, nchains),
      stringsAsFactors = FALSE)
  }))))
}

.study03_fit <- function(method, simulation, stats, Glist, config) {
  profile <- config$profiles[[config$profile]]
  seed <- as.integer(config$seeds$method_offset +
    match(simulation$scenario$architecture, names(config$simulation$architectures)) * 10000L +
    simulation$scenario$replicate * 100L + method$method_index)
  controls <- if (identical(config$profile, "five_replicate_development"))
    .five_replicate_mcmc(method$id) else
    profile[c("nit", "nburn", "nthin", "nchains", "ncores", "convergence")]
  controls$seed <- seed; controls$verbose <- FALSE; controls$h2 <- config$priors$h2
  chain_seeds <- if (controls$nchains == 1L) seed else
    .five_replicate_chain_seeds(seed, controls$nchains)
  if (controls$nchains > 1L) controls$chain_seeds <- chain_seeds
  if (method$prior_class == "BayesR") {
    p <- config$priors$bayesr_active_probability
    controls$pi <- c(1 - p, rep(p / 3, 3)); controls$mixture_var <- config$priors$bayesr_mixture_var
  } else controls$pi_init <- config$priors$bayesc_inclusion_probability
  native <- sblrbench::new_sblr_native_method(method$id, method$id, method$interface,
    method$native_method, capabilities = c("posterior_effects", "pip", "scalar_trait",
      if (method$representation == "BED") "individual_level" else "summary_statistics"))
  input <- if (method$representation == "BED") list(y = simulation$truth$phenotypes,
    Glist = Glist, rows = seq_len(nrow(simulation$truth$phenotypes))) else list(stats = stats, Glist = Glist)
  tryCatch({
    result <- sblrbench::run_sblrbench_method(native, input, controls)
    list(status = "ok", reason = "", method = method, seed = seed, controls = controls,
      result = result, native_fit = result$native_fit, chain_seeds = chain_seeds)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e), method = method,
    seed = seed, controls = controls, result = NULL, native_fit = NULL,
    chain_seeds = chain_seeds))
}

.study03_trace_vector <- function(x, name) {
  if (is.null(x)) stop("Missing posterior component: ", name, call. = FALSE)
  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || ncol(x) != 1L) stop("Study 03 requires one-column trace component: ", name, call. = FALSE)
  as.numeric(x[, 1L])
}

.study03_extract_draws <- function(fit, method, registry, marker_count) {
  if (inherits(fit, "sblrbench_result")) fit <- fit$native_fit
  if (is.null(fit$input) || is.null(fit$input$nburn)) stop("Fit does not record burn-in.", call. = FALSE)
  required <- c("vbs", "vgs", "ves")
  components <- lapply(required, function(n) .study03_trace_vector(fit[[n]], n))
  names(components) <- required
  if (!is.null(fit$pi_trace)) components$pi_trace <- .study03_trace_vector(fit$pi_trace, "pi_trace")
  lengths <- vapply(components, length, integer(1)); if (length(unique(lengths)) != 1L) stop("Posterior traces are not jointly aligned.", call. = FALSE)
  nburn <- as.integer(fit$input$nburn); nit <- as.integer(fit$input$nit); nthin <- as.integer(fit$input$nthin %||% 1L)
  if (lengths[1] == nburn + nit) keep <- seq.int(nburn + 1L, nburn + nit, by = nthin)
  else if (lengths[1] == nit && nburn < nit) keep <- seq.int(nburn + 1L, nit, by = nthin)
  else stop("Posterior trace length is inconsistent with fit metadata.", call. = FALSE)
  z <- lapply(components, `[`, keep)
  values <- list(effect_variance = z$vbs, genetic_variance = z$vgs, residual_variance = z$ves,
    heritability = z$vgs / (z$vgs + z$ves))
  if (!is.null(z$pi_trace)) values <- c(list(causal_proportion = z$pi_trace,
    total_marker_effect_variance = z$vbs * z$pi_trace * marker_count), values)
  do.call(rbind, lapply(names(values), function(id) data.frame(method = method,
    estimand_id = id, draw = seq_along(values[[id]]), value = values[[id]],
    source_component = registry$posterior_source[match(id, registry$estimand_id)],
    status = "ok", reason = "", stringsAsFactors = FALSE)))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.study03_summarise_draws <- function(draws) {
  do.call(rbind, lapply(split(draws, draws$estimand_id), function(x) data.frame(
    method = x$method[1], estimand_id = x$estimand_id[1], posterior_mean = mean(x$value),
    posterior_sd = stats::sd(x$value), posterior_median = stats::median(x$value),
    lower_95 = unname(stats::quantile(x$value, .025)), upper_95 = unname(stats::quantile(x$value, .975)),
    n_posterior_draws = nrow(x), chain_count = length(unique(x$chain %||% 1L)),
    draws_per_chain = if ("chain" %in% names(x)) unique(as.integer(table(x$chain)))[1L] else nrow(x),
    mcse_mean = NA_real_, status = "ok", reason = "",
    stringsAsFactors = FALSE)))
}

.study03_extract_multichain_draws <- function(fit, method, registry, marker_count,
                                               expected_chains = 4L) {
  bundle <- fit$convergence_traces
  if (is.null(bundle$values) || length(dim(bundle$values)) != 3L)
    stop("Native convergence trace bundle is unavailable.", call. = FALSE)
  values <- bundle$values
  if (dim(values)[2L] != expected_chains || dim(values)[1L] != fit$input$nit)
    stop("Multichain trace dimensions disagree with retained-draw metadata.", call. = FALSE)
  groups <- as.character(bundle$quantities$group)
  idx <- match(c("vbs", "vgs", "ves"), groups)
  if (anyNA(idx)) stop("Required scalar traces are absent.", call. = FALSE)
  base <- expand.grid(iteration = seq_len(dim(values)[1L]), chain = seq_len(expected_chains))
  base$vbs <- as.vector(values[, , idx[1L]])
  base$vgs <- as.vector(values[, , idx[2L]])
  base$ves <- as.vector(values[, , idx[3L]])
  if (any(!is.finite(as.matrix(base[c("vbs", "vgs", "ves")]))) ||
      any(as.matrix(base[c("vbs", "vgs", "ves")]) < 0))
    stop("Core multichain draws are invalid.", call. = FALSE)
  values_by_estimand <- list(effect_variance = base$vbs,
    genetic_variance = base$vgs, residual_variance = base$ves,
    heritability = base$vgs / (base$vgs + base$ves))
  pi_values <- NULL
  if (!is.null(fit$chains) && length(fit$chains) == expected_chains) {
    pi_by_chain <- lapply(fit$chains, function(ch) {
      z <- if (is.null(ch$pi_trace)) NULL else as.numeric(ch$pi_trace)
      if (is.null(z)) return(NULL)
      if (length(z) == fit$input$nit) z else if (length(z) == fit$input$nit + fit$input$nburn)
        tail(z, fit$input$nit) else NULL
    })
    if (all(vapply(pi_by_chain, length, integer(1)) == fit$input$nit))
      pi_values <- unlist(pi_by_chain, use.names = FALSE)
  }
  if (!is.null(pi_values)) {
    values_by_estimand$causal_proportion <- pi_values
    values_by_estimand$total_marker_effect_variance <- base$vbs * pi_values * marker_count
  }
  out <- do.call(rbind, lapply(names(values_by_estimand), function(id) data.frame(
    method = method, estimand_id = id, chain = base$chain, iteration = base$iteration,
    value = values_by_estimand[[id]],
    source_component = registry$posterior_source[match(id, registry$estimand_id)],
    status = "ok", reason = "", stringsAsFactors = FALSE)))
  counts <- table(out$estimand_id, out$chain)
  if (any(counts != fit$input$nit)) stop("Retained per-chain draw counts are unequal.", call. = FALSE)
  out
}
