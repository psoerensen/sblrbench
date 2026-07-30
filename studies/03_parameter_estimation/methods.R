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
  do.call(rbind, lapply(specs, function(s) do.call(rbind, lapply(methods, function(m)
    data.frame(architecture = s$architecture, replicate = s$replicate,
      data_selection_seed = config$seeds$data_selection,
      simulation_seed = s$simulation_seed,
      method = m$id, fit_seed = as.integer(config$seeds$method_offset +
        match(s$architecture, names(config$simulation$architectures)) * 10000L +
        s$replicate * 100L + m$method_index), stringsAsFactors = FALSE)))))
}

.study03_fit <- function(method, simulation, stats, Glist, config) {
  profile <- config$profiles[[config$profile]]
  seed <- as.integer(config$seeds$method_offset +
    match(simulation$scenario$architecture, names(config$simulation$architectures)) * 10000L +
    simulation$scenario$replicate * 100L + method$method_index)
  controls <- profile[c("nit", "nburn", "nthin", "nchains", "ncores", "convergence")]
  controls$seed <- seed; controls$verbose <- FALSE; controls$h2 <- config$priors$h2
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
      result = result, native_fit = result$native_fit)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e), method = method,
    seed = seed, controls = controls, result = NULL, native_fit = NULL))
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
    n_posterior_draws = nrow(x), mcse_mean = NA_real_, status = "ok", reason = "",
    stringsAsFactors = FALSE)))
}
