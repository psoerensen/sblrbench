.study06_chain_vector <- function(chain, name, nit) {
  x <- chain[[name]]
  if (is.null(x)) return(NULL)
  x <- as.numeric(x)
  if (length(x) < nit) return(NULL)
  tail(x, nit)
}

.study06_extract_draws <- function(run) {
  if (!identical(run$status, "ok"))
    stop("Cannot extract a failed Study 06 fit.", call. = FALSE)
  fit <- run$fit
  values <- fit$convergence_traces$values
  quantities <- fit$convergence_traces$quantities
  nit <- as.integer(fit$input$nit)
  if (length(dim(values)) != 3L || dim(values)[1L] != nit ||
      dim(values)[2L] != 4L || nrow(quantities) != dim(values)[3L])
    stop("Study 06 convergence trace dimensions are invalid.",
      call. = FALSE)
  core <- list()
  for (group in c("vbs", "vgs", "ves", "vle", "vld")) {
    index <- which(quantities$group == group)
    if (length(index) == 1L)
      core[[group]] <- values[, , index, drop = TRUE]
  }
  if (!all(c("vbs", "vgs", "ves") %in% names(core)))
    stop("Required Study 06 variance traces are absent.", call. = FALSE)
  core$heritability <- core$vgs / (core$vgs + core$ves)
  pis <- lapply(fit$chains, .study06_chain_vector,
    name = "pis", nit = nit)
  if (all(vapply(pis, length, integer(1)) == nit))
    core$global_nonnull_proportion <- do.call(cbind, pis)
  out <- do.call(rbind, lapply(names(core), function(name) {
    x <- core[[name]]
    data.frame(architecture = run$architecture,
      replicate = run$replicate,
      configuration = run$method$configuration,
      method = run$method$native_method,
      estimand = name,
      iteration = rep(seq_len(nrow(x)), times = ncol(x)),
      chain = rep(seq_len(ncol(x)), each = nrow(x)),
      value = as.vector(x), stringsAsFactors = FALSE)
  }))
  if (any(!is.finite(out$value)) ||
      !identical(sort(unique(out$chain)), 1:4) ||
      any(table(out$estimand, out$chain) != nit))
    stop("Study 06 retained chain draws are invalid.", call. = FALSE)
  out
}

.study06_parameter_estimates <- function(draws, simulation) {
  groups <- split(draws, draws$estimand)
  available <- do.call(rbind, lapply(groups, function(x) {
    truth <- if (x$estimand[1L] %in% names(simulation$truth))
      unname(simulation$truth[x$estimand[1L]]) else NA_real_
    mean_value <- mean(x$value)
    lower <- unname(stats::quantile(x$value, .025))
    upper <- unname(stats::quantile(x$value, .975))
    data.frame(architecture = x$architecture[1L],
      replicate = x$replicate[1L],
      configuration = x$configuration[1L],
      method = x$method[1L], estimand = x$estimand[1L],
      available = TRUE, unavailable_reason = "",
      posterior_mean = mean_value,
      posterior_sd = stats::sd(x$value),
      posterior_median = stats::median(x$value),
      lower_95 = lower, upper_95 = upper, truth = truth,
      bias = if (is.finite(truth)) mean_value - truth else NA_real_,
      squared_error = if (is.finite(truth))
        (mean_value - truth)^2 else NA_real_,
      interval_coverage = if (is.finite(truth))
        truth >= lower && truth <= upper else NA,
      interpretation = if (is.finite(truth))
        "parameter_recovery" else "posterior_behaviour",
      chain_count = length(unique(x$chain)),
      draws_per_chain = unique(as.integer(table(x$chain)))[1L],
      stringsAsFactors = FALSE)
  }))
  required <- c("effect_variance", "genetic_variance",
    "residual_variance", "heritability", "global_nonnull_proportion",
    "vld", "vle", paste0("component_proportion_", 0:3))
  missing <- setdiff(required, available$estimand)
  if (length(missing)) {
    unavailable <- do.call(rbind, lapply(missing, function(id)
      data.frame(architecture = simulation$architecture,
        replicate = simulation$replicate,
        configuration = draws$configuration[1L],
        method = draws$method[1L], estimand = id, available = FALSE,
        unavailable_reason =
          "not returned as identifiable per-chain draws by installed model",
        posterior_mean = NA_real_, posterior_sd = NA_real_,
        posterior_median = NA_real_, lower_95 = NA_real_,
        upper_95 = NA_real_, truth = NA_real_, bias = NA_real_,
        squared_error = NA_real_, interval_coverage = NA,
        interpretation = "unavailable", chain_count = 4L,
        draws_per_chain = 0L, stringsAsFactors = FALSE)))
    available <- rbind(available, unavailable)
  }
  available[order(available$estimand), , drop = FALSE]
}
