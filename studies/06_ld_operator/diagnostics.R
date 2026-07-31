.study06_diagnostic_one <- function(x, thresholds) {
  wide <- reshape(x[c("iteration", "chain", "value")],
    idvar = "iteration", timevar = "chain", direction = "wide")
  m <- as.matrix(wide[-1L])
  sd_all <- stats::sd(as.numeric(m))
  rhat <- posterior::rhat(m)
  bulk <- posterior::ess_bulk(m)
  tail <- posterior::ess_tail(m)
  mcse <- posterior::mcse_mean(m)
  relative <- if (is.finite(sd_all) && sd_all > 0)
    mcse / sd_all else NA_real_
  flags <- c(is.finite(rhat) && rhat <= thresholds$rhat,
    is.finite(bulk) && bulk >= thresholds$ess_bulk,
    is.finite(tail) && tail >= thresholds$ess_tail,
    is.finite(relative) && relative <= thresholds$relative_mcse,
    ncol(m) == thresholds$chain_count,
    all(is.finite(m)) && is.finite(sd_all) && sd_all > 0)
  data.frame(rhat = rhat, ess_bulk = bulk, ess_tail = tail,
    mcse_mean = mcse, posterior_sd = sd_all,
    relative_mcse = relative, chain_count = ncol(m),
    draws_per_chain = nrow(m),
    rhat_pass = flags[1L], ess_bulk_pass = flags[2L],
    ess_tail_pass = flags[3L], relative_mcse_pass = flags[4L],
    chain_count_pass = flags[5L], finite_nonconstant_pass = flags[6L],
    overall_pass = all(flags),
    limiting_diagnostic = c("rhat", "ess_bulk", "ess_tail",
      "relative_mcse", "chain_count", "finite_nonconstant")[
        which(!flags)[1L] %||% 1L],
    stringsAsFactors = FALSE)
}

.study06_diagnostics <- function(draws, burnin, retained, config) {
  keep <- draws$iteration > burnin &
    draws$iteration <= burnin + retained
  z <- draws[keep, , drop = FALSE]
  if (!nrow(z)) stop("Empty Study 06 diagnostic window.", call. = FALSE)
  groups <- split(z, z$estimand)
  do.call(rbind, lapply(groups, function(x) {
    d <- .study06_diagnostic_one(x,
      config$convergence$thresholds)
    cbind(x[1L, c("architecture", "replicate",
      "configuration", "method", "estimand")],
      burnin_candidate = burnin,
      retained_draw_candidate = retained, d,
      stringsAsFactors = FALSE)
  }))
}

.study06_select_recommendation <- function(draws, config) {
  grid <- expand.grid(
    nburn = config$convergence$burnin_candidates,
    nit = config$convergence$retained_candidates,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[grid$nburn + grid$nit <=
      config$convergence$maximum_nit, , drop = FALSE]
  grid <- grid[order(grid$nburn + grid$nit, grid$nit,
    grid$nburn), , drop = FALSE]
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(grid)),
    function(i) .study06_diagnostics(draws, grid$nburn[i],
      grid$nit[i], config)))
  key <- interaction(diagnostics$burnin_candidate,
    diagnostics$retained_draw_candidate, drop = TRUE)
  candidates <- do.call(rbind, lapply(split(diagnostics, key),
    function(x) data.frame(
      architecture = x$architecture[1L],
      configuration = x$configuration[1L],
      nburn = x$burnin_candidate[1L],
      nit = x$retained_draw_candidate[1L],
      failed_estimands = sum(!x$overall_pass),
      maximum_rhat = max(x$rhat),
      minimum_bulk_ess = min(x$ess_bulk),
      minimum_tail_ess = min(x$ess_tail),
      maximum_relative_mcse = max(x$relative_mcse),
      limiting_estimand = x$estimand[which.max(pmax(
        x$rhat / config$convergence$thresholds$rhat,
        config$convergence$thresholds$ess_bulk / x$ess_bulk,
        config$convergence$thresholds$ess_tail / x$ess_tail,
        x$relative_mcse /
          config$convergence$thresholds$relative_mcse))],
      pass = all(x$overall_pass), stringsAsFactors = FALSE)))
  candidates <- candidates[order(candidates$nburn + candidates$nit,
    candidates$nit, candidates$nburn), , drop = FALSE]
  passed <- candidates[candidates$pass, , drop = FALSE]
  if (!nrow(passed))
    stop("No supported Study 06 convergence setting for ",
      draws$architecture[1L], " / ", draws$configuration[1L],
      call. = FALSE)
  selected <- passed[1L, , drop = FALSE]
  recommendation <- data.frame(
    architecture = selected$architecture,
    configuration = selected$configuration,
    recommendation_status = "available",
    nburn = selected$nburn, nit = selected$nit,
    nthin = config$convergence$nthin,
    nchains = config$convergence$nchains,
    ncores = config$convergence$ncores,
    limiting_estimand = selected$limiting_estimand,
    maximum_rhat = selected$maximum_rhat,
    minimum_bulk_ess = selected$minimum_bulk_ess,
    minimum_tail_ess = selected$minimum_tail_ess,
    maximum_relative_mcse = selected$maximum_relative_mcse,
    reason = "earliest stable operational pass",
    stringsAsFactors = FALSE)
  list(diagnostics = diagnostics, candidates = candidates,
    recommendation = recommendation)
}

`%||%` <- function(x, y) {
  missing_scalar <- length(x) == 1L && is.na(x)
  if (is.null(x) || !length(x) || missing_scalar) y else x
}
