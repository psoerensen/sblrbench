.study07_diagnostic_one <- function(x, thresholds) {
  wide <- reshape(x[c("iteration", "chain", "value")],
    idvar = "iteration", timevar = "chain", direction = "wide")
  m <- as.matrix(wide[-1L]); sd_all <- stats::sd(as.numeric(m))
  rhat <- posterior::rhat(m); bulk <- posterior::ess_bulk(m)
  tail <- posterior::ess_tail(m); mcse <- posterior::mcse_mean(m)
  relative <- if (is.finite(sd_all) && sd_all > 0) mcse / sd_all else NA_real_
  flags <- c(is.finite(rhat) && rhat <= thresholds$rhat,
    is.finite(bulk) && bulk >= thresholds$ess_bulk,
    is.finite(tail) && tail >= thresholds$ess_tail,
    is.finite(relative) && relative <= thresholds$relative_mcse,
    ncol(m) == thresholds$chain_count,
    all(is.finite(m)) && is.finite(sd_all) && sd_all > 0)
  data.frame(rhat = rhat, ess_bulk = bulk, ess_tail = tail,
    mcse_mean = mcse, posterior_sd = sd_all, relative_mcse = relative,
    chain_count = ncol(m), draws_per_chain = nrow(m),
    rhat_pass = flags[[1L]], ess_bulk_pass = flags[[2L]],
    ess_tail_pass = flags[[3L]], relative_mcse_pass = flags[[4L]],
    chain_count_pass = flags[[5L]], finite_nonconstant_pass = flags[[6L]],
    overall_pass = all(flags),
    limiting_diagnostic = if (all(flags)) "none" else
      c("rhat", "ess_bulk", "ess_tail", "relative_mcse",
        "chain_count", "finite_nonconstant")[which(!flags)[[1L]]],
    stringsAsFactors = FALSE)
}

.study07_diagnostics <- function(draws, burnin, retained, config) {
  z <- draws[draws$iteration > burnin &
    draws$iteration <= burnin + retained, , drop = FALSE]
  if (!nrow(z)) stop("Empty Study 07 diagnostic window.", call. = FALSE)
  do.call(rbind, lapply(split(z, z$estimand), function(x) cbind(
    x[1L, c("architecture", "replicate", "implementation", "estimand")],
    nburn = burnin, nit = retained,
    .study07_diagnostic_one(x, config$convergence$thresholds),
    stringsAsFactors = FALSE)))
}

.study07_select_recommendation <- function(draws, config) {
  grid <- expand.grid(nburn = config$convergence$burnin_candidates,
    nit = config$convergence$retained_candidates,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- grid[grid$nburn + grid$nit <= config$convergence$maximum_nit, ]
  grid <- grid[order(grid$nburn + grid$nit, grid$nit, grid$nburn), ]
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    .study07_diagnostics(draws, grid$nburn[[i]], grid$nit[[i]], config)))
  candidates <- do.call(rbind, lapply(split(diagnostics,
    interaction(diagnostics$nburn, diagnostics$nit, drop = TRUE)), function(x) {
      score <- pmax(x$rhat / config$convergence$thresholds$rhat,
        config$convergence$thresholds$ess_bulk / x$ess_bulk,
        config$convergence$thresholds$ess_tail / x$ess_tail,
        x$relative_mcse / config$convergence$thresholds$relative_mcse)
      data.frame(implementation = x$implementation[[1L]],
        nburn = x$nburn[[1L]], nit = x$nit[[1L]],
        failed_estimands = sum(!x$overall_pass),
        maximum_rhat = max(x$rhat), minimum_bulk_ess = min(x$ess_bulk),
        minimum_tail_ess = min(x$ess_tail),
        maximum_relative_mcse = max(x$relative_mcse),
        limiting_estimand = x$estimand[[which.max(score)]],
        limiting_diagnostic = x$limiting_diagnostic[[which.max(score)]],
        pass = all(x$overall_pass), stringsAsFactors = FALSE)
    }))
  candidates <- candidates[order(candidates$nburn + candidates$nit,
    candidates$nit, candidates$nburn), ]
  passed <- candidates[candidates$pass, , drop = FALSE]
  if (!nrow(passed)) stop("No supported Study 07 convergence setting for ",
    draws$implementation[[1L]], call. = FALSE)
  selected <- passed[1L, ]
  list(diagnostics = diagnostics, candidates = candidates,
    recommendation = data.frame(implementation = selected$implementation,
      recommendation_status = "available", nburn = selected$nburn,
      nit = selected$nit, nthin = config$convergence$nthin,
      nchains = 4L, ncores = 4L,
      limiting_estimand = selected$limiting_estimand,
      limiting_diagnostic = selected$limiting_diagnostic,
      maximum_rhat = selected$maximum_rhat,
      minimum_bulk_ess = selected$minimum_bulk_ess,
      minimum_tail_ess = selected$minimum_tail_ess,
      maximum_relative_mcse = selected$maximum_relative_mcse,
      stringsAsFactors = FALSE))
}
