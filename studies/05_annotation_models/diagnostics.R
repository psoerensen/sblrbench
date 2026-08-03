.study05_one_diagnostic <- function(z, thresholds) {
  wide <- reshape(z[c("iteration", "chain", "value")],
    idvar = "iteration", timevar = "chain", direction = "wide")
  m <- as.matrix(wide[-1L])
  sd_all <- stats::sd(as.numeric(m))
  rhat <- posterior::rhat(m)
  bulk <- posterior::ess_bulk(m)
  tail <- posterior::ess_tail(m)
  mcse <- posterior::mcse_mean(m)
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
    rhat_pass = flags[1L], ess_bulk_pass = flags[2L],
    ess_tail_pass = flags[3L], relative_mcse_pass = flags[4L],
    chain_count_pass = flags[5L], finite_nonconstant_pass = flags[6L],
    overall_pass = all(flags), stringsAsFactors = FALSE)
}

.study05_diagnostics <- function(draws, burnin, retained, config) {
  z <- .study05_chain_window(draws, burnin, retained)
  groups <- split(z, z$quantity)
  out <- do.call(rbind, lapply(groups, function(x) cbind(
    data.frame(scenario = x$scenario[1L], replicate = x$replicate[1L],
      method = x$method[1L], quantity = x$quantity[1L],
      parameter_name = x$parameter_name[1L],
      annotation_name = x$annotation_name[1L],
      stick_name = x$stick_name[1L], burnin = burnin,
      retained_draws = retained, stringsAsFactors = FALSE),
    .study05_one_diagnostic(x, config$convergence_pilot$thresholds))))
  rownames(out) <- NULL
  out
}

.study05_candidate_grid <- function(config) {
  x <- expand.grid(burnin = config$convergence_pilot$burnin_candidates,
    retained_draws = config$convergence_pilot$retained_candidates)
  x <- x[x$burnin + x$retained_draws <= config$convergence_pilot$nit, ]
  x[order(x$burnin + x$retained_draws, x$burnin, x$retained_draws), ]
}

.study05_select_recommendation <- function(draws, config) {
  grid <- .study05_candidate_grid(config)
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    .study05_diagnostics(draws, grid$burnin[i], grid$retained_draws[i], config)))
  key <- interaction(diagnostics$burnin, diagnostics$retained_draws, drop = TRUE)
  candidates <- do.call(rbind, lapply(split(diagnostics, key), function(z) {
    failed <- z[!z$overall_pass, , drop = FALSE]
    data.frame(method = z$method[1L], burnin = z$burnin[1L],
      retained_draws = z$retained_draws[1L],
      all_estimands_pass = all(z$overall_pass),
      limiting_estimand = if (nrow(failed)) failed$quantity[1L] else "none",
      limiting_diagnostic = if (nrow(failed)) {
        checks <- c(rhat = failed$rhat_pass[1L],
          ess_bulk = failed$ess_bulk_pass[1L],
          ess_tail = failed$ess_tail_pass[1L],
          relative_mcse = failed$relative_mcse_pass[1L],
          finite_nonconstant = failed$finite_nonconstant_pass[1L])
        names(checks)[which(!checks)[1L]]
      } else "none",
      maximum_rhat = max(z$rhat), minimum_bulk_ess = min(z$ess_bulk),
      minimum_tail_ess = min(z$ess_tail),
      maximum_relative_mcse = max(z$relative_mcse),
      stringsAsFactors = FALSE)
  }))
  candidates <- candidates[order(candidates$burnin + candidates$retained_draws,
    candidates$burnin, candidates$retained_draws), ]
  pass <- candidates[candidates$all_estimands_pass, , drop = FALSE]
  supported <- nrow(pass) > 0L
  selected <- if (supported) pass[1L, , drop = FALSE] else
    candidates[nrow(candidates), , drop = FALSE]
  recommendation <- data.frame(method = selected$method,
    recommendation_status = if (supported) "available" else "unsupported",
    nchains = 4L, ncores = 4L,
    nburn = selected$burnin, nit = selected$retained_draws, nthin = 1L,
    limiting_estimand = selected$limiting_estimand,
    limiting_diagnostic = selected$limiting_diagnostic,
    source = if (supported) "Study 05 maximum-history convergence pilot" else
      "Study 05 maximum-history convergence pilot; no candidate passed",
    stringsAsFactors = FALSE)
  list(diagnostics = diagnostics, candidates = candidates,
    recommendation = recommendation)
}

.study05_validation_summary <- function(diagnostics) {
  groups <- split(diagnostics, interaction(diagnostics$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    rep_pass <- tapply(z$overall_pass, z$replicate, all)
    failed <- z[!z$overall_pass, , drop = FALSE]
    data.frame(method = z$method[1L], replicate_count = length(rep_pass),
      pass_count = sum(rep_pass), pass_proportion = mean(rep_pass),
      maximum_rhat = max(z$rhat), minimum_bulk_ess = min(z$ess_bulk),
      minimum_tail_ess = min(z$ess_tail),
      maximum_relative_mcse = max(z$relative_mcse),
      limiting_estimand = if (nrow(failed)) failed$quantity[1L] else "none",
      limiting_diagnostic = if (nrow(failed)) {
        checks <- c(rhat = failed$rhat_pass[1L],
          ess_bulk = failed$ess_bulk_pass[1L],
          ess_tail = failed$ess_tail_pass[1L],
          relative_mcse = failed$relative_mcse_pass[1L])
        names(checks)[which(!checks)[1L]]
      } else "none",
      limiting_replicate = if (nrow(failed)) failed$replicate[1L] else NA_integer_,
      recommendation_validation_status =
        if (all(rep_pass)) "supported_in_all_replicates"
        else if (any(rep_pass)) "partially_supported" else "not_supported",
      stringsAsFactors = FALSE)
  }))
}
