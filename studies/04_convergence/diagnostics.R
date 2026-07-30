.study04_diagnostics <- function(draws, burnin, retained, registry, thresholds) {
  z <- .study04_window(draws, burnin, retained); core <- registry$estimand[registry$required]
  do.call(rbind, lapply(core, function(id) {
    wide <- reshape(z[c("raw_iteration", "chain", id)], idvar = "raw_iteration", timevar = "chain", direction = "wide")
    m <- as.matrix(wide[-1]); sd_all <- stats::sd(as.numeric(m))
    rhat <- posterior::rhat(m); eb <- posterior::ess_bulk(m); et <- posterior::ess_tail(m); mc <- posterior::mcse_mean(m)
    rel <- if (is.finite(sd_all) && sd_all > 0) mc / sd_all else NA_real_
    flags <- c(rhat <= thresholds$rhat, eb >= thresholds$ess_bulk,
      et >= thresholds$ess_tail, is.finite(rel) && rel <= thresholds$relative_mcse,
      ncol(m) == thresholds$chain_count, all(is.finite(m)) && sd_all > 0)
    data.frame(architecture = z$architecture[1], method = z$method[1], burnin_candidate = burnin,
      retained_draw_candidate = retained, estimand = id, rhat = rhat, ess_bulk = eb,
      ess_tail = et, mcse_mean = mc, posterior_sd = sd_all, relative_mcse = rel,
      chain_count = ncol(m), draws_per_chain = nrow(m), total_draws = length(m),
      rhat_pass = flags[1], ess_bulk_pass = flags[2], ess_tail_pass = flags[3],
      relative_mcse_pass = flags[4], chain_count_pass = flags[5], finite_pass = flags[6],
      overall_pass = all(flags), status = if (all(flags)) "pass" else if (!flags[6]) "indeterminate" else "fail",
      reason = if (all(flags)) "" else "one or more operational thresholds failed", stringsAsFactors = FALSE)
  }))
}
