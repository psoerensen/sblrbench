.study04_burnin_stability <- function(draws, diagnostics, config) {
  ref <- .study04_window(draws, 1000L, 2000L); core <- .study04_registry()$estimand[.study04_registry()$required]
  do.call(rbind, lapply(config$burnin_candidates, function(b) do.call(rbind, lapply(core, function(id) {
    cand <- .study04_window(draws, b, 2000L); rsd <- stats::sd(ref[[id]])
    shift <- if (is.finite(rsd) && rsd > 0) abs(mean(cand[[id]]) - mean(ref[[id]])) / rsd else NA_real_
    data.frame(architecture = draws$architecture[1], method = draws$method[1], burnin_candidate = b,
      estimand = id, candidate_mean = mean(cand[[id]]), reference_mean = mean(ref[[id]]),
      reference_posterior_sd = rsd, standardized_mean_shift = shift,
      stable = is.finite(shift) && shift <= config$thresholds$standardized_mean_shift,
      status = if (is.finite(shift)) "ok" else "indeterminate")
  }))))
}

.study04_recommend <- function(diagnostics, stability, config) {
  methods <- unique(diagnostics[c("architecture", "method")]); out <- list()
  for (i in seq_len(nrow(methods))) {
    d <- diagnostics[diagnostics$architecture == methods$architecture[i] & diagnostics$method == methods$method[i], ]
    s <- stability[stability$architecture == methods$architecture[i] & stability$method == methods$method[i], ]
    chosen <- NULL
    for (b in config$burnin_candidates) if (all(s$stable[s$burnin_candidate == b])) {
      for (n in config$retained_draw_candidates) {
        later <- d$overall_pass[d$burnin_candidate == b & d$retained_draw_candidate >= n]
        if (length(later) && all(later)) { chosen <- c(b, n); break }
      }; if (!is.null(chosen)) break
    }
    final_burnin <- if (is.null(chosen)) 1000L else chosen[1]
    final_retained <- if (is.null(chosen)) 2000L else chosen[2]
    final <- d[d$burnin_candidate == final_burnin &
      d$retained_draw_candidate == final_retained, ]
    limiting <- if (nrow(final)) {
      final$estimand[which.max(ifelse(is.finite(final$rhat), final$rhat, Inf))][1]
    } else NA_character_
    out[[i]] <- data.frame(method = methods$method[i], architecture = methods$architecture[i],
      recommendation_status = if (is.null(chosen)) "unavailable" else "available",
      recommended_nchains = 4L, recommended_ncores = 4L,
      recommended_nburn = if (is.null(chosen)) NA_integer_ else chosen[1],
      recommended_post_burnin_draws = if (is.null(chosen)) NA_integer_ else chosen[2],
      # In the verified sblr contract, nit is the number of post-burn-in draws
      # and nburn is an additional warm-up count.
      recommended_nit_argument = if (is.null(chosen)) NA_integer_ else chosen[2], recommended_nthin = 1L,
      limiting_estimand = limiting, limiting_diagnostic = "maximum_rhat",
      maximum_rhat = if (nrow(final)) max(final$rhat, na.rm = TRUE) else NA_real_,
      minimum_ess_bulk = if (nrow(final)) min(final$ess_bulk, na.rm = TRUE) else NA_real_,
      minimum_ess_tail = if (nrow(final)) min(final$ess_tail, na.rm = TRUE) else NA_real_,
      maximum_relative_mcse = if (nrow(final)) max(final$relative_mcse, na.rm = TRUE) else NA_real_,
      reason = if (is.null(chosen)) "no stable checkpoint passed all core thresholds" else "earliest stable operational pass")
  }
  do.call(rbind, out)
}
