# Internal scalar-chain mechanics shared across convergence studies.

benchmark_chain_window <- function(draws, burnin, retained,
                                   iteration_col = "iteration",
                                   chain_col = "chain",
                                   group_cols = character(),
                                   expected_chains = 4L) {
  needed <- c(iteration_col, chain_col, group_cols)
  if (!all(needed %in% names(draws)))
    stop("Chain draws are missing required columns.", call. = FALSE)
  keep <- draws[[iteration_col]] > burnin &
    draws[[iteration_col]] <= burnin + retained
  out <- draws[keep, , drop = FALSE]
  if (!nrow(out)) stop("Chain window is empty.", call. = FALSE)
  groups <- if (length(group_cols)) interaction(out[group_cols],
    drop = TRUE, lex.order = TRUE) else factor(rep("all", nrow(out)))
  counts <- table(groups, out[[chain_col]])
  if (ncol(counts) != expected_chains || any(counts != retained))
    stop("Chain window lacks equal retained draws.", call. = FALSE)
  out
}

benchmark_trace_array_long <- function(values, quantity_ids,
                                       metadata = list(),
                                       expected_chains = 4L) {
  if (length(dim(values)) != 3L ||
      dim(values)[2L] != expected_chains ||
      length(quantity_ids) != dim(values)[3L])
    stop("Convergence trace array dimensions are invalid.", call. = FALSE)
  base <- expand.grid(iteration = seq_len(dim(values)[1L]),
    chain = seq_len(dim(values)[2L]))
  rows <- lapply(seq_along(quantity_ids), function(i) {
    out <- data.frame(base, quantity = quantity_ids[[i]],
      value = as.vector(values[, , i]), stringsAsFactors = FALSE)
    for (name in names(metadata)) out[[name]] <- metadata[[name]]
    out
  })
  do.call(rbind, rows)
}

benchmark_scalar_diagnostics <- function(draws, thresholds,
                                         iteration_col = "iteration",
                                         chain_col = "chain",
                                         value_col = "value") {
  required_thresholds <- c("rhat", "ess_bulk", "ess_tail",
    "relative_mcse", "chain_count")
  if (!all(required_thresholds %in% names(thresholds)))
    stop("Convergence thresholds are incomplete.", call. = FALSE)
  wide <- reshape(draws[c(iteration_col, chain_col, value_col)],
    idvar = iteration_col, timevar = chain_col, direction = "wide")
  values <- as.matrix(wide[-1L])
  posterior_sd <- stats::sd(as.numeric(values))
  rhat <- posterior::rhat(values)
  bulk <- posterior::ess_bulk(values)
  tail <- posterior::ess_tail(values)
  mcse <- posterior::mcse_mean(values)
  relative <- if (is.finite(posterior_sd) && posterior_sd > 0)
    mcse / posterior_sd else NA_real_
  flags <- c(is.finite(rhat) && rhat <= thresholds$rhat,
    is.finite(bulk) && bulk >= thresholds$ess_bulk,
    is.finite(tail) && tail >= thresholds$ess_tail,
    is.finite(relative) && relative <= thresholds$relative_mcse,
    ncol(values) == thresholds$chain_count,
    all(is.finite(values)) && is.finite(posterior_sd) && posterior_sd > 0)
  names(flags) <- c("rhat", "ess_bulk", "ess_tail", "relative_mcse",
    "chain_count", "finite_nonconstant")
  data.frame(rhat = rhat, ess_bulk = bulk, ess_tail = tail,
    mcse_mean = mcse, posterior_sd = posterior_sd,
    relative_mcse = relative, chain_count = ncol(values),
    draws_per_chain = nrow(values), rhat_pass = flags[[1L]],
    ess_bulk_pass = flags[[2L]], ess_tail_pass = flags[[3L]],
    relative_mcse_pass = flags[[4L]], chain_count_pass = flags[[5L]],
    finite_nonconstant_pass = flags[[6L]], overall_pass = all(flags),
    limiting_diagnostic = if (all(flags)) "none" else
      names(flags)[which(!flags)[[1L]]], stringsAsFactors = FALSE)
}
