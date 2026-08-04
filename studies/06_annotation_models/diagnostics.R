.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-convergence.R"), local = TRUE)

.study06_one_diagnostic <- function(z, thresholds) {
  out <- benchmark_scalar_diagnostics(z, thresholds)
  out[setdiff(names(out), "limiting_diagnostic")]
}

.study06_diagnostics <- function(draws, burnin, retained, config) {
  z <- .study06_chain_window(draws, burnin, retained)
  groups <- split(z, z$quantity)
  out <- do.call(rbind, lapply(groups, function(x) cbind(
    data.frame(scenario = x$scenario[1L], replicate = x$replicate[1L],
      method = x$method[1L], quantity = x$quantity[1L],
      parameter_name = x$parameter_name[1L],
      annotation_name = x$annotation_name[1L],
      stick_name = x$stick_name[1L], burnin = burnin,
      retained_draws = retained, stringsAsFactors = FALSE),
    .study06_one_diagnostic(x, config$convergence_pilot$thresholds))))
  rownames(out) <- NULL
  out
}

.study06_candidate_grid <- function(config) {
  x <- expand.grid(burnin = config$convergence_pilot$burnin_candidates,
    retained_draws = config$convergence_pilot$retained_candidates)
  x <- x[x$burnin + x$retained_draws <= config$convergence_pilot$nit, ]
  x[order(x$burnin + x$retained_draws, x$burnin, x$retained_draws), ]
}

.study06_select_recommendation <- function(draws, config) {
  grid <- .study06_candidate_grid(config)
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
    .study06_diagnostics(draws, grid$burnin[i], grid$retained_draws[i], config)))
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
    source = if (supported) "Study 06 maximum-history convergence pilot" else
      "Study 06 maximum-history convergence pilot; no candidate passed",
    stringsAsFactors = FALSE)
  list(diagnostics = diagnostics, candidates = candidates,
    recommendation = recommendation)
}

.study06_validation_summary <- function(diagnostics) {
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
