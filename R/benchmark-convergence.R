# Shared true-chain convergence mechanics. These functions never reconstruct
# traces from posterior means, final states, or pooled draws.

#' Build the matched Study 04 coordinate grid
#'
#' @param spec A validated convergence specification.
#' @param profile Either `"workshop"` or `"benchmark"`.
#' @return A deterministic data frame with explicit selection/validation stage.
#' @export
benchmark_convergence_coordinates <- function(spec, profile = "benchmark") {
  validate_benchmark_spec(spec)
  resolved <- resolve_benchmark_profile(spec, profile)
  matched <- benchmark_matched_spec(spec)
  grid <- spec$matched_grid
  if (any(!grid$scenario %in% names(matched$scenarios)) ||
      any(!grid$method %in% names(matched$methods)))
    stop("Study 04 matched grid is not present in the authoritative Study 03 specification.",
      call. = FALSE)
  rows <- list()
  if ("selection" %in% resolved$stages) {
    selection <- grid
    selection$stage <- "selection"
    selection$replicate <- 1L
    rows[[length(rows) + 1L]] <- selection
  }
  if ("validation" %in% resolved$stages) {
    validation <- do.call(rbind, lapply(seq_len(resolved$replicate_count),
      function(replicate) transform(grid, stage = "validation",
        replicate = as.integer(replicate))))
    rows[[length(rows) + 1L]] <- validation
  }
  out <- do.call(rbind, rows)
  out <- out[c("stage", "scenario", "replicate", "method")]
  out$profile <- profile
  rownames(out) <- NULL
  out
}

#' Add the historical Study 04 seeds to convergence coordinates
#'
#' @inheritParams benchmark_convergence_coordinates
#' @return The coordinate grid with simulation, fit, and four chain seeds.
#' @export
benchmark_convergence_seeds <- function(spec, profile = "benchmark") {
  coordinates <- benchmark_convergence_coordinates(spec, profile)
  matched <- benchmark_matched_spec(spec)
  architecture_index <- match(coordinates$scenario, names(matched$scenarios))
  method_index <- match(coordinates$method, names(matched$methods))
  coordinates$architecture_seed <- as.integer(
    matched$seeds$simulation_base + architecture_index * 1000L)
  coordinates$simulation_seed <- as.integer(
    coordinates$architecture_seed + coordinates$replicate)
  seeds <- lapply(seq_len(nrow(coordinates)), function(i)
    as.integer(spec$seeds$fit_base + architecture_index[[i]] * 1000000L +
      coordinates$replicate[[i]] * 10000L + method_index[[i]] * 1000L +
      seq_len(spec$diagnostics$thresholds$chain_count) *
        spec$seeds$chain_stride))
  coordinates$fit_seed <- vapply(seeds, `[[`, integer(1), 1L)
  coordinates$chain_seeds <- I(seeds)
  coordinates
}

#' Tabulate the Study 04 diagnostic design
#'
#' @param spec A validated convergence specification.
#' @return A list of ordinary data frames for quantities, candidates,
#'   thresholds, and recommendation rules.
#' @export
benchmark_convergence_design <- function(spec) {
  validate_benchmark_spec(spec)
  thresholds <- spec$diagnostics$thresholds
  list(
    quantities = spec$diagnostics$registry,
    candidates = list(
      burnin = data.frame(burnin_candidate =
        spec$diagnostics$burnin_candidates),
      retained = data.frame(retained_draw_candidate =
        spec$diagnostics$retained_draw_candidates)),
    thresholds = data.frame(
      diagnostic = c("rhat", "ess_bulk", "ess_tail", "relative_mcse",
        "chain_count", "standardized_mean_shift"),
      operator = c("<=", ">=", ">=", "<=", "==", "<="),
      threshold = unname(unlist(thresholds)), stringsAsFactors = FALSE),
    recommendation_rules = data.frame(
      rule = names(spec$diagnostics$recommendation_rules),
      definition = vapply(spec$diagnostics$recommendation_rules,
        function(x) paste(x, collapse = "; "), character(1)),
      stringsAsFactors = FALSE))
}

benchmark_chain_window <- function(draws, burnin, retained,
                                   iteration_col = "iteration",
                                   chain_col = "chain",
                                   group_cols = character(),
                                   expected_chains = 4L) {
  needed <- c(iteration_col, chain_col, group_cols)
  if (!all(needed %in% names(draws)))
    stop("Chain draws are missing required columns.", call. = FALSE)
  if (anyDuplicated(draws[c(chain_col, iteration_col, group_cols)]))
    stop("Chain draws contain duplicate chain/iteration identities.",
      call. = FALSE)
  chains <- sort(unique(draws[[chain_col]]))
  if (!identical(as.integer(chains), seq_len(as.integer(expected_chains))))
    stop("Chain identity does not match the required identifiable chains.",
      call. = FALSE)
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
  wide <- wide[order(wide[[iteration_col]]), , drop = FALSE]
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

benchmark_convergence_diagnostics <- function(draws, burnin, retained,
                                               registry, thresholds) {
  required <- registry$quantity[registry$required]
  window <- benchmark_chain_window(draws, burnin, retained,
    group_cols = "quantity", expected_chains = thresholds$chain_count)
  rows <- lapply(required, function(quantity) {
    value <- window[window$quantity == quantity, , drop = FALSE]
    result <- benchmark_scalar_diagnostics(value, thresholds)
    metadata <- value[1L, intersect(c("stage", "scenario", "replicate",
      "method"), names(value)), drop = FALSE]
    cbind(metadata, burnin_candidate = as.integer(burnin),
      retained_draw_candidate = as.integer(retained),
      quantity = quantity, result,
      status = if (result$overall_pass) "pass" else if (
        !result$finite_nonconstant_pass) "indeterminate" else "fail",
      reason = if (result$overall_pass) "" else
        "one or more operational thresholds failed",
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

benchmark_convergence_candidate_grid <- function(draws, spec) {
  maximum <- max(draws$iteration)
  rows <- list()
  for (burnin in spec$diagnostics$burnin_candidates)
    for (retained in spec$diagnostics$retained_draw_candidates)
      if (burnin + retained <= maximum)
        rows[[length(rows) + 1L]] <- benchmark_convergence_diagnostics(
          draws, burnin, retained, spec$diagnostics$registry,
          spec$diagnostics$thresholds)
  do.call(rbind, rows)
}

benchmark_burnin_stability <- function(draws, spec) {
  rules <- spec$diagnostics$recommendation_rules
  thresholds <- spec$diagnostics$thresholds
  required <- spec$diagnostics$registry$quantity[
    spec$diagnostics$registry$required]
  reference <- benchmark_chain_window(draws, rules$reference_burnin,
    rules$reference_retained_draws, group_cols = "quantity",
    expected_chains = thresholds$chain_count)
  rows <- list()
  for (burnin in spec$diagnostics$burnin_candidates) {
    candidate <- benchmark_chain_window(draws, burnin,
      rules$reference_retained_draws, group_cols = "quantity",
      expected_chains = thresholds$chain_count)
    for (quantity in required) {
      ref <- reference$value[reference$quantity == quantity]
      value <- candidate$value[candidate$quantity == quantity]
      reference_sd <- stats::sd(ref)
      shift <- if (is.finite(reference_sd) && reference_sd > 0)
        abs(mean(value) - mean(ref)) / reference_sd else NA_real_
      metadata <- draws[1L, intersect(c("scenario", "method"), names(draws)),
        drop = FALSE]
      rows[[length(rows) + 1L]] <- cbind(metadata,
        burnin_candidate = burnin, quantity = quantity,
        candidate_mean = mean(value), reference_mean = mean(ref),
        reference_posterior_sd = reference_sd,
        standardized_mean_shift = shift,
        stable = is.finite(shift) &&
          shift <= thresholds$standardized_mean_shift,
        status = if (is.finite(shift)) "ok" else "indeterminate",
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

benchmark_convergence_recommendations <- function(diagnostics, stability,
                                                   spec) {
  grid <- unique(diagnostics[c("scenario", "method")])
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    d <- diagnostics[diagnostics$scenario == grid$scenario[[i]] &
      diagnostics$method == grid$method[[i]], , drop = FALSE]
    s <- stability[stability$scenario == grid$scenario[[i]] &
      stability$method == grid$method[[i]], , drop = FALSE]
    chosen <- NULL
    for (burnin in spec$diagnostics$burnin_candidates) {
      stable <- s$stable[s$burnin_candidate == burnin]
      if (length(stable) && all(stable))
        for (retained in spec$diagnostics$retained_draw_candidates) {
          later <- d$overall_pass[d$burnin_candidate == burnin &
            d$retained_draw_candidate >= retained]
          if (length(later) && all(later)) {
            chosen <- c(burnin, retained)
            break
          }
        }
      if (!is.null(chosen)) break
    }
    final_burnin <- if (is.null(chosen))
      spec$diagnostics$recommendation_rules$reference_burnin else chosen[[1L]]
    final_retained <- if (is.null(chosen))
      spec$diagnostics$recommendation_rules$reference_retained_draws else
        chosen[[2L]]
    final <- d[d$burnin_candidate == final_burnin &
      d$retained_draw_candidate == final_retained, , drop = FALSE]
    limiting <- if (nrow(final)) final$quantity[which.max(ifelse(
      is.finite(final$rhat), final$rhat, Inf))][[1L]] else NA_character_
    data.frame(method = grid$method[[i]], scenario = grid$scenario[[i]],
      architecture = grid$scenario[[i]],
      recommendation_status = if (is.null(chosen)) "unavailable" else
        "available", recommended_nchains = 4L, recommended_ncores = 4L,
      recommended_nburn = if (is.null(chosen)) NA_integer_ else chosen[[1L]],
      recommended_post_burnin_draws = if (is.null(chosen)) NA_integer_ else
        chosen[[2L]],
      recommended_nit_argument = if (is.null(chosen)) NA_integer_ else
        chosen[[2L]], recommended_nthin = 1L,
      limiting_estimand = limiting,
      limiting_diagnostic = "maximum_rhat",
      maximum_rhat = if (nrow(final)) max(final$rhat, na.rm = TRUE) else
        NA_real_,
      minimum_ess_bulk = if (nrow(final)) min(final$ess_bulk,
        na.rm = TRUE) else NA_real_,
      minimum_ess_tail = if (nrow(final)) min(final$ess_tail,
        na.rm = TRUE) else NA_real_,
      maximum_relative_mcse = if (nrow(final)) max(final$relative_mcse,
        na.rm = TRUE) else NA_real_,
      reason = if (is.null(chosen))
        spec$diagnostics$recommendation_rules$unavailable else
        "earliest stable operational pass", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

benchmark_convergence_validation_status <- function(successful, passing,
                                                    required = 5L) {
  if (successful < required) return("indeterminate")
  if (passing == required) return("supported_in_all_replicates")
  if (passing >= ceiling(required / 2)) return("supported_in_most_replicates")
  if (passing >= 1L) return("mixed")
  "not_supported"
}

benchmark_convergence_validation_summary <- function(status, diagnostics,
                                                      required = 5L) {
  keys <- interaction(status$scenario, status$method, drop = TRUE)
  do.call(rbind, lapply(split(status, keys), function(rows) {
    d <- diagnostics[diagnostics$scenario == rows$scenario[[1L]] &
      diagnostics$method == rows$method[[1L]], , drop = FALSE]
    successful <- sum(rows$status == "ok")
    passing <- sum(rows$all_core_quantities_pass & rows$status == "ok")
    limiting <- d$quantity[!d$overall_pass]
    limiting <- if (!length(limiting)) "none" else names(sort(table(limiting),
      decreasing = TRUE))[[1L]]
    data.frame(scenario = rows$scenario[[1L]], method = rows$method[[1L]],
      replicate_count = required, successful_replicates = successful,
      replicates_passing_all_core_quantities = passing,
      pass_proportion = passing / required,
      maximum_rhat = if (nrow(d)) max(d$rhat, na.rm = TRUE) else NA_real_,
      minimum_bulk_ess = if (nrow(d)) min(d$ess_bulk, na.rm = TRUE) else
        NA_real_,
      minimum_tail_ess = if (nrow(d)) min(d$ess_tail, na.rm = TRUE) else
        NA_real_,
      maximum_relative_mcse = if (nrow(d)) max(d$relative_mcse,
        na.rm = TRUE) else NA_real_,
      most_frequent_limiting_quantity = limiting,
      recommendation_validation_status =
        benchmark_convergence_validation_status(successful, passing,
          required), stringsAsFactors = FALSE)
  }))
}
