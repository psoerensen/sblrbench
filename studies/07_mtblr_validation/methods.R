.study07_method_spec <- function(id) {
  map <- list(
    mt_bed_bayesc = list(interface = "mtblr_bed", method = "bayesc",
      operator = "packed_bed"),
    mt_csr_sbayesc = list(interface = "mtblr_csr", method = "sbayesc",
      operator = "full_csr"),
    mt_block_csr_sbayesc = list(interface = "mtblr_csr", method = "sbayesc",
      operator = "runtime_matched_block_csr"),
    mt_block_eigen_sbayesc = list(interface = "mtblr_block_eigen",
      method = "sbayesc", operator = "block_eigen"))
  out <- map[[id]]
  if (is.null(out)) stop("Unknown Study 07 implementation.", call. = FALSE)
  c(list(id = id), out)
}

.study07_fit_seeds <- function(config, architecture, replicate,
                               implementation, nchains = 4L) {
  ai <- match(architecture, config$contract_architectures)
  mi <- match(implementation, unique(c(config$runtime_implementations,
    config$implementations)))
  if (anyNA(c(ai, mi))) stop("Invalid Study 07 seed coordinates.")
  fit <- as.integer(config$seeds$fit_base +
    ai * config$seeds$architecture_stride +
    replicate * config$seeds$replicate_stride +
    mi * config$seeds$implementation_stride)
  chains <- as.integer(fit + seq_len(nchains) * config$seeds$chain_stride)
  list(fit_seed = fit, chain_seeds = chains)
}

.study07_extended_controls <- function(keep_traces = TRUE) list(
  warn = FALSE, keep_traces = keep_traces,
  extended_groups = c("covariance", "probability"),
  max_trace_gb = 1, allow_large_traces = FALSE)

.study07_pause_message <- function() paste(
  "Study 07 MT block-eigen execution is paused.",
  "The current mtblr_block_eigen() backend is reconstructed dense.",
  "Resume this phase only after the retained low-rank MT operator",
  "has been implemented and validated in sblr."
)

.study07_assert_phase_allowed <- function(phase) {
  if (phase %in% c("runtime", "convergence", "benchmark", "stress",
      "aggregate", "all"))
    stop(.study07_pause_message(), call. = FALSE)
  invisible(TRUE)
}

.study07_assert_execution_allowed <- function(implementation, phase = "fit") {
  if (identical(implementation, "mt_block_eigen_sbayesc"))
    stop(.study07_pause_message(), call. = FALSE)
  invisible(TRUE)
}

.study07_fit <- function(implementation, simulation, stats, Glist,
                         runtime_glist, rows, blocks, config, controls) {
  .study07_assert_execution_allowed(implementation, "fit")
  spec <- .study07_method_spec(implementation)
  seeds <- .study07_fit_seeds(config, simulation$architecture,
    simulation$replicate, implementation, controls$nchains)
  common <- c(controls, list(
    method = spec$method, h2 = unname(config$simulation$h2),
    models = .study07_state_models(config$trait_names),
    pimodels = unname(config$model_probabilities),
    residual_covariance = if (spec$interface == "mtblr_bed")
      "diagonal" else NULL,
    seed = seeds$fit_seed, chain_seeds = seeds$chain_seeds,
    keep_chains = TRUE, convergence = "extended",
    convergence_control = .study07_extended_controls(TRUE),
    verbose = FALSE))
  common <- common[!vapply(common, is.null, logical(1))]
  input <- switch(implementation,
    mt_bed_bayesc = list(y = simulation$phenotype[rows, , drop = FALSE],
      Glist = Glist, rows = as.integer(rows),
      cls = stats$cls),
    mt_csr_sbayesc = list(stats = stats, Glist = Glist),
    mt_block_csr_sbayesc = list(stats = stats, Glist = runtime_glist),
    mt_block_eigen_sbayesc = list(stats = stats, Glist = Glist,
      block_start = blocks$start, operator_sharing = "shared",
      eigen_filter = "ridge_fixed", eigen_eta = 0))
  started <- Sys.time(); elapsed <- proc.time()[["elapsed"]]
  warnings <- character()
  result <- tryCatch(withCallingHandlers(
    do.call(getExportedValue("sblr", spec$interface), c(input, common)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = identity)
  runtime <- proc.time()[["elapsed"]] - elapsed
  if (inherits(result, "error")) return(list(status = "failed",
    error = conditionMessage(result), fit = NULL, implementation = spec,
    runtime = runtime, warnings = warnings, seeds = seeds,
    controls = controls, started_at = started, finished_at = Sys.time()))
  validation <- try(.study07_validate_fit_contract(result,
    length(simulation$marker_ids), config$trait_names, controls$nchains),
    silent = TRUE)
  if (inherits(validation, "try-error")) return(list(status = "failed",
    error = as.character(validation), fit = result, implementation = spec,
    runtime = runtime, warnings = warnings, seeds = seeds,
    controls = controls, started_at = started, finished_at = Sys.time()))
  list(status = "ok", error = "", fit = result,
    implementation = spec, runtime = runtime, warnings = warnings,
    seeds = seeds, controls = controls, started_at = started,
    finished_at = Sys.time())
}

.study07_validate_fit_contract <- function(fit, marker_count, trait_names,
                                            chain_count) {
  reported_chain_count <- if (!is.null(fit$nchains)) fit$nchains else
    fit$input$nchains
  chain_ids <- if (is.list(fit$chains)) vapply(fit$chains, function(x) {
    value <- if (!is.null(x$chain_index)) x$chain_index else x$chain
    if (length(value) != 1L || !is.finite(value)) NA_integer_ else
      as.integer(value)
  }, integer(1)) else integer()
  if (!inherits(fit, "mtblr_fit") ||
      !identical(dim(fit$bm), c(marker_count, 2L)) ||
      !identical(colnames(fit$bm), trait_names) ||
      !identical(dim(fit$dm), c(marker_count, 2L)) ||
      any(!is.finite(fit$bm)) || any(!is.finite(fit$dm)) ||
      any(fit$dm < 0 | fit$dm > 1) ||
      !identical(as.integer(reported_chain_count), as.integer(chain_count)) ||
      length(fit$chains) != chain_count ||
      anyNA(chain_ids) || !identical(unname(sort(chain_ids)),
        seq_len(chain_count)) ||
      is.null(fit$convergence_traces$values) ||
      dim(fit$convergence_traces$values)[2L] != chain_count)
    stop("MTBLR fit dimensions, labels, chains, or trace contract failed.",
      call. = FALSE)
  for (field in c("cov_b_mean", "cov_g_mean", "cov_e_mean")) {
    x <- fit[[field]]
    if (!identical(dim(x), c(2L, 2L)) || any(!is.finite(x)) ||
        max(abs(x - t(x))) > 1e-10)
      stop("MTBLR covariance output contract failed: ", field,
        call. = FALSE)
  }
  if (is.null(fit$pi_mean) || length(fit$pi_mean) != 4L ||
      any(!is.finite(fit$pi_mean)) ||
      abs(sum(fit$pi_mean) - 1) > 1e-10)
    stop("MTBLR global joint-state probability contract failed.",
      call. = FALSE)
  invisible(TRUE)
}
