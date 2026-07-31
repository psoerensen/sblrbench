.study06_simulation_specs <- function(config) {
  out <- list()
  for (architecture in config$architectures)
    for (replicate in seq_len(config$replicate_count))
      out[[length(out) + 1L]] <- list(
        architecture = architecture, replicate = replicate,
        simulation_seed = .study06_seed(config, architecture, replicate))
  out
}

.study06_simulate <- function(spec, Z_all, config) {
  architecture <- config$architecture_specs[[spec$architecture]]
  if (is.null(architecture))
    stop("Unknown Study 06 architecture.", call. = FALSE)
  set.seed(spec$simulation_seed)
  causal_index <- sort(sample.int(ncol(Z_all),
    config$simulation$n_causal, replace = FALSE))
  if (architecture$effect_distribution == "single_normal") {
    component <- rep("single_normal", length(causal_index))
    raw_effect <- stats::rnorm(length(causal_index))
  } else if (architecture$effect_distribution == "variance_mixture") {
    k <- sample.int(length(architecture$mixture_var),
      length(causal_index), replace = TRUE,
      prob = architecture$mixture_prob)
    component <- paste0("variance_", architecture$mixture_var[k])
    raw_effect <- stats::rnorm(length(k),
      sd = sqrt(architecture$mixture_var[k]))
  } else stop("Unsupported Study 06 effect distribution.", call. = FALSE)
  effects <- matrix(0, ncol(Z_all), 1L,
    dimnames = list(colnames(Z_all), config$trait))
  effects[causal_index, 1L] <- raw_effect
  genetic <- Z_all %*% effects
  target_vg <- config$simulation$h2 / (1 - config$simulation$h2)
  effect_scale <- sqrt(target_vg / stats::var(genetic[, 1L]))
  effects[, 1L] <- effects[, 1L] * effect_scale
  genetic <- Z_all %*% effects
  set.seed(spec$simulation_seed + config$seeds$residual_offset)
  residual <- stats::rnorm(nrow(Z_all))
  residual <- (residual - mean(residual)) / stats::sd(residual)
  residual <- matrix(residual, ncol = 1L,
    dimnames = list(rownames(Z_all), config$trait))
  phenotype <- genetic + residual
  vg <- stats::var(genetic[, 1L])
  ve <- stats::var(residual[, 1L])
  h2 <- vg / (vg + ve)
  list(
    architecture = spec$architecture,
    replicate = as.integer(spec$replicate),
    simulation_seed = as.integer(spec$simulation_seed),
    marker_ids = colnames(Z_all), sample_ids = rownames(Z_all),
    causal_index = causal_index,
    causal_marker_ids = colnames(Z_all)[causal_index],
    effects = effects, genetic_values = genetic,
    residuals = residual, phenotype = phenotype,
    component = component, raw_effect = raw_effect,
    effect_scale = effect_scale,
    truth = c(causal_proportion = length(causal_index) / ncol(Z_all),
      effect_variance = mean(effects[causal_index, 1L]^2),
      genetic_variance = vg, residual_variance = ve,
      heritability = h2))
}

.study06_validate_simulation <- function(x, Z_all, config) {
  if (!identical(x$marker_ids, colnames(Z_all)) ||
      !identical(x$sample_ids, rownames(Z_all)) ||
      length(x$causal_index) != config$simulation$n_causal ||
      anyDuplicated(x$causal_index) ||
      any(!is.finite(x$effects)) ||
      any(!is.finite(x$phenotype)) ||
      abs(unname(x$truth["heritability"]) -
        config$simulation$h2) > 0.015)
    stop("Study 06 simulation truth validation failed.", call. = FALSE)
  reconstructed <- as.numeric(Z_all %*% x$effects)
  if (!isTRUE(all.equal(reconstructed,
      as.numeric(x$genetic_values), tolerance = 1e-10)))
    stop("Study 06 genetic-value oracle failed.", call. = FALSE)
  invisible(TRUE)
}

.study06_simulation_summary <- function(x) data.frame(
  architecture = x$architecture, replicate = x$replicate,
  simulation_seed = x$simulation_seed,
  marker_count = length(x$marker_ids),
  causal_count = length(x$causal_index),
  effect_scale = x$effect_scale,
  genetic_variance = unname(x$truth["genetic_variance"]),
  residual_variance = unname(x$truth["residual_variance"]),
  phenotypic_variance = stats::var(x$phenotype[, 1L]),
  realized_heritability = unname(x$truth["heritability"]),
  stringsAsFactors = FALSE)
