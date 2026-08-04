.study06_simulation_seeds <- function(scenario, replicate, config) {
  s <- match(scenario, config$scenarios)
  base <- config$seeds$scenario_base + s * 100000L +
    as.integer(replicate) * config$seeds$replicate_stride
  c(component_allocation = base + config$seeds$component_offset,
    effect_generation = base + config$seeds$effect_offset,
    phenotype_residuals = base + config$seeds$residual_offset)
}

.study06_sample_components <- function(probability, seed) {
  set.seed(seed)
  u <- stats::runif(nrow(probability))
  cumulative <- t(apply(probability, 1L, cumsum))
  component <- rowSums(cumulative < u) + 1L
  if (any(component < 1L | component > ncol(probability)))
    stop("Component allocation failed.", call. = FALSE)
  component
}

.study06_simulate <- function(scenario, replicate, Z_all, train_rows, A,
                              alpha_bundle, config) {
  if (!scenario %in% config$scenarios || !replicate %in% seq_len(5L))
    stop("Invalid Study 06 simulation specification.", call. = FALSE)
  if (!identical(colnames(Z_all), rownames(A)))
    stop("Simulation annotation-marker alignment failed.", call. = FALSE)
  seeds <- .study06_simulation_seeds(scenario, replicate, config)
  alpha <- alpha_bundle[[scenario]]
  prior_probability <- .study06_marker_probabilities(A, alpha, config$mixture_var)
  component <- .study06_sample_components(prior_probability,
    seeds[["component_allocation"]])
  set.seed(seeds[["effect_generation"]])
  raw_effect <- numeric(ncol(Z_all))
  active <- component > 1L
  raw_effect[active] <- stats::rnorm(sum(active),
    sd = sqrt(config$mixture_var[component[active]]))
  if (sum(active) < config$simulation$nonnull_sanity_range[1L] ||
      sum(active) > config$simulation$nonnull_sanity_range[2L])
    stop("Realized non-null count outside prespecified sanity range: ",
      sum(active), call. = FALSE)
  genetic_raw <- as.numeric(Z_all %*% raw_effect)
  target_vg <- config$simulation$h2 / (1 - config$simulation$h2)
  raw_vg <- stats::var(genetic_raw[train_rows])
  if (!is.finite(raw_vg) || raw_vg <= 0)
    stop("Non-finite training genetic variance.", call. = FALSE)
  effect_scale <- sqrt(target_vg / raw_vg)
  effect <- raw_effect * effect_scale
  genetic_value <- as.numeric(Z_all %*% effect)
  set.seed(seeds[["phenotype_residuals"]])
  residual <- stats::rnorm(nrow(Z_all))
  residual <- residual - mean(residual[train_rows])
  residual <- residual / stats::sd(residual[train_rows])
  phenotype <- genetic_value + residual
  vg <- stats::var(genetic_value[train_rows])
  ve <- stats::var(residual[train_rows])
  vp <- stats::var(phenotype[train_rows])
  h2 <- vg / (vg + ve)
  if (any(!is.finite(c(vg, ve, vp, h2))) ||
      abs(h2 - config$simulation$h2) > 1e-10)
    stop("Simulation variance truth contract failed.", call. = FALSE)
  enriched <- A[, "enriched_binary"] == 1
  marker_truth <- data.frame(
    marker_id = colnames(Z_all),
    component_index = component - 1L,
    component_variance = config$mixture_var[component],
    scaled_component_variance = config$mixture_var[component] * effect_scale^2,
    raw_effect = raw_effect, effect = effect,
    true_nonnull = active,
    enriched_binary = as.integer(enriched),
    stringsAsFactors = FALSE)
  for (j in seq_len(ncol(prior_probability)))
    marker_truth[[paste0("true_prior_component_", j - 1L)]] <- prior_probability[, j]
  list(
    scenario = scenario, replicate = as.integer(replicate),
    seeds = seeds, alpha = alpha, prior_probability = prior_probability,
    marker_truth = marker_truth,
    effect = stats::setNames(effect, colnames(Z_all)),
    genetic_value = stats::setNames(genetic_value, rownames(Z_all)),
    residual = stats::setNames(residual, rownames(Z_all)),
    phenotype = matrix(phenotype, ncol = 1L,
      dimnames = list(rownames(Z_all), config$trait)),
    truth = list(target_h2 = config$simulation$h2, genetic_variance = vg,
      residual_variance = ve, phenotypic_variance = vp, heritability = h2,
      effect_variance = stats::var(effect),
      global_nonnull_proportion = mean(active),
      global_component_proportions = tabulate(component,
        nbins = length(config$mixture_var)) / length(component),
      nonnull_count = sum(active),
      component_counts = tabulate(component, nbins = length(config$mixture_var)),
      effect_scale = effect_scale,
      realized_enriched_causal_fraction =
        if (sum(active)) mean(enriched[active]) else NA_real_,
      true_enriched_nonnull_probability_fraction =
        sum((1 - prior_probability[, 1L])[enriched]) /
          sum(1 - prior_probability[, 1L]))
  )
}

.study06_simulation_summary <- function(x) data.frame(
  scenario = x$scenario, replicate = x$replicate,
  component_seed = x$seeds[["component_allocation"]],
  effect_seed = x$seeds[["effect_generation"]],
  residual_seed = x$seeds[["phenotype_residuals"]],
  target_h2 = x$truth$target_h2, realized_h2 = x$truth$heritability,
  genetic_variance = x$truth$genetic_variance,
  residual_variance = x$truth$residual_variance,
  phenotypic_variance = x$truth$phenotypic_variance,
  effect_variance = x$truth$effect_variance,
  effect_scale = x$truth$effect_scale,
  nonnull_count = x$truth$nonnull_count,
  component_counts = paste(x$truth$component_counts, collapse = ";"),
  realized_enriched_causal_fraction = x$truth$realized_enriched_causal_fraction,
  true_enriched_nonnull_probability_fraction =
    x$truth$true_enriched_nonnull_probability_fraction,
  stringsAsFactors = FALSE)
