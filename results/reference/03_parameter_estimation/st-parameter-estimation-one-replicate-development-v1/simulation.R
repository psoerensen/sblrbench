.study03_replicate_specs <- function(config) {
  profile <- config$profiles[[config$profile]]
  if (is.null(profile)) stop("Unknown Study 03 profile.", call. = FALSE)
  n <- profile$replicate_count
  if (config$profile == "development" && n != 1L) stop("Development benchmark requires exactly one replicate.", call. = FALSE)
  out <- list()
  for (a in names(config$simulation$architectures)) for (i in seq_len(n)) {
    out[[length(out) + 1L]] <- list(architecture = a, replicate = i,
      simulation_seed = as.integer(config$simulation$base_seed +
        match(a, names(config$simulation$architectures)) * config$seeds$architecture_offset + i))
  }
  out
}

.study03_simulate <- function(spec, Z, config) {
  architecture <- config$simulation$architectures[[spec$architecture]]
  set.seed(spec$simulation_seed)
  causal_index <- sort(sample.int(ncol(Z), config$simulation$n_causal))
  if (architecture$effect_distribution == "single_normal") {
    component <- rep("single_normal", length(causal_index))
    raw_effect <- stats::rnorm(length(causal_index))
  } else {
    k <- sample.int(length(architecture$mixture_var), length(causal_index), replace = TRUE,
      prob = architecture$mixture_prob)
    component <- paste0("variance_", architecture$mixture_var[k])
    raw_effect <- stats::rnorm(length(k), sd = sqrt(architecture$mixture_var[k]))
  }
  B <- matrix(0, ncol(Z), 1L, dimnames = list(colnames(Z), config$trait))
  B[causal_index, 1L] <- raw_effect
  G <- Z %*% B
  target_vg <- config$simulation$h2 / (1 - config$simulation$h2)
  scale <- sqrt(target_vg / stats::var(G[, 1L]))
  B[, 1L] <- B[, 1L] * scale
  G <- Z %*% B
  E <- matrix(stats::rnorm(nrow(Z)), ncol = 1L, dimnames = list(rownames(Z), config$trait))
  E[, 1L] <- (E[, 1L] - mean(E[, 1L])) / stats::sd(E[, 1L])
  Y <- G + E
  h2 <- stats::var(G[, 1L]) / (stats::var(G[, 1L]) + stats::var(E[, 1L]))
  causal_ids <- colnames(Z)[causal_index]
  raw <- list(y = Y, W = Z, B = B, G = G, E = E,
    causal = list(shared = causal_ids,
      specific = stats::setNames(list(character()), config$trait), all = causal_ids),
    rsids = colnames(Z), ids = rownames(Z), h2_target = config$simulation$h2,
    h2_observed = h2, shared_idx = causal_index,
    specific_idx = stats::setNames(list(integer()), config$trait), causal_rsids = causal_ids)
  sim <- sblrbench::as_sblrbench_simulation(raw, study = config$study,
    architecture = spec$architecture, replicate = spec$replicate, seed = spec$simulation_seed)
  sim$extras$effect_components <- data.frame(marker = causal_ids, component = component,
    raw_effect = raw_effect, final_effect = B[causal_index, 1L], stringsAsFactors = FALSE)
  sim$extras$effect_scale <- scale
  sim$extras$effect_distribution <- architecture$effect_distribution
  sblrbench::validate_sblrbench_simulation(sim)
  sim
}

.study03_truth <- function(simulation, config) {
  B <- simulation$truth$effects[, 1L]
  causal <- B != 0
  vg <- stats::var(simulation$truth$genetic_values[, 1L])
  ve <- stats::var(simulation$truth$residuals[, 1L])
  values <- c(sum(causal) / length(B), mean(B[causal]^2), sum(B^2), vg, ve, vg / (vg + ve))
  ids <- .study03_estimand_registry()$estimand_id
  data.frame(architecture = simulation$scenario$architecture,
    replicate = simulation$scenario$replicate, estimand_id = ids, truth = values,
    truth_type = "realized_quantity",
    truth_definition = .study03_estimand_registry()$truth_source,
    source = "validated simulation on the fitted analysis scale", status = "ok", reason = "",
    stringsAsFactors = FALSE)
}
