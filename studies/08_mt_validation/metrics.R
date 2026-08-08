.study08_safe_cor <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2L ||
      any(!is.finite(c(x, y))) || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(NA_real_)
  stats::cor(x, y)
}

.study08_auc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || all(truth) || any(!is.finite(score))) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - sum(seq_len(sum(truth)))) /
    (sum(truth) * sum(!truth))
}

.study08_auprc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || any(!is.finite(score))) return(NA_real_)
  o <- order(score, decreasing = TRUE); y <- truth[o]
  precision <- cumsum(y) / seq_along(y)
  mean(precision[y])
}

.study08_prediction_metrics <- function(run, simulation, Z_test, test_rows) {
  B <- run$fit$bm; pred <- Z_test %*% B
  y <- simulation$phenotype[test_rows, , drop = FALSE]
  g <- simulation$genetic_values[test_rows, , drop = FALSE]
  do.call(rbind, lapply(seq_len(2L), function(t) {
    slope <- if (stats::var(pred[, t]) > 0)
      unname(stats::coef(stats::lm(y[, t] ~ pred[, t]))[[2L]]) else NA_real_
    data.frame(architecture = simulation$architecture,
      replicate = simulation$replicate,
      implementation = run$implementation$id,
      trait = colnames(B)[[t]],
      metric = c("phenotype_prediction_correlation",
        "phenotype_prediction_rmse", "prediction_slope",
        "genetic_value_correlation", "genetic_value_rmse"),
      value = c(.study08_safe_cor(pred[, t], y[, t]),
        sqrt(mean((pred[, t] - y[, t])^2)), slope,
        .study08_safe_cor(pred[, t], g[, t]),
        sqrt(mean((pred[, t] - g[, t])^2))),
      stringsAsFactors = FALSE)
  }))
}

.study08_marker_metrics <- function(run, simulation, Z_all) {
  est <- run$fit$bm; truth <- simulation$effects; pip <- run$fit$dm
  do.call(rbind, lapply(seq_len(2L), function(t) {
    causal <- simulation$state[, t] == 1L
    values <- c(marker_effect_correlation = .study08_safe_cor(est[, t], truth[, t]),
      marker_effect_rmse = sqrt(mean((est[, t] - truth[, t])^2)),
      causal_marker_effect_correlation = .study08_safe_cor(est[causal, t],
        truth[causal, t]), noncausal_marker_rmse = sqrt(mean(est[!causal, t]^2)),
      pip_auroc = .study08_auc(pip[, t], causal),
      pip_auprc = .study08_auprc(pip[, t], causal),
      pip_brier = mean((pip[, t] - causal)^2),
      genetic_value_reconstruction_correlation = .study08_safe_cor(
        Z_all %*% est[, t], simulation$genetic_values[, t]))
    data.frame(architecture = simulation$architecture,
      replicate = simulation$replicate,
      implementation = run$implementation$id,
      trait = colnames(est)[[t]], metric = names(values),
      value = unname(values), stringsAsFactors = FALSE)
  }))
}

.study08_internal_consistency <- function(run, simulation) {
  fit <- run$fit
  rg <- fit$cov_g_mean[1L, 2L] /
    sqrt(fit$cov_g_mean[1L, 1L] * fit$cov_g_mean[2L, 2L])
  checks <- c(probability_sum = abs(sum(fit$pi_mean) - 1) <= 1e-10,
    covariance_symmetry = max(abs(fit$cov_g_mean - t(fit$cov_g_mean))) <= 1e-10,
    residual_covariance_diagonal = abs(fit$cov_e_mean[1L, 2L]) <= 1e-10,
    covariance_psd = min(eigen(fit$cov_g_mean, symmetric = TRUE,
      only.values = TRUE)$values) >= -1e-8,
    genetic_correlation_identity = is.finite(rg) && abs(rg) <= 1 + 1e-8,
    effect_dimensions = identical(dim(fit$bm), dim(simulation$effects)),
    pip_dimensions = identical(dim(fit$dm), dim(simulation$state)),
    trait_labels = identical(colnames(fit$bm), colnames(simulation$effects)),
    marker_labels = identical(rownames(fit$bm), rownames(simulation$effects)))
  data.frame(architecture = simulation$architecture,
    replicate = simulation$replicate,
    implementation = run$implementation$id,
    identity = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE)
}

.study08_comparison_registry <- function() data.frame(
  comparison_id = c("full_csr_minus_bed", "block_eigen_minus_full_csr",
    "block_eigen_minus_runtime_block_csr"),
  focal = c("mt_csr_sbayesc", "mt_block_eigen_sbayesc",
    "mt_block_eigen_sbayesc"),
  reference = c("mt_bed_bayesc", "mt_csr_sbayesc",
    "mt_block_csr_sbayesc"),
  interpretation = c("likelihood_and_representation",
    "block_approximation", "runtime_operator_equivalence"),
  stringsAsFactors = FALSE)

.study08_paired_differences <- function(table) {
  registry <- .study08_comparison_registry(); out <- list()
  for (i in seq_len(nrow(registry))) for (architecture in unique(table$architecture)) {
    a <- table[table$architecture == architecture &
      table$implementation == registry$focal[[i]], ]
    b <- table[table$architecture == architecture &
      table$implementation == registry$reference[[i]], ]
    keys <- intersect(names(a), names(b))
    keys <- intersect(c("replicate", "trait", "metric"), keys)
    if (!nrow(a) || !nrow(b)) next
    z <- merge(a[c(keys, "value")], b[c(keys, "value")], by = keys,
      suffixes = c("_focal", "_reference"), all = TRUE)
    z$comparison_id <- registry$comparison_id[[i]]
    z$architecture <- architecture
    z$difference <- z$value_focal - z$value_reference
    z$complete_pair <- is.finite(z$difference)
    z$interpretation <- registry$interpretation[[i]]
    out[[length(out) + 1L]] <- z
  }
  do.call(rbind, out)
}

.study08_paired_summary <- function(x) {
  keys <- interaction(x$comparison_id, x$architecture,
    if ("trait" %in% names(x)) x$trait else "joint", x$metric, drop = TRUE)
  do.call(rbind, lapply(split(x, keys), function(z) {
    y <- z$difference[z$complete_pair]
    data.frame(comparison_id = z$comparison_id[[1L]],
      architecture = z$architecture[[1L]],
      trait = if ("trait" %in% names(z)) z$trait[[1L]] else "joint",
      metric = z$metric[[1L]], complete_paired_replicates = length(y),
      mean_difference = mean(y), sd_difference = if (length(y) > 1L)
        stats::sd(y) else NA_real_, median_difference = stats::median(y),
      minimum_difference = min(y), maximum_difference = max(y),
      interpretation = z$interpretation[[1L]], stringsAsFactors = FALSE)
  }))
}
