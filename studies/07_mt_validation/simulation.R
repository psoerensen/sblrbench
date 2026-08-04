.study07_simulation_seed <- function(config, architecture, replicate) {
  i <- match(architecture, config$contract_architectures)
  if (is.na(i) || replicate < 1L) stop("Invalid simulation coordinates.")
  as.integer(config$seeds$simulation_base +
    i * config$seeds$architecture_stride +
    replicate * config$seeds$replicate_stride)
}

.study07_allocate_states <- function(marker_count, architecture, seed,
                                     config) {
  counts <- config$simulation$state_counts[[architecture]]
  if (is.null(counts) || sum(counts) > marker_count)
    stop("Invalid Study 07 causal-state specification.", call. = FALSE)
  set.seed(seed + config$seeds$state_offset)
  selected <- sample.int(marker_count, sum(counts), replace = FALSE)
  state <- integer(marker_count)
  cursor <- 0L
  for (label in names(counts)) if (counts[[label]] > 0L) {
    idx <- selected[cursor + seq_len(counts[[label]])]
    state[idx] <- match(label,
      c("trait1_only", "trait2_only", "both"))
    cursor <- cursor + counts[[label]]
  }
  state
}

.study07_orthogonal_residuals <- function(n, seed, trait_names) {
  set.seed(seed)
  e1 <- scale(stats::rnorm(n), center = TRUE, scale = TRUE)[, 1L]
  raw <- stats::rnorm(n); raw <- raw - mean(raw)
  e2 <- raw - sum(raw * e1) / sum(e1^2) * e1
  e2 <- e2 / stats::sd(e2)
  out <- cbind(e1, e2)
  colnames(out) <- trait_names
  out
}

.study07_simulate <- function(Z, architecture, replicate, config) {
  marker_ids <- colnames(Z); sample_ids <- rownames(Z)
  seed <- .study07_simulation_seed(config, architecture, replicate)
  state_id <- .study07_allocate_states(ncol(Z), architecture, seed, config)
  state <- .study07_state_inclusion(state_id, config$trait_names)
  rownames(state) <- marker_ids
  set.seed(seed + config$seeds$effect_offset)
  B <- matrix(0, ncol(Z), 2L,
    dimnames = list(marker_ids, config$trait_names))
  t1 <- which(state[, 1L] == 1L & state[, 2L] == 0L)
  t2 <- which(state[, 1L] == 0L & state[, 2L] == 1L)
  shared <- which(rowSums(state) == 2L)
  if (length(t1)) B[t1, 1L] <- stats::rnorm(length(t1))
  if (length(t2)) B[t2, 2L] <- stats::rnorm(length(t2))
  if (length(shared)) {
    common <- stats::rnorm(length(shared))
    B[shared, 1L] <- common
    B[shared, 2L] <- if (architecture == "fully_shared_negative")
      -common else if (architecture == "partially_shared")
        0.7 * common + sqrt(1 - 0.7^2) * stats::rnorm(length(shared)) else
          common
  }
  target_vg <- config$simulation$h2 / (1 - config$simulation$h2)
  G <- Z %*% B
  scales <- sqrt(target_vg / apply(G, 2L, stats::var))
  B <- sweep(B, 2L, scales, `*`)
  G <- Z %*% B
  E <- .study07_orthogonal_residuals(nrow(Z),
    seed + config$seeds$residual_offset, config$trait_names)
  Y <- G + E
  rownames(G) <- rownames(E) <- rownames(Y) <- sample_ids
  cov_g <- stats::cov(G); cov_e <- stats::cov(E)
  h2 <- diag(cov_g) / (diag(cov_g) + diag(cov_e))
  list(
    architecture = architecture, replicate = as.integer(replicate),
    seed = seed, marker_ids = marker_ids, sample_ids = sample_ids,
    state_id = state_id, state = state, effects = B,
    genetic_values = G, residuals = E, phenotype = Y,
    truth = list(cov_g = cov_g, cov_e = cov_e,
      genetic_correlation = cov2cor(cov_g)[1L, 2L],
      heritability = h2,
      state_probabilities = table(factor(state_id, levels = 0:3)) /
        length(state_id), state_counts = table(factor(state_id, levels = 0:3)),
      effect_scales = scales))
}

.study07_validate_simulation <- function(x, Z, config, tolerance = 0.015) {
  .study07_assert_matrix_alignment(x$effects, colnames(Z),
    config$trait_names, "effect")
  .study07_assert_matrix_alignment(x$phenotype, rownames(Z),
    config$trait_names, "phenotype")
  if (any(abs(x$truth$heritability - config$simulation$h2) > tolerance) ||
      max(abs(stats::cov(x$genetic_values) - x$truth$cov_g)) > 1e-12 ||
      max(abs(stats::cov(x$residuals) - x$truth$cov_e)) > 1e-12 ||
      abs(x$truth$cov_e[1L, 2L]) > 1e-10 ||
      !identical(.study07_state_id(x$state, config$trait_names), x$state_id))
    stop("Study 07 simulation/covariance contract failed.", call. = FALSE)
  reconstructed <- Z %*% x$effects
  if (!isTRUE(all.equal(unname(reconstructed),
      unname(x$genetic_values), tolerance = 1e-10)))
    stop("Study 07 genetic-value identity failed.", call. = FALSE)
  invisible(TRUE)
}

.study07_simulation_summary <- function(x) data.frame(
  architecture = x$architecture, replicate = x$replicate,
  seed = x$seed, marker_count = length(x$marker_ids),
  trait1_only_count = sum(x$state_id == 1L),
  trait2_only_count = sum(x$state_id == 2L),
  shared_count = sum(x$state_id == 3L),
  genetic_variance_trait1 = x$truth$cov_g[1L, 1L],
  genetic_variance_trait2 = x$truth$cov_g[2L, 2L],
  genetic_covariance = x$truth$cov_g[1L, 2L],
  genetic_correlation = x$truth$genetic_correlation,
  residual_covariance = x$truth$cov_e[1L, 2L],
  heritability_trait1 = x$truth$heritability[[1L]],
  heritability_trait2 = x$truth$heritability[[2L]],
  stringsAsFactors = FALSE)
