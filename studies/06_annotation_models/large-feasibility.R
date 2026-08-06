# Study 06 large information-scale feasibility experiment.
#
# This file contains deterministic design, audit, simulation, registry, and
# analysis helpers. It does not execute a fit when sourced.

study06_large_spec <- function() {
  list(
    schema = "sblrbench-study06-large-feasibility-v1",
    profile = "large_information_scale_feasibility",
    status = "single registered replicate; not a formal qualification",
    source = list(
      qgdata_repository = "psoerensen/qgdata",
      qgdata_sha = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
      chromosome = 1L,
      sample_count = 5000L,
      glist = "results/local/06_annotation_models/checkpoints/data/human_glist.rds",
      qc = list(excludeMAF = 0.05, excludeMISS = 0.05,
        excludeCGAT = TRUE, excludeINDEL = TRUE, excludeDUPS = TRUE,
        excludeHWE = 1e-12, excludeMHC = FALSE)
    ),
    block = list(max_size = 500L, representation = "low_rank",
      eigen_policy = "cumulative_positive_mass",
      eigen_prop = 1 - .Machine$double.eps,
      positive_tolerance = 1e-10),
    seeds = list(experiment = 760000L, annotation = 760101L,
      truth = 760202L, trace_panel = 760303L, fit = 760020L,
      chain = c(760121L, 760222L, 760323L, 760424L),
      smoke = c(761121L, 761222L, 761323L, 761424L)),
    annotation = list(prevalence = 0.20, binary_continuous_correlation = 0.20,
      columns = c("intercept", "enriched_binary", "continuous_signal",
        "null_annotation")),
    mixture = list(gamma = c(0, 0.01, 0.1, 1),
      target_pi = c(null = 0.970, small = 0.015, medium = 0.010,
        large = 0.005),
      nonintercept_alpha = rbind(
        enriched_binary = c(1.00, 0.60, 0.40),
        continuous_signal = c(0.40, 0.25, 0.15),
        null_annotation = c(0, 0, 0))),
    truth = list(h2 = 0.50, active_min = 900L, active_max = 1400L,
      stick2_eligible_min = 450L, stick3_eligible_min = 150L,
      large_min = 50L),
    mcmc = list(nchains = 4L, nit = 12000L, nburn = 3000L, nthin = 1L,
      retained = 9000L, ncores = 4L, ordinary_allocation_sweeps = 1L,
      ordinary_hierarchy_updates = 1L),
    prior = list(h2 = 0.50, pi = c(0.970, 0.015, 0.010, 0.005),
      mixture_var = c(0, 0.01, 0.1, 1), sigmaSqAlpha_init = c(1, 1, 1),
      sigmaSqAlpha_a = 2, sigmaSqAlpha_b = 2, pi_floor = 1e-12,
      intercept_flat = FALSE, intercept_prior = "proper package default",
      alpha_update_every = 1L, updateB = TRUE, updateE = TRUE),
    trace = list(marker_target = 300L,
      quantities = c("b", "d", "component"), max_trace_gb = 4),
    convergence = list(rhat_max = 1.01, bulk_ess_min = 400,
      tail_ess_min = 400, relative_mcse_max = 0.05),
    recovery = list(alpha_interval_min = 5L, alpha_median_abs_max = 0.25,
      alpha_max_abs_max = 0.50, active_relative_max = 0.10,
      component_relative_max = 0.15, h2_absolute_max = 0.05)
  )
}

study06_large_hash <- function(x) digest::digest(x, algo = "sha256",
  serialize = TRUE)

study06_large_registry <- function(spec = study06_large_spec()) {
  out <- data.frame(
    fit_id = c("E0", "B0", "E2", "B2", "E1", "B1"),
    route = rep(c("bed", "block_eigen"), 3L),
    model_class = rep(c("baseline", "fixed_true_alpha", "learned_alpha"),
      each = 2L),
    method = c("bayesr", "sbayesr", "bayesrc", "sbayesrc", "bayesrc",
      "sbayesrc"),
    annotation_aware = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    update_alpha = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    execution_order = seq_len(6L), stringsAsFactors = FALSE)
  out$config_hash <- vapply(seq_len(nrow(out)), function(i)
    study06_large_hash(list(row = out[i, 1:7], mcmc = spec$mcmc,
      prior = spec$prior, block = spec$block, seeds = spec$seeds)), character(1))
  out
}

study06_large_annotations <- function(marker_ids, block_id,
                                      spec = study06_large_spec()) {
  stopifnot(length(marker_ids) == length(block_id), !anyDuplicated(marker_ids))
  set.seed(spec$seeds$annotation)
  enriched <- integer(length(marker_ids))
  for (b in unique(block_id)) {
    idx <- which(block_id == b)
    enriched[idx[order(runif(length(idx)))[seq_len(max(1L,
      round(length(idx) * spec$annotation$prevalence)))]]] <- 1L
  }
  z_binary <- as.numeric(scale(enriched))
  independent <- as.numeric(scale(rnorm(length(marker_ids))))
  rho <- spec$annotation$binary_continuous_correlation
  continuous <- as.numeric(scale(rho * z_binary + sqrt(1 - rho^2) * independent))
  null_raw <- rnorm(length(marker_ids))
  null <- as.numeric(scale(stats::residuals(stats::lm(null_raw ~ enriched +
    continuous))))
  A <- cbind(intercept = 1, enriched_binary = enriched,
    continuous_signal = continuous, null_annotation = null)
  rownames(A) <- marker_ids
  correlation <- stats::cor(A[, -1L, drop = FALSE])
  blocks_with_both <- vapply(split(enriched, block_id), function(x)
    length(unique(x)) == 2L, logical(1))
  audit <- list(prevalence = mean(enriched), correlation = correlation,
    rank = qr(A)$rank, condition_number = kappa(A), finite = all(is.finite(A)),
    blocks_with_binary_variation = sum(blocks_with_both),
    block_count = length(blocks_with_both), seed = spec$seeds$annotation,
    annotation_hash = study06_large_hash(A))
  if (!all(A[, 1L] == 1) || audit$rank != ncol(A) ||
      !audit$finite || abs(audit$prevalence - spec$annotation$prevalence) > .01 ||
      abs(correlation[1L, 2L]) < .15 || abs(correlation[1L, 2L]) > .30 ||
      max(abs(correlation[c(1, 2), 3])) > .02 || !all(blocks_with_both))
    stop("Large-feasibility annotation audit failed.", call. = FALSE)
  list(matrix = A, audit = audit)
}

study06_large_component_probabilities <- function(A, alpha) {
  q <- stats::pnorm(A %*% alpha)
  cbind(null = 1 - q[, 1L], small = q[, 1L] * (1 - q[, 2L]),
    medium = q[, 1L] * q[, 2L] * (1 - q[, 3L]),
    large = q[, 1L] * q[, 2L] * q[, 3L])
}

study06_large_calibrate_alpha <- function(A, spec = study06_large_spec()) {
  alpha <- rbind(intercept = 0, spec$mixture$nonintercept_alpha)
  target <- spec$mixture$target_pi
  solve_intercept <- function(stick, target_continuation) {
    objective <- function(value) {
      alpha[1L, stick] <- value
      q <- stats::pnorm(A %*% alpha)
      continuation <- q[, stick]
      if (stick > 1L) continuation <- continuation * apply(q[, seq_len(stick - 1L),
        drop = FALSE], 1L, prod)
      mean(continuation) - target_continuation
    }
    stats::uniroot(objective, c(-12, 12), tol = 1e-13)$root
  }
  alpha[1L, 1L] <- solve_intercept(1L, 1 - target[["null"]])
  alpha[1L, 2L] <- solve_intercept(2L, target[["medium"]] + target[["large"]])
  alpha[1L, 3L] <- solve_intercept(3L, target[["large"]])
  colnames(alpha) <- paste0("component_", 0:2, "_stick")
  probability <- study06_large_component_probabilities(A, alpha)
  expected <- colMeans(probability)
  error <- max(abs(expected - unname(target)))
  if (error > 1e-6 || any(abs(rowSums(probability) - 1) > 1e-12))
    stop("Alpha-intercept calibration failed.", call. = FALSE)
  list(alpha = alpha, probability = probability, expected_pi = expected,
    expected_count = colSums(probability), expected_active = sum(probability[, -1L]),
    maximum_absolute_error = error, alpha_hash = study06_large_hash(alpha))
}

study06_large_blocks <- function(marker_ids, positions,
                                 spec = study06_large_spec()) {
  if (is.unsorted(positions, strictly = FALSE))
    stop("Canonical marker positions are not ordered.", call. = FALSE)
  block_id <- ceiling(seq_along(marker_ids) / spec$block$max_size)
  panel <- data.frame(marker_index = seq_along(marker_ids), marker_id = marker_ids,
    chromosome = spec$source$chromosome, position_bp = positions,
    block_id = as.integer(block_id), stringsAsFactors = FALSE)
  panel$within_block <- ave(panel$marker_index, panel$block_id, FUN = seq_along)
  starts <- which(!duplicated(panel$block_id))
  if (max(table(panel$block_id)) > spec$block$max_size ||
      !identical(panel$marker_id, marker_ids))
    stop("Large-feasibility block construction failed.", call. = FALSE)
  list(panel = panel, block_start = as.integer(starts),
    block_hash = study06_large_hash(panel),
    marker_to_block_hash = study06_large_hash(panel[c("marker_id", "block_id")]))
}

study06_large_truth <- function(A, calibrated, glist, marker_ids, sample_ids,
                                block_id, spec = study06_large_spec()) {
  set.seed(spec$seeds$truth)
  u <- runif(length(marker_ids))
  cumulative <- t(apply(calibrated$probability, 1L, cumsum))
  component <- rowSums(u > cumulative[, 1:3, drop = FALSE])
  raw_effect <- rnorm(length(marker_ids)) * sqrt(spec$mixture$gamma[component + 1L])
  genetic_raw <- numeric(length(sample_ids))
  for (b in unique(block_id)) {
    idx <- which(block_id == b)
    X <- benchmark_extract_scaled_genotypes(glist, spec$source$chromosome,
      sample_ids, marker_ids[idx])
    genetic_raw <- genetic_raw + as.numeric(X %*% raw_effect[idx])
  }
  target_vg <- spec$truth$h2 / (1 - spec$truth$h2)
  effect_scale <- sqrt(target_vg / stats::var(genetic_raw))
  effect <- raw_effect * effect_scale
  genetic <- genetic_raw * effect_scale
  residual <- rnorm(length(sample_ids)); residual <- residual - mean(residual)
  residual <- residual / stats::sd(residual)
  target_ve <- stats::var(genetic) * (1 - spec$truth$h2) / spec$truth$h2
  residual <- residual * sqrt(target_ve)
  phenotype <- genetic + residual
  realized_h2 <- stats::var(genetic) / (stats::var(genetic) + stats::var(residual))
  marker_truth <- data.frame(marker_id = marker_ids, block_id = block_id,
    enriched_binary = A[, "enriched_binary"], component = component,
    active = component > 0L, raw_effect = raw_effect, effect = effect,
    true_prior_nonnull = 1 - calibrated$probability[, 1L],
    stringsAsFactors = FALSE)
  list(marker_truth = marker_truth, genetic_value = genetic, residual = residual,
    phenotype = phenotype, component_count = tabulate(component + 1L, nbins = 4L),
    effect_scale = effect_scale, genetic_variance = stats::var(genetic),
    residual_variance = stats::var(residual), realized_h2 = realized_h2)
}

study06_large_truth_gates <- function(truth, A, calibrated,
                                      spec = study06_large_spec()) {
  component <- truth$marker_truth$component
  enriched <- A[, "enriched_binary"] == 1
  eligible <- list(stick1 = rep(TRUE, length(component)),
    stick2 = component >= 1L, stick3 = component >= 2L)
  continuation <- list(stick1 = component >= 1L, stick2 = component >= 2L,
    stick3 = component >= 3L)
  subset_rows <- list(); pass <- TRUE; reasons <- character()
  for (s in seq_along(eligible)) {
    e <- eligible[[s]]; z <- continuation[[s]]
    rank <- qr(A[e, , drop = FALSE])$rank
    cells <- table(enriched[e], z[e])
    ok <- all(dim(cells) == c(2L, 2L)) && all(cells > 0L) && rank == ncol(A)
    subset_rows[[s]] <- data.frame(stick = s, eligible = sum(e),
      continuation = sum(z[e]), stopping = sum(!z[e]), design_rank = rank,
      minimum_binary_outcome_cell = if (length(cells) == 4L) min(cells) else 0L,
      pass = ok)
    if (!ok) reasons <- c(reasons, paste0("stick_", s, " subset variation/rank"))
  }
  count <- truth$component_count
  block_large <- table(truth$marker_truth$block_id[component == 3L])
  gates <- c(finite = all(is.finite(c(truth$phenotype, truth$genetic_value,
      truth$residual, truth$marker_truth$effect))),
    h2 = abs(truth$realized_h2 - spec$truth$h2) <= 1e-8,
    active = sum(component > 0L) >= spec$truth$active_min &&
      sum(component > 0L) <= spec$truth$active_max,
    stick2 = sum(component >= 1L) >= spec$truth$stick2_eligible_min,
    stick3 = sum(component >= 2L) >= spec$truth$stick3_eligible_min,
    all_stick_subsets = all(vapply(subset_rows, `[[`, logical(1), "pass")),
    large_count = count[4L] >= spec$truth$large_min,
    large_enrichment_groups = all(table(enriched[component == 3L]) > 0L),
    no_dominant_large_block = !length(block_large) || max(block_large) /
      sum(block_large) < .25)
  if (any(!gates)) reasons <- c(reasons, names(gates)[!gates])
  list(pass = all(gates), gates = gates, stick = do.call(rbind, subset_rows),
    reasons = unique(reasons), component_count = count,
    expected_count = calibrated$expected_count)
}

study06_large_trace_panel <- function(marker_truth, A, positions,
                                      spec = study06_large_spec()) {
  set.seed(spec$seeds$trace_panel)
  m <- nrow(marker_truth); target <- spec$trace$marker_target
  magnitude <- cut(rank(abs(marker_truth$effect), ties.method = "first"),
    breaks = c(0, .25, .5, .75, 1) * m, include.lowest = TRUE, labels = FALSE)
  continuous <- cut(A[, "continuous_signal"], unique(stats::quantile(
    A[, "continuous_signal"], seq(0, 1, .2))), include.lowest = TRUE,
    labels = FALSE)
  stratum <- interaction(marker_truth$component, marker_truth$enriched_binary,
    magnitude, continuous, drop = TRUE)
  selected <- unlist(lapply(split(seq_len(m), stratum), function(idx)
    idx[order(runif(length(idx)))[seq_len(min(length(idx), 2L))]]),
    use.names = FALSE)
  if (length(selected) < target) {
    remaining <- setdiff(seq_len(m), selected)
    selected <- c(selected, remaining[order(runif(length(remaining)))[
      seq_len(target - length(selected))]])
  }
  if (length(selected) > target)
    selected <- selected[order(abs(marker_truth$effect[selected]), decreasing = TRUE,
      method = "radix")[seq_len(target)]]
  selected <- sort(unique(selected))
  panel <- cbind(marker_truth[selected, ], position_bp = positions[selected],
    continuous_signal = A[selected, "continuous_signal"])
  if (nrow(panel) != target) stop("Selected-marker panel size is invalid.",
    call. = FALSE)
  list(panel = panel, hash = study06_large_hash(panel), seed = spec$seeds$trace_panel)
}

study06_large_controls <- function(fit_row, alpha_truth,
                                   spec = study06_large_spec(), smoke = FALSE) {
  controls <- list(nit = if (smoke) 12L else spec$mcmc$nit,
    nburn = if (smoke) 4L else spec$mcmc$nburn, nthin = spec$mcmc$nthin,
    seed = if (smoke) spec$seeds$smoke[[1L]] else spec$seeds$fit,
    nchains = spec$mcmc$nchains, ncores = spec$mcmc$ncores,
    chain_seeds = if (smoke) spec$seeds$smoke else spec$seeds$chain,
    keep_chains = TRUE, convergence = "extended", verbose = FALSE,
    h2 = spec$prior$h2, convergence_control = list(warn = FALSE,
      extended_groups = if (fit_row$annotation_aware)
        c("annotations", "probability") else "probability",
      keep_traces = TRUE, max_trace_gb = spec$trace$max_trace_gb,
      allow_large_traces = FALSE), updateB = TRUE, updateE = TRUE)
  if (fit_row$annotation_aware) {
    controls$alpha_init <- if (fit_row$model_class == "fixed_true_alpha")
      alpha_truth else {
        neutral <- alpha_truth * 0
        target <- spec$mixture$target_pi
        neutral[1L, ] <- c(stats::qnorm(1 - target[["null"]]),
          stats::qnorm((target[["medium"]] + target[["large"]]) /
            (1 - target[["null"]])),
          stats::qnorm(target[["large"]] /
            (target[["medium"]] + target[["large"]])))
        neutral
      }
    controls$sigmaSqAlpha_init <- spec$prior$sigmaSqAlpha_init
    controls$sigmaSqAlpha_a <- spec$prior$sigmaSqAlpha_a
    controls$sigmaSqAlpha_b <- spec$prior$sigmaSqAlpha_b
    controls$pi_floor <- spec$prior$pi_floor
    controls$alpha_update_every <- 1L
    controls$updateAlpha <- isTRUE(fit_row$update_alpha)
    controls$add_intercept <- FALSE
    controls$standardize_annotations <- FALSE
    controls$center_binary_annotations <- FALSE
  } else {
    controls$pi <- spec$prior$pi
  }
  controls
}

study06_large_ld_audit <- function(glist, sample_ids, panel,
                                   spec = study06_large_spec()) {
  rows <- vector("list", max(panel$block_id))
  sampled <- vector("list", max(panel$block_id))
  for (b in seq_along(rows)) {
    idx <- which(panel$block_id == b)
    X <- benchmark_extract_scaled_genotypes(glist, spec$source$chromosome,
      sample_ids, panel$marker_id[idx])
    R <- crossprod(X) / nrow(X)
    diagonal_min <- min(diag(R)); diagonal_max <- max(diag(R))
    symmetry_error <- max(abs(R - t(R)))
    eig <- eigen(R, symmetric = TRUE)
    positive <- eig$values > spec$block$positive_tolerance
    retained <- which(cumsum(eig$values[positive]) /
      sum(eig$values[positive]) >= spec$block$eigen_prop)[1L]
    retained <- if (is.na(retained)) sum(positive) else retained
    reconstructed <- tcrossprod(eig$vectors[, positive, drop = FALSE] *
      rep(sqrt(eig$values[positive]), each = nrow(R)))
    score <- crossprod(X, seq_len(nrow(X)) / nrow(X))
    score_reconstructed <- R %*% solve(R, score)
    rows[[b]] <- data.frame(block_id = b, marker_count = length(idx),
      start_marker = panel$marker_id[min(idx)], end_marker = panel$marker_id[max(idx)],
      physical_span_bp = diff(range(panel$position_bp[idx])),
      symmetry_error = symmetry_error, diagonal_min = diagonal_min,
      diagonal_max = diagonal_max,
      minimum_eigenvalue = min(eig$values),
      minimum_positive_eigenvalue = min(eig$values[positive]),
      maximum_eigenvalue = max(eig$values), positive_rank = sum(positive),
      retained_rank = retained,
      reconstruction_error = max(abs(R - reconstructed)),
      transformed_score_error = max(abs(score - score_reconstructed)),
      stringsAsFactors = FALSE)
    take <- unique(round(seq(1, ncol(X), length.out = min(12L, ncol(X)))))
    sampled[[b]] <- X[, take, drop = FALSE]
  }
  audit <- do.call(rbind, rows)
  sample_matrix <- do.call(cbind, sampled)
  sample_block <- rep(seq_along(sampled), vapply(sampled, ncol, integer(1)))
  sample_cor <- abs(stats::cor(sample_matrix)); diag(sample_cor) <- NA_real_
  sample_cor[outer(sample_block, sample_block, `==`)] <- NA_real_
  adjacent_max <- numeric(nrow(audit) - 1L)
  for (b in seq_along(adjacent_max)) {
    left <- which(panel$block_id == b); right <- which(panel$block_id == b + 1L)
    X1 <- benchmark_extract_scaled_genotypes(glist, spec$source$chromosome,
      sample_ids, panel$marker_id[left])
    X2 <- benchmark_extract_scaled_genotypes(glist, spec$source$chromosome,
      sample_ids, panel$marker_id[right])
    adjacent_max[b] <- max(abs(crossprod(X1, X2) / length(sample_ids)))
  }
  block_size_summary <- summary(audit$marker_count)
  block_size_summary <- stats::setNames(as.numeric(block_size_summary),
    names(block_size_summary))
  summary <- list(block_count = nrow(audit),
    block_size = as.list(block_size_summary),
    all_positive_modes_retained = all(audit$retained_rank == audit$positive_rank),
    omitted_nonpositive_modes = sum(audit$marker_count - audit$positive_rank),
    maximum_sampled_cross_block_absolute_correlation = max(sample_cor, na.rm = TRUE),
    maximum_complete_adjacent_block_absolute_correlation = max(adjacent_max),
    maximum_reconstruction_error = max(audit$reconstruction_error),
    maximum_transformed_score_error = max(audit$transformed_score_error))
  if (any(audit$symmetry_error > 1e-10) ||
      any(!is.finite(audit$diagonal_min)) || any(audit$diagonal_min <= 0) ||
      !summary$all_positive_modes_retained ||
      summary$maximum_reconstruction_error > 1e-7 ||
      summary$maximum_transformed_score_error > 1e-6)
    stop("Large-feasibility LD/block audit failed: max symmetry=",
      max(audit$symmetry_error), "; diagonal range=",
      min(audit$diagonal_min), "..", max(audit$diagonal_max),
      "; full positive retention=", summary$all_positive_modes_retained,
      "; max reconstruction=", summary$maximum_reconstruction_error,
      "; max transformed-score=", summary$maximum_transformed_score_error,
      call. = FALSE)
  list(blocks = audit, summary = summary)
}

study06_large_gwas <- function(glist, phenotype, marker_ids,
                               spec = study06_large_spec()) {
  stats <- sblr::make_summary_stats(Glist = glist,
    y = matrix(phenotype, dimnames = list(glist$ids, "trait_1")),
    chr = spec$source$chromosome, rows = seq_along(glist$ids), scale = TRUE,
    nthreads = 1L)
  if (!identical(stats$marker_names, marker_ids) ||
      !identical(as.integer(stats$n), spec$source$sample_count) ||
      any(!is.finite(unlist(stats[c("wy", "ww", "yy")]))))
    stop("Large-feasibility GWAS construction failed.", call. = FALSE)
  list(stats = stats, hash = study06_large_hash(stats), n = stats$n)
}

study06_large_fit_call <- function(fit_row, controls, bundle,
                                   spec = study06_large_spec()) {
  selected <- bundle$trace_panel$panel$marker_id
  controls$convergence_control$selected_markers <- selected
  controls$convergence_control$selected_marker_quantities <- spec$trace$quantities
  common <- controls
  if (fit_row$route == "bed") {
    common$method <- fit_row$method
    common$y <- matrix(bundle$truth$phenotype,
      dimnames = list(bundle$sample_ids, "trait_1"))
    common$Glist <- bundle$glist
    if (fit_row$annotation_aware) {
      common$annotation <- as.data.frame(bundle$annotations)
      common$annot_alpha_init <- common$alpha_init
      common$annot_sigma_sq_alpha_init <- common$sigmaSqAlpha_init
      common$annot_alpha_update_every <- common$alpha_update_every
      common$alpha_init <- NULL; common$sigmaSqAlpha_init <- NULL
      common$alpha_update_every <- NULL
    } else common$mixture_var <- spec$prior$mixture_var
    do.call(sblr::stblr_bed, common)
  } else {
    common$method <- fit_row$method
    common$stats <- bundle$gwas$stats
    common$Glist <- bundle$glist
    common$block_start <- bundle$blocks$block_start
    common$representation <- spec$block$representation
    common$eigen_policy <- spec$block$eigen_policy
    common$eigen_prop <- spec$block$eigen_prop
    if (fit_row$annotation_aware) {
      common$gamma <- spec$prior$mixture_var
      common$annotation <- bundle$annotations
    } else {
      common$mixture_var <- spec$prior$mixture_var
    }
    do.call(sblr::stblr_block_eigen, common)
  }
}

study06_large_scalar_diagnostics <- function(fit, fit_id,
                                             spec = study06_large_spec()) {
  traces <- fit$convergence_traces
  q <- traces$quantities
  values <- traces$values
  rows <- vector("list", nrow(q))
  for (j in seq_len(nrow(q))) {
    matrix_value <- values[, , j, drop = TRUE]
    if (is.vector(matrix_value)) matrix_value <- matrix(matrix_value, ncol = 1L)
    if (all(apply(matrix_value, 2L, function(x) length(unique(x)) == 1L))) {
      rows[[j]] <- data.frame(fit_id = fit_id, quantity = q$quantity_id[j] %||%
        q$group[j], rhat = NA_real_, bulk_ess = NA_real_, tail_ess = NA_real_,
        relative_mcse = NA_real_, constant = TRUE, pass = TRUE)
    } else {
      draws <- posterior::as_draws_matrix(lapply(seq_len(ncol(matrix_value)),
        function(k) matrix_value[, k]))
      rhat <- posterior::rhat(draws)
      bulk <- posterior::ess_bulk(draws); tail <- posterior::ess_tail(draws)
      mcse <- posterior::mcse_mean(draws)
      mean_value <- mean(matrix_value)
      relative <- if (abs(mean_value) > .Machine$double.eps)
        mcse / abs(mean_value) else mcse / stats::sd(as.numeric(matrix_value))
      rows[[j]] <- data.frame(fit_id = fit_id,
        quantity = q$quantity_id[j] %||% q$group[j], rhat = rhat,
        bulk_ess = bulk, tail_ess = tail, relative_mcse = relative,
        constant = FALSE, pass = rhat <= spec$convergence$rhat_max &&
          bulk >= spec$convergence$bulk_ess_min &&
          tail >= spec$convergence$tail_ess_min &&
          relative <= spec$convergence$relative_mcse_max)
    }
  }
  do.call(rbind, rows)
}

study06_large_metrics <- function(pip, effect, marker_truth, genetic_estimate,
                                  genetic_truth, fit_id) {
  causal <- marker_truth$active
  ordering <- order(-pip, seq_along(pip)); ranked <- causal[ordering]
  average_precision <- mean(cumsum(ranked)[ranked] / which(ranked))
  r <- rank(pip, ties.method = "average")
  auroc <- (sum(r[causal]) - sum(seq_len(sum(causal)))) /
    (sum(causal) * sum(!causal))
  rank_value <- rank(-pip, ties.method = "average")
  values <- c(pip_auprc = average_precision, pip_auroc = auroc,
    causal_rank_median = median(rank_value[causal]),
    causal_rank_mean = mean(rank_value[causal]),
    causal_mean_pip = mean(pip[causal]), noncausal_mean_pip = mean(pip[!causal]),
    effect_truth_correlation = cor(effect, marker_truth$effect),
    causal_effect_truth_correlation = cor(effect[causal], marker_truth$effect[causal]),
    effect_rmse = sqrt(mean((effect - marker_truth$effect)^2)),
    genetic_value_correlation = cor(genetic_estimate, genetic_truth),
    genetic_value_slope = coef(lm(genetic_truth ~ genetic_estimate))[2L],
    genetic_value_rmse = sqrt(mean((genetic_estimate - genetic_truth)^2)))
  for (k in c(25L, 50L, 100L, 200L, 500L, 1000L)) {
    selected <- ordering[seq_len(k)]
    values[paste0("precision_", k)] <- mean(causal[selected])
    values[paste0("recall_", k)] <- sum(causal[selected]) / sum(causal)
    values[paste0("causal_count_", k)] <- sum(causal[selected])
  }
  for (level in c(.05, .10)) {
    expected_fdp <- cumsum(1 - pip[ordering]) / seq_along(ordering)
    eligible <- which(expected_fdp <= level)
    selected <- if (length(eligible)) ordering[seq_len(max(eligible))] else integer()
    suffix <- if (level == .05) "05" else "10"
    values[paste0("bayes_fdr_", suffix, "_selected")] <- length(selected)
    values[paste0("bayes_fdr_", suffix, "_true")] <- sum(causal[selected])
    values[paste0("bayes_fdr_", suffix, "_false")] <- sum(!causal[selected])
  }
  data.frame(fit_id = fit_id, metric = names(values), value = as.numeric(values),
    stringsAsFactors = FALSE)
}
