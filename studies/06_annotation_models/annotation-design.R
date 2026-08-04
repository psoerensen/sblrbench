# Study 06-specific annotation construction and annotation-driven simulation.

annotation_zscore <- function(x) {
  z <- (x - mean(x)) / stats::sd(x)
  if (any(!is.finite(z)))
    stop("Annotation standardization failed.", call. = FALSE)
  as.numeric(z)
}

validate_annotation_design <- function(A, marker_ids, spec) {
  required <- spec$annotation_design$columns
  if (!is.matrix(A) || !is.numeric(A))
    stop("Annotation design must be a numeric matrix.", call. = FALSE)
  if (!identical(rownames(A), marker_ids))
    stop("Annotation marker IDs or order differ from the benchmark markers.",
      call. = FALSE)
  if (!identical(colnames(A), required))
    stop("Annotation columns must preserve this exact order: ",
      paste(required, collapse = ", "), ".", call. = FALSE)
  if (any(!is.finite(A)))
    stop("Annotation design contains missing or non-finite values.",
      call. = FALSE)
  if (any(A[, "Intercept"] != 1))
    stop("The annotation intercept must equal one for every marker.",
      call. = FALSE)
  if (!all(A[, "enriched_binary"] %in% c(0, 1)))
    stop("enriched_binary must contain only zero and one.", call. = FALSE)
  expected <- max(1L, round(length(marker_ids) *
    spec$annotation_design$enriched_fraction))
  if (sum(A[, "enriched_binary"]) != expected)
    stop("The enriched binary group size differs from the audited design.",
      call. = FALSE)
  continuous <- A[, c("continuous_signal", "null_annotation"), drop = FALSE]
  if (any(abs(colMeans(continuous)) > 1e-10) ||
      any(abs(apply(continuous, 2L, stats::sd) - 1) > 1e-10))
    stop("Continuous annotations must use the audited sample z-score rule.",
      call. = FALSE)
  if (qr(A)$rank != ncol(A))
    stop("Annotation design is not full column rank.", call. = FALSE)
  invisible(TRUE)
}

construct_annotation_design <- function(marker_ids, spec) {
  if (!is.character(marker_ids) || anyNA(marker_ids) ||
      anyDuplicated(marker_ids))
    stop("marker_ids must be unique non-missing strings.", call. = FALSE)
  set.seed(as.integer(spec$annotation_design$seed))
  m <- length(marker_ids)
  selected <- sort(sample.int(m, max(1L, round(m *
    spec$annotation_design$enriched_fraction))))
  enriched <- integer(m)
  enriched[selected] <- 1L
  A <- cbind(Intercept = 1,
    enriched_binary = enriched,
    continuous_signal = annotation_zscore(stats::rnorm(m)),
    null_annotation = annotation_zscore(stats::rnorm(m)))
  rownames(A) <- marker_ids
  storage.mode(A) <- "double"
  validate_annotation_design(A, marker_ids, spec)
  A
}

annotation_marker_probabilities <- function(A, alpha, mixture_var) {
  out <- sblr::sbayesrc_marker_pi(A, alpha, gamma = mixture_var)
  if (!identical(dim(out), c(nrow(A), length(mixture_var))) ||
      any(!is.finite(out)) || any(out < 0 | out > 1) ||
      any(abs(rowSums(out) - 1) > 1e-12))
    stop("Invalid annotation-implied component probabilities.",
      call. = FALSE)
  rownames(out) <- rownames(A)
  colnames(out) <- paste0("component_", seq_len(ncol(out)) - 1L)
  out
}

annotation_reverse_sticks <- function(component_probability) {
  p <- as.numeric(component_probability)
  if (length(p) < 2L || any(!is.finite(p)) || any(p <= 0) ||
      abs(sum(p) - 1) > 1e-10)
    stop("Invalid marginal component probabilities.", call. = FALSE)
  remaining <- rev(cumsum(rev(p)))
  stats::qnorm(vapply(seq_len(length(p) - 1L),
    function(j) remaining[j + 1L] / remaining[j], numeric(1)))
}

construct_annotation_truth <- function(A, spec) {
  simulation <- spec$controls$simulation
  target <- simulation$target_expected_nonnull / nrow(A)
  nonintercept <- spec$annotation_design$informative_nonintercept_alpha
  active <- simulation$active_component_weights
  objective <- function(intercept) {
    alpha <- rbind(Intercept = c(intercept,
      stats::qnorm(sum(active[-1L])),
      stats::qnorm(active[3L] / sum(active[2:3]))), nonintercept)
    sum(annotation_marker_probabilities(A, alpha,
      simulation$mixture_var)[, -1L, drop = FALSE]) - target * nrow(A)
  }
  intercept <- stats::uniroot(objective, c(-8, -1), tol = 1e-12)$root
  informative <- rbind(Intercept = c(intercept,
    stats::qnorm(sum(active[-1L])),
    stats::qnorm(active[3L] / sum(active[2:3]))), nonintercept)
  informative_pi <- annotation_marker_probabilities(A, informative,
    simulation$mixture_var)
  marginal <- colMeans(informative_pi)
  uninformative <- matrix(0, nrow = ncol(A),
    ncol = length(simulation$mixture_var) - 1L,
    dimnames = list(colnames(A), paste0("step_", 1:3)))
  uninformative["Intercept", ] <- annotation_reverse_sticks(marginal)
  uninformative_pi <- annotation_marker_probabilities(A, uninformative,
    simulation$mixture_var)
  if (max(abs(colMeans(uninformative_pi) - marginal)) > 1e-10)
    stop("Uninformative marginal calibration failed.", call. = FALSE)
  enriched <- A[, "enriched_binary"] == 1
  expected_nonnull <- 1 - informative_pi[, 1L]
  share <- sum(expected_nonnull[enriched]) / sum(expected_nonnull)
  if (sum(expected_nonnull) < 49.5 || sum(expected_nonnull) > 50.5 ||
      share < 0.50 || share > 0.70 ||
      any(uninformative[-1L, , drop = FALSE] != 0))
    stop("Annotation truth does not satisfy the audited design.",
      call. = FALSE)
  list(informative_annotations = informative,
    uninformative_annotations = uninformative,
    marginal_component_probability = marginal,
    expected_nonnull = sum(expected_nonnull),
    enriched_expected_nonnull_share = share)
}

annotation_design_summary <- function(A, truth, spec) {
  do.call(rbind, lapply(names(spec$scenarios), function(scenario) {
    probability <- annotation_marker_probabilities(A, truth[[scenario]],
      spec$controls$simulation$mixture_var)
    enriched <- A[, "enriched_binary"] == 1
    data.frame(scenario = scenario, marker_count = nrow(A),
      annotation_count = ncol(A), enriched_count = sum(enriched),
      enriched_prevalence = mean(enriched),
      continuous_signal_mean = mean(A[, "continuous_signal"]),
      continuous_signal_sd = stats::sd(A[, "continuous_signal"]),
      null_annotation_mean = mean(A[, "null_annotation"]),
      null_annotation_sd = stats::sd(A[, "null_annotation"]),
      expected_nonnull = sum(1 - probability[, 1L]),
      enriched_expected_nonnull_share =
        sum((1 - probability[, 1L])[enriched]) / sum(1 - probability[, 1L]),
      stringsAsFactors = FALSE)
  }))
}

sample_annotation_components <- function(probability, seed) {
  set.seed(as.integer(seed))
  u <- stats::runif(nrow(probability))
  cumulative <- t(apply(probability, 1L, cumsum))
  component <- rowSums(cumulative < u) + 1L
  if (any(component < 1L | component > ncol(probability)))
    stop("Annotation component allocation failed.", call. = FALSE)
  component
}

simulate_annotation_architecture <- function(coordinate, scaled_genotypes,
                                              train_rows, A, truth, spec) {
  scenario <- as.character(coordinate$scenario)
  if (!scenario %in% names(spec$scenarios))
    stop("Unknown annotation scenario: ", scenario, call. = FALSE)
  if (!identical(colnames(scaled_genotypes), rownames(A)))
    stop("Simulation annotation-marker alignment failed.", call. = FALSE)
  mixture <- spec$controls$simulation$mixture_var
  probability <- annotation_marker_probabilities(A, truth[[scenario]], mixture)
  component <- sample_annotation_components(probability,
    coordinate$component_seed)
  set.seed(as.integer(coordinate$effect_seed))
  raw_effect <- numeric(ncol(scaled_genotypes))
  active <- component > 1L
  raw_effect[active] <- stats::rnorm(sum(active),
    sd = sqrt(mixture[component[active]]))
  sanity <- spec$controls$simulation$nonnull_sanity_range
  if (sum(active) < sanity[1L] || sum(active) > sanity[2L])
    stop("Realized active-marker count is outside the audited range.",
      call. = FALSE)
  genetic_raw <- as.numeric(scaled_genotypes %*% raw_effect)
  h2 <- spec$controls$simulation$h2
  scale <- sqrt((h2 / (1 - h2)) /
    stats::var(genetic_raw[train_rows]))
  effects <- raw_effect * scale
  genetic <- as.numeric(scaled_genotypes %*% effects)
  set.seed(as.integer(coordinate$residual_seed))
  residual <- stats::rnorm(nrow(scaled_genotypes))
  residual <- residual - mean(residual[train_rows])
  residual <- residual / stats::sd(residual[train_rows])
  phenotype <- genetic + residual
  trait <- spec$data$trait
  marker_ids <- colnames(scaled_genotypes)
  sample_ids <- rownames(scaled_genotypes)
  effect_matrix <- matrix(effects, ncol = 1L,
    dimnames = list(marker_ids, trait))
  genetic_matrix <- matrix(genetic, ncol = 1L,
    dimnames = list(sample_ids, trait))
  residual_matrix <- matrix(residual, ncol = 1L,
    dimnames = list(sample_ids, trait))
  phenotype_matrix <- matrix(phenotype, ncol = 1L,
    dimnames = list(sample_ids, trait))
  causal_ids <- marker_ids[active]
  raw <- list(y = phenotype_matrix, W = scaled_genotypes, B = effect_matrix,
    G = genetic_matrix, E = residual_matrix,
    causal = list(shared = causal_ids,
      specific = stats::setNames(list(character()), trait), all = causal_ids),
    rsids = marker_ids, ids = sample_ids, h2_target = h2,
    h2_observed = stats::var(genetic[train_rows]) /
      (stats::var(genetic[train_rows]) + stats::var(residual[train_rows])),
    shared_idx = which(active),
    specific_idx = stats::setNames(list(integer()), trait),
    causal_rsids = causal_ids)
  simulation <- as_sblrbench_simulation(raw, study = spec$study,
    architecture = scenario, replicate = as.integer(coordinate$replicate),
    seed = as.integer(coordinate$component_seed))
  simulation$data$train_ids <- sample_ids[train_rows]
  simulation$data$test_ids <- sample_ids[-train_rows]
  marker_truth <- data.frame(marker_id = marker_ids,
    component_index = component - 1L,
    component_variance = mixture[component], raw_effect = raw_effect,
    effect = effects, true_nonnull = active,
    enriched_binary = as.integer(A[, "enriched_binary"]),
    stringsAsFactors = FALSE)
  for (j in seq_len(ncol(probability)))
    marker_truth[[paste0("true_prior_component_", j - 1L)]] <- probability[, j]
  simulation$extras$annotation_matrix <- A
  simulation$extras$true_alpha <- truth[[scenario]]
  simulation$extras$true_marker_prior <- probability
  simulation$extras$marker_truth <- marker_truth
  simulation$extras$component_counts <- tabulate(component,
    nbins = length(mixture))
  simulation$extras$effect_scale <- scale
  validate_sblrbench_simulation(simulation)
  simulation
}

prepare_annotation_simulations <- function(spec, profile, data) {
  A <- construct_annotation_design(data$markers$marker_ids, spec)
  truth <- construct_annotation_truth(A, spec)
  coordinates <- unique(benchmark_annotation_seeds(spec, profile,
    mode = "final")[c("scenario", "replicate", "component_seed",
      "effect_seed", "residual_seed")])
  bundles <- lapply(seq_len(nrow(coordinates)), function(i) {
    coordinate <- as.list(coordinates[i, , drop = FALSE])
    simulation <- simulate_annotation_architecture(coordinate,
      data$scaled$all, data$split$train_rows, A, truth, spec)
    oracle <- check_oracle_genetic_values(simulation,
      tolerance = spec$validation$oracle_tolerance)
    if (!isTRUE(oracle$ok))
      stop("Study 06 simulation oracle failed.", call. = FALSE)
    test_simulation <- subset_sblrbench_simulation_samples(simulation,
      data$split$test_ids)
    stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
      spec$data)
    list(coordinate = coordinate, simulation = simulation,
      test_simulation = test_simulation, stats = stats, oracle = oracle,
      annotations = A, annotation_truth = truth,
      marker_truth = simulation$extras$marker_truth)
  })
  attr(bundles, "annotations") <- A
  attr(bundles, "annotation_truth") <- truth
  bundles
}
