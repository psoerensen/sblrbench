.study05_zscore <- function(x) {
  z <- (x - mean(x)) / stats::sd(x)
  if (any(!is.finite(z))) stop("Annotation standardization failed.", call. = FALSE)
  as.numeric(z)
}

.study05_annotation_design <- function(marker_ids, config) {
  stopifnot(is.character(marker_ids), !anyNA(marker_ids), !anyDuplicated(marker_ids))
  set.seed(config$annotation$seed)
  m <- length(marker_ids)
  annotated <- sort(sample.int(m, max(1L, round(m * config$annotation$enriched_fraction))))
  enriched <- integer(m); enriched[annotated] <- 1L
  A <- cbind(
    Intercept = 1,
    enriched_binary = enriched,
    continuous_signal = .study05_zscore(stats::rnorm(m)),
    null_annotation = .study05_zscore(stats::rnorm(m)))
  rownames(A) <- marker_ids
  storage.mode(A) <- "double"
  .study05_validate_annotation(A, marker_ids, config)
  A
}

.study05_validate_annotation <- function(A, marker_ids, config) {
  required <- c("Intercept", "enriched_binary", "continuous_signal", "null_annotation")
  if (!is.matrix(A) || !is.numeric(A) || !identical(rownames(A), marker_ids) ||
      !identical(colnames(A), required) || any(!is.finite(A)) ||
      any(A[, "Intercept"] != 1) ||
      !all(A[, "enriched_binary"] %in% c(0, 1)) ||
      abs(mean(A[, "enriched_binary"]) - config$annotation$enriched_fraction) > 0.005 ||
      any(abs(colMeans(A[, c("continuous_signal", "null_annotation"), drop = FALSE])) > 1e-10) ||
      any(abs(apply(A[, c("continuous_signal", "null_annotation"), drop = FALSE], 2, stats::sd) - 1) > 1e-10) ||
      qr(A)$rank != ncol(A) || anyDuplicated(as.data.frame(A), MARGIN = 2L))
    stop("Study 05 annotation design contract failed.", call. = FALSE)
  invisible(TRUE)
}

.study05_marker_probabilities <- function(A, alpha, mixture_var = c(0, 0.01, 0.1, 1)) {
  out <- sblr::sbayesrc_marker_pi(A, alpha, gamma = mixture_var)
  if (!identical(dim(out), c(nrow(A), length(mixture_var))) ||
      any(!is.finite(out)) || any(out < 0 | out > 1) ||
      any(abs(rowSums(out) - 1) > 1e-12))
    stop("Invalid annotation-implied prior component probabilities.", call. = FALSE)
  out
}

.study05_reverse_sticks <- function(component_probability) {
  p <- as.numeric(component_probability)
  if (length(p) < 2L || any(!is.finite(p)) || any(p <= 0) ||
      abs(sum(p) - 1) > 1e-10) stop("Invalid marginal component probabilities.", call. = FALSE)
  remaining <- rev(cumsum(rev(p)))
  stats::qnorm(vapply(seq_len(length(p) - 1L),
    function(j) remaining[j + 1L] / remaining[j], numeric(1)))
}

.study05_true_alpha <- function(A, config) {
  target <- config$simulation$target_expected_nonnull / nrow(A)
  nonintercept <- config$simulation$informative_nonintercept_alpha
  objective <- function(intercept) {
    alpha <- rbind(Intercept = c(intercept,
      stats::qnorm(sum(config$simulation$active_component_weights[-1])),
      stats::qnorm(config$simulation$active_component_weights[3] /
        sum(config$simulation$active_component_weights[2:3]))),
      nonintercept)
    sum(.study05_marker_probabilities(A, alpha, config$mixture_var)[, -1L, drop = FALSE]) -
      target * nrow(A)
  }
  intercept1 <- stats::uniroot(objective, c(-8, -1), tol = 1e-12)$root
  informative <- rbind(Intercept = c(intercept1,
    stats::qnorm(sum(config$simulation$active_component_weights[-1])),
    stats::qnorm(config$simulation$active_component_weights[3] /
      sum(config$simulation$active_component_weights[2:3]))),
    nonintercept)
  informative_pi <- .study05_marker_probabilities(A, informative, config$mixture_var)
  marginal <- colMeans(informative_pi)
  uninformative <- matrix(0, nrow = ncol(A), ncol = length(config$mixture_var) - 1L,
    dimnames = list(colnames(A), paste0("step_", seq_len(length(config$mixture_var) - 1L))))
  uninformative["Intercept", ] <- .study05_reverse_sticks(marginal)
  uninformative_pi <- .study05_marker_probabilities(A, uninformative, config$mixture_var)
  if (max(abs(colMeans(uninformative_pi) - marginal)) > 1e-10)
    stop("Uninformative marginal calibration failed.", call. = FALSE)
  enriched <- A[, "enriched_binary"] == 1
  expected_nonnull <- rowSums(informative_pi[, -1L, drop = FALSE])
  enriched_share <- sum(expected_nonnull[enriched]) / sum(expected_nonnull)
  if (sum(expected_nonnull) < 49.5 || sum(expected_nonnull) > 50.5 ||
      enriched_share < 0.50 || enriched_share > 0.70 ||
      any(colSums(informative_pi[, -1L, drop = FALSE]) <= 0))
    stop("Frozen informative alpha does not satisfy design targets.", call. = FALSE)
  list(informative_annotations = informative,
    uninformative_annotations = uninformative,
    marginal_component_probability = marginal,
    expected_nonnull = sum(expected_nonnull),
    enriched_expected_nonnull_share = enriched_share)
}

.study05_annotation_summary <- function(A, alpha_bundle, config) {
  do.call(rbind, lapply(names(alpha_bundle)[names(alpha_bundle) %in% config$scenarios],
    function(scenario) {
      pi <- .study05_marker_probabilities(A, alpha_bundle[[scenario]], config$mixture_var)
      enriched <- A[, "enriched_binary"] == 1
      data.frame(scenario = scenario, marker_count = nrow(A),
        annotation_count = ncol(A), enriched_count = sum(enriched),
        enriched_prevalence = mean(enriched),
        continuous_signal_mean = mean(A[, "continuous_signal"]),
        continuous_signal_sd = stats::sd(A[, "continuous_signal"]),
        null_annotation_mean = mean(A[, "null_annotation"]),
        null_annotation_sd = stats::sd(A[, "null_annotation"]),
        design_rank = qr(A)$rank, expected_nonnull = sum(1 - pi[, 1]),
        enriched_expected_nonnull_share =
          sum((1 - pi[, 1])[enriched]) / sum(1 - pi[, 1]),
        stringsAsFactors = FALSE)
    }))
}
