.study06_operator_metrics <- function(reference, candidate, blocks,
                                      config, operator_id,
                                      reference_id = "runtime_matched_block_csr",
                                      action_vectors = list()) {
  ref <- .study06_dense_blocks(reference)
  alt <- .study06_dense_blocks(candidate)
  if (length(ref) != length(alt) || length(ref) != nrow(blocks))
    stop("Operator block counts differ.", call. = FALSE)
  set.seed(config$seeds$operator_probe)
  probes <- replicate(config$operator_probe_count,
    stats::rnorm(reference$marker_count), simplify = FALSE)
  rows <- lapply(seq_along(ref), function(k) {
    A <- ref[[k]]; B <- alt[[k]]
    delta <- B - A
    denom <- pmax(abs(A), .Machine$double.eps^0.5)
    eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    eig_alt <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
    positive <- sum(pmax(eig, 0))
    data.frame(
      operator_id = operator_id, reference_operator_id = reference_id,
      block_id = blocks$block_id[k], block_size = nrow(A),
      source_matrix_definition =
        "same installed builder packed runtime cross-product",
      filtering_policy = candidate$filter$mode,
      filtering_parameter = if (candidate$filter$mode == "hard_truncate")
        candidate$filter$tau else candidate$filter$eta,
      original_rank = qr(A)$rank,
      retained_rank = if (!is.null(candidate$diagnostics$n_kept))
        candidate$diagnostics$n_kept[k] else qr(B)$rank,
      retained_rank_proportion = if (!is.null(candidate$diagnostics$n_kept))
        candidate$diagnostics$n_kept[k] / nrow(B) else qr(B)$rank / nrow(B),
      minimum_eigenvalue = min(eig_alt),
      maximum_eigenvalue = max(eig_alt),
      negative_eigenvalue_count = sum(eig_alt < -1e-10),
      positive_eigenvalue_mass = positive,
      retained_positive_eigenvalue_mass = sum(pmax(eig_alt, 0)),
      retained_mass_proportion = if (positive > 0)
        sum(pmax(eig_alt, 0)) / positive else NA_real_,
      trace = sum(diag(B)), trace_difference = sum(diag(B)) - sum(diag(A)),
      frobenius_norm = sqrt(sum(B^2)),
      frobenius_relative_error = sqrt(sum(delta^2)) /
        max(sqrt(sum(A^2)), .Machine$double.eps),
      off_diagonal_frobenius_norm = sqrt(sum((B - diag(diag(B)))^2)),
      off_diagonal_norm_retained = sqrt(sum((B - diag(diag(B)))^2)) /
        max(sqrt(sum((A - diag(diag(A)))^2)), .Machine$double.eps),
      condition_number = {
        pos <- eig_alt[eig_alt > max(eig_alt) * .Machine$double.eps^0.5]
        if (length(pos)) max(pos) / min(pos) else Inf
      },
      runtime_diagonal_minimum = min(diag(B)),
      runtime_diagonal_maximum = max(diag(B)),
      runtime_diagonal_mean = mean(diag(B)),
      reconstruction_maximum_absolute_error = max(abs(delta)),
      reconstruction_relative_error = max(abs(delta) / denom),
      matrix_vector_maximum_error = max(vapply(probes, function(x) {
        idx <- blocks$start[k]:blocks$end[k]
        max(abs(as.numeric(B %*% x[idx]) - as.numeric(A %*% x[idx])))
      }, numeric(1))),
      quadratic_form_maximum_error = max(vapply(probes, function(x) {
        idx <- blocks$start[k]:blocks$end[k]
        abs(drop(crossprod(x[idx], B %*% x[idx])) -
          drop(crossprod(x[idx], A %*% x[idx])))
      }, numeric(1))),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (length(action_vectors)) {
    for (name in names(action_vectors)) {
      x <- as.numeric(action_vectors[[name]])
      if (length(x) != reference$marker_count || any(!is.finite(x)))
        stop("Invalid operator action vector: ", name, call. = FALSE)
      yr <- .study06_apply_blocks(ref, x)
      ya <- .study06_apply_blocks(alt, x)
      rel <- sqrt(sum((ya - yr)^2)) / max(sqrt(sum(yr^2)),
        .Machine$double.eps)
      cosine <- sum(yr * ya) / max(sqrt(sum(yr^2) * sum(ya^2)),
        .Machine$double.eps)
      out[[paste0(name, "_action_relative_error")]] <- rel
      out[[paste0(name, "_action_cosine_similarity")]] <- cosine
    }
  }
  out
}

.study06_equivalence_gate <- function(reference, candidate, blocks, config,
                                      action_vectors = list()) {
  metrics <- .study06_operator_metrics(reference, candidate, blocks, config,
    operator_id = "block_eigen_unfiltered",
    action_vectors = action_vectors)
  tol <- config$operator_tolerance
  mapping_ok <- identical(reference$marker_count, candidate$marker_count) &&
    identical(reference$block_start, candidate$block_start) &&
    identical(reference$block_size, candidate$block_size) &&
    identical(reference$block_of, candidate$block_of) &&
    identical(reference$local_of, candidate$local_of)
  diag_error <- max(abs(reference$diagonal - candidate$diagonal))
  pass <- mapping_ok && all(is.finite(candidate$diagonal)) &&
    all(candidate$diagonal > 0) &&
    diag_error <= tol$absolute &&
    max(metrics$reconstruction_maximum_absolute_error) <= tol$absolute &&
    max(metrics$reconstruction_relative_error) <= tol$relative &&
    max(metrics$matrix_vector_maximum_error) <= tol$product_absolute &&
    max(metrics$quadratic_form_maximum_error) <= tol$quadratic_absolute
  summary <- data.frame(
    gate = "runtime_matched_block_csr_vs_unfiltered_block_eigen",
    pass = pass, marker_mapping_pass = mapping_ok,
    block_coverage_pass = blocks$start[1L] == 1L &&
      blocks$end[nrow(blocks)] == reference$marker_count &&
      identical(blocks$start[-1L],
        blocks$end[-nrow(blocks)] + 1L),
    diagonal_maximum_absolute_error = diag_error,
    reconstruction_maximum_absolute_error =
      max(metrics$reconstruction_maximum_absolute_error),
    reconstruction_maximum_relative_error =
      max(metrics$reconstruction_relative_error),
    matrix_vector_maximum_error =
      max(metrics$matrix_vector_maximum_error),
    quadratic_form_maximum_error =
      max(metrics$quadratic_form_maximum_error),
    matrix_error_block = metrics$block_id[
      which.max(metrics$reconstruction_maximum_absolute_error)],
    relative_error_block = metrics$block_id[
      which.max(metrics$reconstruction_relative_error)],
    matrix_vector_error_block = metrics$block_id[
      which.max(metrics$matrix_vector_maximum_error)],
    quadratic_form_error_block = metrics$block_id[
      which.max(metrics$quadratic_form_maximum_error)],
    float32_reconstruction_contribution =
      max(metrics$reconstruction_maximum_absolute_error),
    absolute_tolerance = tol$absolute,
    relative_tolerance = tol$relative,
    product_absolute_tolerance = tol$product_absolute,
    quadratic_absolute_tolerance = tol$quadratic_absolute,
    rationale = paste("float32 packed runtime storage with",
      "double-precision validation products"),
    stringsAsFactors = FALSE)
  if (!pass) {
    first <- metrics[which.max(metrics$reconstruction_maximum_absolute_error), ]
    stop("Study 06 operator-equivalence gate failed at ",
      first$block_id, "; maximum absolute error = ",
      signif(first$reconstruction_maximum_absolute_error, 6),
      call. = FALSE)
  }
  list(summary = summary, block_metrics = metrics)
}

.study06_filter_dense_matrix <- function(A,
    mode = c("ridge_fixed", "hard_truncate", "ridge_lw"),
    tau = 0, eta = 0, lw_shrinkage = NULL) {
  mode <- match.arg(mode)
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || any(!is.finite(A)) ||
      max(abs(A - t(A))) > 1e-10 || any(diag(A) <= 0))
    stop("Synthetic Study 06 operator must be finite, symmetric, and positive-diagonal.",
      call. = FALSE)
  d <- diag(A)
  if (mode == "hard_truncate") {
    C <- A / sqrt(outer(d, d))
    e <- eigen((C + t(C)) / 2, symmetric = TRUE)
    effective <- max(tau, 0.01)
    keep <- e$values >= effective
    if (!any(keep)) keep[which.max(e$values)] <- TRUE
    V <- e$vectors[, keep, drop = FALSE]
    B <- (sqrt(d) * V) %*% diag(e$values[keep], nrow = sum(keep)) %*%
      t(sqrt(d) * V)
    return(list(matrix = B, eigenvalues = e$values,
      retained = e$values[keep], retained_rank = sum(keep),
      effective_threshold = effective,
      shrinkage_weight = 1 - sum(keep) / length(keep)))
  }
  a <- if (mode == "ridge_fixed") eta / (1 + eta) else lw_shrinkage
  if (length(a) != 1L || !is.finite(a) || a < 0 || a > 1)
    stop("Synthetic ridge shrinkage must lie in [0, 1].", call. = FALSE)
  B <- (1 - a) * A
  diag(B) <- d
  list(matrix = B, retained_rank = qr(B)$rank,
    effective_threshold = NA_real_, shrinkage_weight = a)
}

.study06_known_spectrum_matrix <- function(values) {
  p <- length(values)
  raw <- outer(seq_len(p), seq_len(p), function(i, j)
    sin(i * (j + 0.5)) + cos((i + 0.25) * j))
  Q <- qr.Q(qr(raw))
  Q %*% diag(values, nrow = p) %*% t(Q)
}

.study06_filter_symmetric_spectrum <- function(A, tau) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  effective <- max(tau, 0.01)
  keep <- e$values >= effective
  if (!any(keep)) keep[which.max(e$values)] <- TRUE
  V <- e$vectors[, keep, drop = FALSE]
  list(matrix = V %*% diag(e$values[keep], nrow = sum(keep)) %*% t(V),
    eigenvalues = e$values, retained = e$values[keep],
    retained_rank = sum(keep), effective_threshold = effective)
}

.study06_synthetic_filter_validation <- function(tolerance = 1e-10) {
  spectra <- list(
    well_conditioned = c(5, 3, 2, 1, 0.5),
    controlled_small = c(5, 2, 1, 0.05, 0.005),
    near_singular = c(5, 2, 1, 0.05, 1e-8))
  rows <- list()
  for (id in names(spectra)) {
    A <- .study06_known_spectrum_matrix(spectra[[id]])
    zero <- .study06_filter_dense_matrix(A, "ridge_fixed", eta = 0)
    for (tau in c(0.001, 0.01, 0.10)) {
      hard <- .study06_filter_symmetric_spectrum(A, tau)
      e <- eigen(A, symmetric = TRUE)
      keep <- e$values >= max(tau, 0.01)
      V <- e$vectors[, keep, drop = FALSE]
      direct <- V %*% diag(e$values[keep], nrow = sum(keep)) %*% t(V)
      x <- seq_len(nrow(A)) / nrow(A)
      rows[[length(rows) + 1L]] <- data.frame(
        matrix_id = id, requested_threshold = tau,
        effective_threshold = max(tau, 0.01),
        original_rank = nrow(A), retained_rank = hard$retained_rank,
        retained_mass_proportion = sum(hard$retained) /
          sum(pmax(hard$eigenvalues, 0)),
        zero_ridge_maximum_error = max(abs(zero$matrix - A)),
        filtered_reconstruction_maximum_error = max(abs(hard$matrix - direct)),
        matrix_vector_maximum_error = max(abs(hard$matrix %*% x - direct %*% x)),
        quadratic_form_error = abs(drop(crossprod(x, hard$matrix %*% x)) -
          drop(crossprod(x, direct %*% x))),
        minimum_diagonal = min(diag(hard$matrix)),
        pass = max(abs(zero$matrix - A)) <= tolerance &&
          max(abs(hard$matrix - direct)) <= tolerance &&
          all(diag(hard$matrix) > 0), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (!all(out$pass)) stop("Synthetic Study 06 filter validation failed.",
    call. = FALSE)
  out
}
