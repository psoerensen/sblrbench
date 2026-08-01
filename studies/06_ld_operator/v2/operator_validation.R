.study06v2_low_rank_identity_gate <- function(inspect, beta, yy,
    source_blocks = NULL, source_score = NULL,
    tolerance = list(absolute = 3e-5, relative = 3e-6,
      sse = 1e-8, update = 3e-6, source_matrix = 3e-5,
      source_score = 3e-5)) {
  Qs <- inspect$factor
  ws <- inspect$transformed_score
  sizes <- vapply(Qs, ncol, integer(1))
  if (!length(Qs) || sum(sizes) != length(beta) ||
      any(vapply(Qs, nrow, integer(1)) <= 0L) ||
      any(!is.finite(beta)) || length(yy) != 1L || !is.finite(yy))
    stop("Invalid Study 06 v2 low-rank deterministic inputs.", call. = FALSE)
  starts <- cumsum(c(1L, head(sizes, -1L)))
  rows <- vector("list", length(Qs))
  projected_score <- numeric(length(beta))
  diagonal <- numeric(length(beta))
  quadratic <- 0
  score_dot <- 0
  residual_norm <- 0
  transformed_norm <- 0
  for (b in seq_along(Qs)) {
    Q <- Qs[[b]]
    w <- as.numeric(ws[[b]][1L, ])
    idx <- starts[b] + seq_len(sizes[b]) - 1L
    bb <- beta[idx]
    A <- crossprod(Q)
    s <- as.numeric(crossprod(Q, w))
    r <- as.numeric(w - Q %*% bb)
    d <- colSums(Q * Q)
    u <- as.numeric(crossprod(Q, r)) + d * bb
    direct_u <- s - as.numeric(A %*% bb) + d * bb
    delta <- (seq_along(bb) - mean(seq_along(bb))) * 1e-4
    residual_updated <- r
    for (j in seq_along(delta)) residual_updated <-
      residual_updated - Q[, j] * delta[j]
    residual_rebuilt <- as.numeric(w - Q %*% (bb + delta))
    source_matrix_error <- NA_real_
    source_score_error <- NA_real_
    if (!is.null(source_blocks))
      source_matrix_error <- max(abs(A - source_blocks[[b]]))
    if (!is.null(source_score))
      source_score_error <- max(abs(s - source_score[idx]))
    rows[[b]] <- data.frame(block = b, block_start = starts[b],
      block_size = sizes[b], retained_rank = nrow(Q),
      score_error = max(abs(s - as.numeric(inspect$projected_score[1L, idx]))),
      diagonal_error = max(abs(d - as.numeric(inspect$diagonal[1L, idx]))),
      residual_error = max(abs(r - as.numeric(inspect$residual[1L,
        inspect$residual_offset[b] + seq_len(nrow(Q))]))),
      marker_residual_error = max(abs(as.numeric(crossprod(Q, r)) -
        as.numeric(inspect$marker_residual[1L, idx]))),
      corrected_score_error = max(abs(u - direct_u)),
      residual_update_error = max(abs(residual_updated - residual_rebuilt)),
      source_matrix_error = source_matrix_error,
      source_score_error = source_score_error,
      runtime_diagonal_minimum = min(d), stringsAsFactors = FALSE)
    projected_score[idx] <- s
    diagonal[idx] <- d
    quadratic <- quadratic + sum((Q %*% bb)^2)
    score_dot <- score_dot + sum(bb * s)
    residual_norm <- residual_norm + sum(r^2)
    transformed_norm <- transformed_norm + sum(w^2)
  }
  block_table <- do.call(rbind, rows)
  sse_residual <- yy - transformed_norm + residual_norm
  sse_quadratic <- yy - 2 * score_dot + quadratic
  summary <- data.frame(
    maximum_score_error = max(block_table$score_error),
    maximum_diagonal_error = max(block_table$diagonal_error),
    maximum_residual_error = max(block_table$residual_error),
    maximum_marker_residual_error = max(block_table$marker_residual_error),
    maximum_corrected_score_error = max(block_table$corrected_score_error),
    maximum_residual_update_error = max(block_table$residual_update_error),
    quadratic_form_error = abs(quadratic - inspect$quadratic_form),
    projected_score_dot_error = abs(score_dot - inspect$projected_score_dot),
    residual_norm_error = abs(residual_norm - inspect$residual_norm_squared),
    transformed_score_norm_error = abs(transformed_norm -
      inspect$transformed_score_norm_squared),
    projected_sse_residual = sse_residual,
    projected_sse_quadratic = sse_quadratic,
    projected_sse_identity_error = abs(sse_residual - sse_quadratic),
    maximum_source_matrix_error = if (is.null(source_blocks)) NA_real_ else
      max(block_table$source_matrix_error),
    maximum_source_score_error = if (is.null(source_score)) NA_real_ else
      max(block_table$source_score_error), stringsAsFactors = FALSE)
  finite <- all(is.finite(unlist(summary[setdiff(names(summary),
    c("maximum_source_matrix_error", "maximum_source_score_error"))])))
  core_errors <- unlist(summary[c("maximum_score_error",
    "maximum_diagonal_error", "maximum_residual_error",
    "maximum_marker_residual_error", "maximum_corrected_score_error",
    "maximum_residual_update_error", "quadratic_form_error",
    "projected_score_dot_error", "residual_norm_error",
    "transformed_score_norm_error")])
  pass <- finite && all(core_errors <= tolerance$absolute) &&
    summary$projected_sse_identity_error <= tolerance$sse &&
    summary$projected_sse_residual >= -tolerance$sse &&
    min(block_table$runtime_diagonal_minimum) > 0
  if (!is.null(source_blocks)) pass <- pass &&
    summary$maximum_source_matrix_error <= tolerance$source_matrix
  if (!is.null(source_score)) pass <- pass &&
    summary$maximum_source_score_error <= tolerance$source_score
  summary$pass <- pass
  if (!isTRUE(pass))
    stop("Study 06 v2 deterministic low-rank identity gate failed.",
      call. = FALSE)
  list(summary = summary, blocks = block_table,
    projected_score = projected_score, diagonal = diagonal)
}

.study06v2_operator_diagnostics <- function(inspect, eigen_prop,
                                             configuration) {
  d <- inspect$diagnostics
  d$configuration <- configuration
  d$representation <- "low_rank"
  d$operator_contract <- "block_low_rank_v1"
  d$eigen_prop <- eigen_prop
  d$rank_fraction <- sum(d$retained_rank) / sum(d$block_size)
  d$work_storage_proxy <- sum(d$block_size * d$retained_rank) /
    sum(d$block_size^2)
  d
}

.study06v2_validate_full_positive_rank <- function(inspect, eigen_prop) {
  d <- inspect$diagnostics
  if (!is.finite(eigen_prop) || eigen_prop >= 1 || eigen_prop <= 0 ||
      any(d$retained_rank != d$positive_rank))
    stop("Near-full eigen_prop did not retain every positive eigenvalue.",
      call. = FALSE)
  invisible(TRUE)
}
