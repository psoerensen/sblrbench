# Reusable numerical comparisons for reference and approximate LD operators.

#' Compare two finite square operators
#'
#' @param reference Reference matrix.
#' @param candidate Approximate matrix with the same dimensions.
#' @param tolerance Near-zero eigenvalue tolerance.
#' @return A one-row data frame of operator comparison metrics.
#' @export
operator_matrix_metrics <- function(reference, candidate, tolerance = 1e-10) {
  valid <- function(x) is.matrix(x) && nrow(x) == ncol(x) && all(is.finite(x))
  if (!valid(reference) || !valid(candidate) || !identical(dim(reference), dim(candidate)))
    stop("reference and candidate must be equally sized finite square matrices.",
      call. = FALSE)
  difference <- candidate - reference
  scale <- sqrt(sum(reference^2))
  eigenvalues <- eigen((candidate + t(candidate)) / 2,
    symmetric = TRUE, only.values = TRUE)$values
  positive <- eigenvalues[eigenvalues > tolerance]
  data.frame(
    relative_frobenius_error = if (scale == 0) NA_real_ else
      sqrt(sum(difference^2)) / scale,
    maximum_absolute_error = max(abs(difference)),
    diagonal_error = max(abs(diag(difference))),
    symmetry_error = max(abs(candidate - t(candidate))),
    trace_error = abs(sum(diag(candidate)) - sum(diag(reference))),
    negative_eigenvalue_count = sum(eigenvalues < -tolerance),
    negative_eigenvalue_mass = sum(abs(eigenvalues[eigenvalues < -tolerance])),
    near_zero_eigenvalue_count = sum(abs(eigenvalues) <= tolerance),
    effective_rank = length(positive),
    condition_number = if (length(positive)) max(positive) / min(positive) else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Compare operator products and quadratic forms
#'
#' @param reference Reference matrix.
#' @param candidate Approximate matrix.
#' @param vectors Numeric matrix whose columns are deterministic probes.
#' @return A one-row data frame.
#' @export
operator_action_metrics <- function(reference, candidate, vectors) {
  if (!is.matrix(vectors) || nrow(vectors) != nrow(reference))
    stop("vectors must have one row per operator marker.", call. = FALSE)
  reference_action <- reference %*% vectors
  candidate_action <- candidate %*% vectors
  reference_quadratic <- colSums(vectors * reference_action)
  candidate_quadratic <- colSums(vectors * candidate_action)
  data.frame(
    product_maximum_absolute_error = max(abs(candidate_action - reference_action)),
    quadratic_form_maximum_absolute_error =
      max(abs(candidate_quadratic - reference_quadratic)),
    stringsAsFactors = FALSE
  )
}
