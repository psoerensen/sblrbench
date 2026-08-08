.study08_assert_matrix_alignment <- function(x, row_ids, column_ids,
                                               name = "matrix") {
  if (!is.matrix(x) || !identical(rownames(x), row_ids) ||
      !identical(colnames(x), column_ids) || any(!is.finite(x)))
    stop(name, " alignment or finite-value contract failed.", call. = FALSE)
  invisible(TRUE)
}

.study08_marker_subset <- function(marker_ids, marker_count) {
  if (marker_count > length(marker_ids) || marker_count < 1L)
    stop("Invalid nested marker count.", call. = FALSE)
  marker_ids[seq_len(marker_count)]
}

.study08_permutation_contract <- function(Z, Y, B, seed) {
  set.seed(seed)
  marker_order <- sample.int(ncol(Z))
  sample_order <- sample.int(nrow(Z))
  trait_order <- rev(seq_len(ncol(Y)))
  base_g <- Z %*% B
  marker_g <- Z[, marker_order, drop = FALSE] %*%
    B[marker_order, , drop = FALSE]
  sample_g <- Z[sample_order, , drop = FALSE] %*% B
  trait_g <- Z %*% B[, trait_order, drop = FALSE]
  data.frame(
    contract = c("marker_permutation", "sample_permutation",
      "trait_permutation", "sample_phenotype_permutation"),
    maximum_absolute_error = c(max(abs(base_g - marker_g)),
      max(abs(base_g[sample_order, ] - sample_g)),
      max(abs(base_g[, trait_order] - trait_g)),
      max(abs(Y[sample_order, ] - Y[sample_order, , drop = FALSE]))),
    passed = TRUE, stringsAsFactors = FALSE)
}

.study08_validate_stats <- function(stats, marker_ids, trait_names,
                                    sample_count) {
  if (!identical(stats$marker_names, marker_ids) ||
      !identical(stats$trait_names, trait_names) ||
      !(length(stats$n) %in% c(1L, 2L) &&
        all(as.integer(stats$n) == as.integer(sample_count))) ||
      any(!is.finite(unlist(stats[c("wy", "ww", "yy")],
        recursive = TRUE, use.names = FALSE))))
    stop("Study 08 summary-statistic alignment failed.", call. = FALSE)
  invisible(TRUE)
}
