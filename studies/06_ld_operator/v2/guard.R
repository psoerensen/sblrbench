.study06v2_low_rank_configuration <- function(configuration, config) {
  prop <- switch(configuration,
    low_rank_full = config$eigen_prop_full,
    low_rank_0999 = unname(config$eigen_props[["low_rank_0999"]]),
    low_rank_0995 = unname(config$eigen_props[["low_rank_0995"]]),
    NULL)
  if (is.null(prop)) return(list())
  list(representation = "low_rank", eigen_prop = prop)
}

.study06v2_assert_fit_spec <- function(configuration, controls, config,
                                       phase = NULL) {
  low_rank <- startsWith(configuration, "low_rank_")
  if (identical(configuration, "dense_reconstructed_unfiltered") &&
      !identical(phase, "operator-pilot"))
    stop("The reconstructed-dense comparator is validation-only and may run only in the operator pilot.",
      call. = FALSE)
  dense_validation <- identical(configuration,
    "dense_reconstructed_unfiltered") && identical(phase, "operator-pilot")
  forbidden_names <- intersect(names(controls),
    c("eigen_filter", "eigen_tau", "eigen_eta"))
  if (length(forbidden_names) && !dense_validation)
    stop("Study 06 v2 scientific fits prohibit legacy eigen controls: ",
      paste(forbidden_names, collapse = ", "), call. = FALSE)
  has_dense <- any(vapply(controls, function(x)
    is.character(x) && any(x == "dense_reconstructed"), logical(1)))
  if (has_dense && !dense_validation)
    stop("Study 06 v2 prohibits dense_reconstructed execution.", call. = FALSE)
  if (dense_validation && (!identical(controls$representation,
      "dense_reconstructed") || !identical(controls$eigen_policy,
      "ridge_fixed") || !identical(controls$eigen_eta, 0) ||
      !identical(controls$eigen_tau, 0)))
    stop("Validation-only dense execution must be unfiltered ridge_fixed with eta zero.",
      call. = FALSE)
  if (low_rank && (!identical(controls$representation, "low_rank") ||
      length(controls$eigen_prop) != 1L ||
      !is.finite(controls$eigen_prop) || controls$eigen_prop <= 0 ||
      controls$eigen_prop >= 1))
    stop("Every Study 06 v2 low-rank fit requires explicit representation and eigen_prop.",
      call. = FALSE)
  if (!low_rank && !dense_validation &&
      any(c("representation", "eigen_prop") %in% names(controls)))
    stop("Non-low-rank Study 06 v2 fits must not receive eigen controls.",
      call. = FALSE)
  invisible(TRUE)
}

.study06v2_validate_grid <- function(config) {
  grid <- expand.grid(architecture = config$architectures,
    replicate = seq_len(config$replicate_count),
    configuration = config$configurations, stringsAsFactors = FALSE)
  if (nrow(grid) != config$expected_fit_count || anyDuplicated(grid))
    stop("Study 06 v2 expected 60-fit grid contract failed.", call. = FALSE)
  for (configuration in config$configurations) {
    controls <- .study06v2_low_rank_configuration(configuration, config)
    .study06v2_assert_fit_spec(configuration, controls, config)
  }
  grid
}
