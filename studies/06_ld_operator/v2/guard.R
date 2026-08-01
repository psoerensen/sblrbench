.study06v2_low_rank_configuration <- function(configuration, config) {
  prop <- switch(configuration,
    low_rank_full = config$eigen_prop_full,
    low_rank_0999 = unname(config$eigen_props[["low_rank_0999"]]),
    low_rank_0995 = unname(config$eigen_props[["low_rank_0995"]]),
    NULL)
  if (is.null(prop)) return(list())
  list(representation = "low_rank", eigen_prop = prop)
}

.study06v2_assert_fit_spec <- function(configuration, controls, config) {
  low_rank <- startsWith(configuration, "low_rank_")
  forbidden_names <- intersect(names(controls),
    c("eigen_filter", "eigen_tau", "eigen_eta"))
  if (length(forbidden_names))
    stop("Study 06 v2 scientific fits prohibit legacy eigen controls: ",
      paste(forbidden_names, collapse = ", "), call. = FALSE)
  if (any(vapply(controls, function(x)
      is.character(x) && any(x == "dense_reconstructed"), logical(1))))
    stop("Study 06 v2 prohibits dense_reconstructed execution.", call. = FALSE)
  if (low_rank && (!identical(controls$representation, "low_rank") ||
      length(controls$eigen_prop) != 1L ||
      !is.finite(controls$eigen_prop) || controls$eigen_prop <= 0 ||
      controls$eigen_prop >= 1))
    stop("Every Study 06 v2 low-rank fit requires explicit representation and eigen_prop.",
      call. = FALSE)
  if (!low_rank && any(c("representation", "eigen_prop") %in% names(controls)))
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
