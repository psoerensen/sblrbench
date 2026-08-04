.study07_state_table <- function(trait_names = c("trait1", "trait2")) {
  stopifnot(length(trait_names) == 2L, !anyDuplicated(trait_names))
  data.frame(
    state_id = 0:3,
    state_label = c("neither", "trait1_only", "trait2_only", "both"),
    trait1 = c(0L, 1L, 0L, 1L),
    trait2 = c(0L, 0L, 1L, 1L),
    internal_key = c("0_0", "1_0", "0_1", "1_1"),
    trait1_name = trait_names[[1L]], trait2_name = trait_names[[2L]],
    stringsAsFactors = FALSE)
}

.study07_state_models <- function(trait_names = c("trait1", "trait2")) {
  x <- .study07_state_table(trait_names)
  out <- as.matrix(x[, c("trait1", "trait2")])
  storage.mode(out) <- "integer"
  rownames(out) <- x$state_label
  colnames(out) <- trait_names
  out
}

.study07_state_id <- function(inclusion,
                              trait_names = c("trait1", "trait2")) {
  x <- as.matrix(inclusion)
  if (ncol(x) != 2L || any(!x %in% 0:1))
    stop("State inclusion must be a binary marker-by-two-trait matrix.",
      call. = FALSE)
  key <- apply(x, 1L, paste, collapse = "_")
  map <- .study07_state_table(trait_names)
  id <- map$state_id[match(key, map$internal_key)]
  if (anyNA(id)) stop("Unknown MTBLR state pattern.", call. = FALSE)
  as.integer(id)
}

.study07_state_inclusion <- function(state_id,
                                     trait_names = c("trait1", "trait2")) {
  map <- .study07_state_table(trait_names)
  if (any(!state_id %in% map$state_id))
    stop("Unknown MTBLR state ID.", call. = FALSE)
  out <- as.matrix(map[match(state_id, map$state_id),
    c("trait1", "trait2"), drop = FALSE])
  storage.mode(out) <- "integer"
  colnames(out) <- trait_names
  out
}

.study07_trait_pip_from_state_probabilities <- function(probabilities,
                                                         models) {
  p <- as.matrix(probabilities); models <- as.matrix(models)
  if (ncol(p) != nrow(models) || ncol(models) != 2L ||
      any(!is.finite(p)) || any(p < 0) ||
      any(abs(rowSums(p) - 1) > 1e-10))
    stop("Invalid joint-state probability matrix.", call. = FALSE)
  out <- p %*% models
  colnames(out) <- colnames(models)
  out
}

.study07_permute_states <- function(probabilities) {
  p <- as.matrix(probabilities)
  if (ncol(p) != 4L) stop("Expected four two-trait state probabilities.")
  p[, c(1L, 3L, 2L, 4L), drop = FALSE]
}

.study07_validate_state_contract <- function(config) {
  expected <- .study07_state_models(config$trait_names)
  if (!identical(unname(expected), unname(config$models)) ||
      !identical(unname(apply(expected, 1L, paste, collapse = "_")),
        c("0_0", "1_0", "0_1", "1_1")))
    stop("Two-trait joint-state ordering differs from the verified contract.",
      call. = FALSE)
  ids <- .study07_state_id(expected, config$trait_names)
  if (!identical(ids, 0:3) ||
      !identical(unname(.study07_state_inclusion(ids,
        config$trait_names)), unname(expected)))
    stop("Joint-state round-trip contract failed.", call. = FALSE)
  invisible(TRUE)
}
