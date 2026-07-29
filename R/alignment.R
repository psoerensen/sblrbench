#' Strictly align markers
#' @param x A named vector, matrix, or data frame.
#' @param marker_ids Canonical marker order.
#' @param margin Matrix/data-frame margin containing marker IDs.
#' @param allow_extra Whether identifiers outside the requested set are allowed.
#' @export
align_markers <- function(x, marker_ids, margin = 1L, allow_extra = FALSE) .align_axis(x, marker_ids, margin, allow_extra, "marker")

#' Strictly align samples
#' @inheritParams align_markers
#' @param sample_ids Canonical sample order.
#' @export
align_samples <- function(x, sample_ids, margin = 1L, allow_extra = FALSE) .align_axis(x, sample_ids, margin, allow_extra, "sample")

#' Strictly align traits
#' @inheritParams align_markers
#' @param trait_names Canonical trait order.
#' @param margin Matrix/data-frame margin containing trait names (default columns).
#' @export
align_traits <- function(x, trait_names, margin = 2L, allow_extra = FALSE) .align_axis(x, trait_names, margin, allow_extra, "trait")

.align_axis <- function(x, requested, margin, allow_extra, what) {
  requested <- .assert_ids(requested, paste("requested", what, "IDs"))
  if (!is.logical(allow_extra) || length(allow_extra) != 1L || is.na(allow_extra)) stop("allow_extra must be TRUE or FALSE.", call.=FALSE)
  if (is.atomic(x) && is.null(dim(x))) {
    ids <- names(x); if (is.null(ids)) stop("Unnamed object cannot be identity-aligned by ", what, ".", call.=FALSE)
  } else if (is.matrix(x) || is.data.frame(x)) {
    if (!margin %in% c(1L,2L)) stop("margin must be 1 or 2.", call.=FALSE)
    ids <- dimnames(x)[[margin]]; if (is.null(ids)) stop("Unnamed ", what, " dimension cannot be identity-aligned.", call.=FALSE)
  } else stop("x must be a vector, matrix, or data frame.", call.=FALSE)
  .assert_ids(ids, paste("object", what, "IDs"))
  missing <- requested[is.na(match(requested, ids))]
  if (length(missing)) stop("Missing ", what, " identifiers: ", paste(missing, collapse=", "), call.=FALSE)
  extra <- ids[is.na(match(ids, requested))]
  if (length(extra) && !allow_extra) stop("Unexpected extra ", what, " identifiers: ", paste(extra, collapse=", "), call.=FALSE)
  idx <- match(requested, ids)
  if (is.atomic(x) && is.null(dim(x))) return(x[idx])
  if (margin == 1L) x[idx, , drop=FALSE] else x[, idx, drop=FALSE]
}
