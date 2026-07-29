`%||%` <- function(x, y) if (is.null(x)) y else x

.assert_ids <- function(x, what) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop(what, " must be unique, non-missing, non-empty character identifiers.", call. = FALSE)
  x
}

.assert_scalar_string <- function(x, what) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop(what, " must be one non-empty string.", call. = FALSE)
  x
}

.as_named_matrix <- function(x, row_ids, col_ids, what) {
  if (is.null(x)) return(NULL)
  if (is.vector(x) && !is.list(x)) {
    if (is.null(names(x)) || length(col_ids) != 1L) stop(what, " cannot be converted unambiguously.", call. = FALSE)
    x <- matrix(x, ncol = 1L, dimnames = list(names(x), col_ids))
  }
  if (!is.matrix(x) || !is.numeric(x)) stop(what, " must be a numeric matrix.", call. = FALSE)
  if (!identical(dim(x), c(length(row_ids), length(col_ids)))) stop(what, " has incompatible dimensions.", call. = FALSE)
  if (!identical(rownames(x), row_ids) || !identical(colnames(x), col_ids)) stop(what, " must have the canonical row and column names.", call. = FALSE)
  if (any(!is.finite(x))) stop(what, " contains non-finite values.", call. = FALSE)
  x
}

.scenario_fields <- function(simulation) {
  list(study = simulation$scenario$study %||% NA_character_,
       scenario = simulation$scenario$architecture %||% NA_character_,
       replicate = simulation$scenario$replicate)
}

.metric_row <- function(simulation, result, trait, metric, value = NA_real_, status = "ok", reason = "") {
  s <- .scenario_fields(simulation)
  data.frame(study=s$study, scenario=s$scenario, replicate=as.integer(s$replicate),
             method_id=result$method_id, trait=trait, metric=metric,
             value=as.numeric(value), status=status, reason=reason,
             stringsAsFactors=FALSE, check.names=FALSE)
}
