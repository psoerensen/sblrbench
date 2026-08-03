# Internal checkpoint mechanics. Scientific identity construction remains in studies.

benchmark_hash_object <- function(x, algorithm = "sha256") {
  digest::digest(x, algo = algorithm, serialize = TRUE)
}

benchmark_atomic_save_rds <- function(x, path, compress = FALSE,
                                      temporary_prefix = ".benchmark-") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(temporary_prefix, dirname(path), ".rds")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(x, temporary, compress = compress)
  if (!file.rename(temporary, path))
    stop("Atomic checkpoint replacement failed: ", path, call. = FALSE)
  invisible(path)
}

benchmark_load_checkpoint <- function(path, expected_hash = NULL,
                                      hash_field = "input_hash",
                                      validator = NULL) {
  if (!file.exists(path))
    return(list(value = NULL, reused = FALSE, reason = "missing"))
  value <- try(readRDS(path), silent = TRUE)
  if (inherits(value, "try-error"))
    stop("Checkpoint cannot be read: ", path, call. = FALSE)
  if (!is.null(expected_hash) &&
      !identical(value[[hash_field]], expected_hash))
    stop("Checkpoint input hash differs; refusing reuse: ", path,
      call. = FALSE)
  if (!is.null(validator) && !isTRUE(validator(value)))
    stop("Checkpoint validation failed; refusing reuse: ", path,
      call. = FALSE)
  list(value = value, reused = TRUE, reason = "validated")
}
