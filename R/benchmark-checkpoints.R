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

benchmark_semantic_checkpoint_identity <- function(diagnostic_id,
                                                    scientific_inputs) {
  .benchmark_scalar_string(diagnostic_id, "diagnostic_id")
  if (!is.list(scientific_inputs) || is.null(names(scientific_inputs)) ||
      any(!nzchar(names(scientific_inputs))) ||
      anyDuplicated(names(scientific_inputs)))
    stop("scientific_inputs must be a uniquely named list.", call. = FALSE)
  forbidden <- c("source_hash", "source_file_hash", "source_files",
    "script_hash", "driver_hash", "script_path", "source_path",
    "report_path", "documentation_path", "working_directory", "timestamp")
  field_names <- function(x) {
    if (!is.list(x)) return(character())
    c(names(x), unlist(lapply(x, field_names), use.names = FALSE))
  }
  found <- intersect(unique(field_names(scientific_inputs)), forbidden)
  if (length(found))
    stop("Semantic checkpoint identity cannot contain: ",
      paste(found, collapse = ", "), ".", call. = FALSE)
  list(checkpoint_schema = "sblrbench-semantic-v2",
    diagnostic_id = diagnostic_id, scientific_inputs = scientific_inputs)
}

benchmark_semantic_checkpoint_hash <- function(identity) {
  if (!is.list(identity) ||
      !identical(identity$checkpoint_schema, "sblrbench-semantic-v2"))
    stop("Semantic checkpoint identity must use sblrbench-semantic-v2.",
      call. = FALSE)
  benchmark_hash_object(identity)
}

benchmark_load_semantic_checkpoint <- function(path, expected_hash,
                                                validator = NULL) {
  if (!file.exists(path))
    return(list(value = NULL, reused = FALSE, reason = "missing"))
  value <- try(readRDS(path), silent = TRUE)
  if (inherits(value, "try-error"))
    stop("Checkpoint cannot be read: ", path, call. = FALSE)
  if (!identical(value$checkpoint_schema, "sblrbench-semantic-v2"))
    stop(paste("Legacy source-hashed diagnostic checkpoint detected.",
      "This checkpoint schema has been retired and is not reusable under",
      "the shared semantic checkpoint framework."), call. = FALSE)
  benchmark_load_checkpoint(path, expected_hash, hash_field = "semantic_hash",
    validator = validator)
}
