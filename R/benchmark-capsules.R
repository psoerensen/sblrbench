# Internal mechanics for compact reference-capsule assembly.

benchmark_capsule_checksums <- function(path, files = NULL) {
  if (is.null(files))
    files <- sort(setdiff(list.files(path, recursive = FALSE), "checksums.csv"))
  paths <- file.path(path, files)
  info <- file.info(paths)
  data.frame(file = files, size_bytes = info$size,
    md5 = unname(benchmark_canonical_md5(paths)), stringsAsFactors = FALSE)
}

benchmark_capsule_staging_directory <- function(root, name) {
  path <- file.path(root, "promotion_staging",
    paste0(name, "-", Sys.getpid()))
  if (dir.exists(path))
    stop("Capsule promotion staging already exists: ", path,
      call. = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

benchmark_copy_capsule_files <- function(source, destination, files) {
  if (!dir.exists(destination))
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  source_paths <- file.path(source, files)
  if (any(!file.exists(source_paths)) ||
      !all(file.copy(source_paths, file.path(destination, files),
        overwrite = FALSE)))
    stop("Capsule source copying failed.", call. = FALSE)
  invisible(file.path(destination, files))
}

benchmark_promote_capsule <- function(staging, destination, validator) {
  if (!is.function(validator)) stop("validator must be a function.", call. = FALSE)
  validator(staging)
  if (dir.exists(destination))
    stop("Capsule destination already exists: ", destination,
      call. = FALSE)
  if (!file.rename(staging, destination))
    stop("Atomic capsule promotion failed: ", destination,
      call. = FALSE)
  invisible(destination)
}
