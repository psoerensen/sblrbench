# Reusable structural checks for compact reference capsules.

benchmark_validate_capsule_checksums <- function(path, required,
                                                 allow_extra = FALSE) {
  if (!dir.exists(path)) stop("Capsule directory does not exist: ", path,
    call. = FALSE)
  missing <- setdiff(required, list.files(path, recursive = FALSE))
  if (length(missing))
    stop("Capsule is missing required files: ", paste(missing,
      collapse = ", "), call. = FALSE)
  inventory <- utils::read.csv(file.path(path, "checksums.csv"),
    stringsAsFactors = FALSE)
  expected_names <- c("file", "size_bytes", "md5")
  expected_files <- setdiff(required, "checksums.csv")
  invalid <- !identical(names(inventory), expected_names) ||
    anyNA(inventory$file) || anyDuplicated(inventory$file) ||
    any(inventory$file != basename(inventory$file)) ||
    any(grepl("(^[A-Za-z]:|^[/\\\\]|(^|[/\\\\])[.][.]([/\\\\]|$))",
      inventory$file)) || any(!grepl("^[0-9a-f]{32}$", inventory$md5))
  if (!allow_extra) invalid <- invalid ||
    !setequal(inventory$file, expected_files)
  if (invalid) stop("Invalid capsule checksum inventory.", call. = FALSE)
  paths <- file.path(path, inventory$file)
  if (any(!file.exists(paths)) ||
      any(unname(benchmark_canonical_md5(paths)) != inventory$md5))
    stop("Capsule checksum validation failed.", call. = FALSE)
  invisible(TRUE)
}
