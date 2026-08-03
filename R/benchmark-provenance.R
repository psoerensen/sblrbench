# Internal provenance and hashing mechanics shared by benchmark studies.

benchmark_git_sha <- function(path = ".", warn = TRUE) {
  normalized <- normalizePath(path, mustWork = FALSE)
  out <- tryCatch(suppressWarnings(system2("git",
    c("-C", normalized, "rev-parse", "HEAD"), stdout = TRUE,
    stderr = TRUE)), error = identity)
  failed <- inherits(out, "error") || !is.null(attr(out, "status")) ||
    length(out) != 1L || !grepl("^[0-9a-f]{40}$", out)
  if (!failed) return(unname(out))
  if (isTRUE(warn)) warning("Git commit metadata is unavailable for ", path,
    ".", call. = FALSE)
  NA_character_
}

#' Inspect installed package provenance
#'
#' @param package Package name.
#' @param lib.loc Optional library location.
#' @return A list containing package, version, source SHA, path, and status.
#' @export
benchmark_package_provenance <- function(package, lib.loc = NULL) {
  description <- utils::packageDescription(package, lib.loc = lib.loc)
  sha <- description[["RemoteSha"]]
  if (is.null(sha) || !length(sha) || is.na(sha) || !nzchar(sha))
    sha <- description[["GithubSHA1"]]
  if (is.null(sha) || !length(sha) || is.na(sha) || !nzchar(sha))
    sha <- NA_character_
  list(package = package, version = description[["Version"]], sha = sha,
    path = find.package(package, lib.loc = lib.loc),
    source_status = if (is.na(sha))
      "installed package; source commit unavailable" else
      "installed package; commit from package metadata")
}

benchmark_assert_package_sha <- function(package, expected_sha,
                                         lib.loc = NULL) {
  provenance <- benchmark_package_provenance(package, lib.loc = lib.loc)
  if (!identical(provenance$sha, expected_sha))
    stop("Installed ", package, " SHA mismatch: ", provenance$sha,
      "; expected ", expected_sha, ".", call. = FALSE)
  invisible(provenance)
}

benchmark_canonical_md5 <- function(paths) {
  paths <- as.character(paths)
  if (!length(paths)) return(setNames(character(), character()))
  if (anyNA(paths) || any(!file.exists(paths)))
    stop("All checksum paths must identify existing files.", call. = FALSE)
  text_extensions <- c("r", "qmd", "md", "csv", "json", "txt", "yml", "yaml")
  output <- character(length(paths))
  names(output) <- paths
  for (i in seq_along(paths)) {
    path <- paths[[i]]
    if (!tolower(tools::file_ext(path)) %in% text_extensions) {
      output[[i]] <- unname(tools::md5sum(path))
      next
    }
    bytes <- readBin(path, what = "raw", n = file.info(path)$size)
    values <- as.integer(bytes)
    if (length(values)) {
      is_cr <- values == 13L
      cr_before_lf <- is_cr & c(values[-1L] == 10L, FALSE)
      values[is_cr] <- 10L
      bytes <- as.raw(values[!cr_before_lf])
    }
    temporary <- tempfile("sblrbench-canonical-", fileext = ".bin")
    on.exit(unlink(temporary), add = TRUE)
    connection <- file(temporary, open = "wb")
    tryCatch(writeBin(bytes, connection), finally = close(connection))
    output[[i]] <- unname(tools::md5sum(temporary))
    unlink(temporary)
  }
  output
}

benchmark_file_sha256 <- function(paths) {
  paths <- as.character(paths)
  if (anyNA(paths) || any(!file.exists(paths)))
    stop("All checksum paths must identify existing files.", call. = FALSE)
  stats::setNames(vapply(paths, function(path) digest::digest(file = path,
    algo = "sha256", serialize = FALSE), character(1)), paths)
}

benchmark_session_information <- function() capture.output(utils::sessionInfo())
