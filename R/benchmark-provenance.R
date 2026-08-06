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
  path <- find.package(package, lib.loc = lib.loc)
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE,
    all.files = TRUE, no.. = TRUE))
  files <- files[file.info(files)$isdir %in% FALSE]
  relative <- substring(normalizePath(files, winslash = "/"),
    nchar(normalizePath(path, winslash = "/")) + 2L)
  hashes <- if (length(files)) unname(benchmark_file_sha256(files)) else
    character()
  tree_sha <- digest::digest(paste(relative, hashes, sep = "=", collapse = "\n"),
    algo = "sha256", serialize = FALSE)
  list(package = package, version = description[["Version"]], sha = sha,
    remote_sha_available = !is.na(sha), installed_tree_sha256 = tree_sha,
    path = path,
    source_status = if (is.na(sha))
      "installed package; source commit unavailable" else
      "installed package; commit from package metadata")
}

#' Inspect a local source checkout used for an installed package
#'
#' @param path Source repository path.
#' @return A list containing the resolved path, branch, SHA, short status, and
#'   cleanliness flag.
#' @export
benchmark_source_provenance <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  git <- function(arguments) tryCatch(suppressWarnings(system2(
    "git", c("-C", normalized, arguments), stdout = TRUE, stderr = TRUE)),
    error = identity)
  sha <- git(c("rev-parse", "HEAD"))
  branch <- git(c("branch", "--show-current"))
  status <- git(c("status", "--short"))
  valid <- function(value) !inherits(value, "error") &&
    is.null(attr(value, "status"))
  if (!valid(sha) || length(sha) != 1L ||
      !grepl("^[0-9a-f]{40}$", sha))
    stop("Package source identity is unavailable at ", normalized, ".",
      call. = FALSE)
  if (!valid(branch) || length(branch) > 1L || !valid(status))
    stop("Package source Git status is unavailable at ", normalized, ".",
      call. = FALSE)
  list(path = normalized, branch = unname(branch %||% ""),
    sha = unname(sha), status = unname(status), clean = length(status) == 0L)
}

#' Resolve installed and local-source package provenance
#'
#' A local `R CMD INSTALL` commonly omits `RemoteSha`. In that case this
#' function records the clean source checkout identity and the installed tree
#' hash without inventing package metadata.
#'
#' @param package Package name.
#' @param source_path Local source checkout.
#' @param lib.loc Optional library location.
#' @param installation_command Exact local installation command when known.
#' @return Combined installed/source provenance.
#' @export
benchmark_resolve_package_provenance <- function(
    package, source_path, lib.loc = NULL, installation_command = NA_character_) {
  installed <- benchmark_package_provenance(package, lib.loc)
  source <- benchmark_source_provenance(source_path)
  list(package = package, version = installed$version,
    installed_path = installed$path,
    installed_tree_sha256 = installed$installed_tree_sha256,
    remote_sha = installed$sha,
    remote_sha_available = installed$remote_sha_available,
    source_path = source$path, source_branch = source$branch,
    source_sha = source$sha, source_clean = source$clean,
    source_status = source$status,
    installation_command = as.character(installation_command))
}

#' Validate a package identity for a new local experiment
#'
#' @param provenance Result of [benchmark_resolve_package_provenance()].
#' @param expected_source_sha Required clean source SHA.
#' @return `provenance`, invisibly.
#' @export
benchmark_assert_experiment_package <- function(provenance,
                                                expected_source_sha) {
  if (is.null(provenance$source_sha) || !nzchar(provenance$source_sha))
    stop("Package source identity is missing.", call. = FALSE)
  if (!isTRUE(provenance$source_clean))
    stop("Package source repository is dirty.", call. = FALSE)
  if (!identical(provenance$source_sha, expected_source_sha))
    stop("Package source SHA mismatch: ", provenance$source_sha,
      "; expected ", expected_source_sha, ".", call. = FALSE)
  if (isTRUE(provenance$remote_sha_available) &&
      !identical(provenance$remote_sha, expected_source_sha))
    stop("Installed package RemoteSha mismatch: ", provenance$remote_sha,
      "; expected ", expected_source_sha, ".", call. = FALSE)
  if (is.null(provenance$installed_tree_sha256) ||
      !grepl("^[0-9a-f]{64}$", provenance$installed_tree_sha256))
    stop("Installed package tree identity is missing.", call. = FALSE)
  invisible(provenance)
}

#' Validate a package SHA stored in historical evidence
#'
#' @param recorded_sha SHA stored in the artifact.
#' @param expected_sha Historical experiment pin.
#' @return `recorded_sha`, invisibly.
#' @export
benchmark_assert_recorded_package_sha <- function(recorded_sha,
                                                  expected_sha) {
  if (!identical(recorded_sha, expected_sha))
    stop("Recorded package SHA mismatch: ", recorded_sha, "; expected ",
      expected_sha, ".", call. = FALSE)
  invisible(recorded_sha)
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
