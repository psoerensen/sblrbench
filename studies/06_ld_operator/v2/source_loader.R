.study06v2_load_pinned_sblr <- function(config, recompile = FALSE) {
  if (isTRUE(recompile))
    stop("Study 06 v2 uses the installed package; runtime recompilation is prohibited.",
      call. = FALSE)
  library_path <- normalizePath(config$installed_library, winslash = "/",
    mustWork = TRUE)
  .libPaths(unique(c(library_path, .libPaths())))
  if ("sblr" %in% loadedNamespaces()) {
    loaded_path <- normalizePath(getNamespaceInfo(asNamespace("sblr"),
      "path"), winslash = "/", mustWork = TRUE)
    expected_path <- normalizePath(file.path(library_path, "sblr"),
      winslash = "/", mustWork = TRUE)
    if (!identical(loaded_path, expected_path))
      stop("A non-pinned sblr namespace was loaded before Study 06 v2.",
        call. = FALSE)
  }
  if (!requireNamespace("sblr", quietly = TRUE))
    stop("Pinned installed sblr package is unavailable.", call. = FALSE)
  namespace_path <- getNamespaceInfo(asNamespace("sblr"), "path")
  api <- getExportedValue("sblr", "stblr_block_eigen")
  description <- packageDescription("sblr")
  if (!identical(as.character(packageVersion("sblr")),
      config$required_sblr_version) ||
      !identical(description[["RemoteSha"]], config$required_sblr_sha) ||
      !identical(normalizePath(namespace_path, winslash = "/"),
        normalizePath(file.path(library_path, "sblr"), winslash = "/")) ||
      !all(c("representation", "eigen_prop",
        "low_rank_residual_rebuild_every") %in% names(formals(api))))
    stop("Pinned installed sblr 0.2.0 contract failed.", call. = FALSE)
  invisible(list(path = normalizePath(namespace_path, winslash = "/"),
    library = library_path, version = config$required_sblr_version,
    sha = config$required_sblr_sha))
}
