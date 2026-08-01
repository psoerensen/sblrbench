.study06v2_load_pinned_sblr <- function(config, recompile = FALSE) {
  snapshot <- normalizePath(config$source_snapshot, winslash = "/",
    mustWork = TRUE)
  pkgload::load_all(snapshot, recompile = recompile, quiet = TRUE,
    export_all = FALSE, helpers = FALSE)
  namespace_path <- getNamespaceInfo(asNamespace("sblr"), "path")
  api <- getExportedValue("sblr", "stblr_block_eigen")
  if (!identical(as.character(packageVersion("sblr")),
      config$required_sblr_version) ||
      !identical(normalizePath(namespace_path, winslash = "/"), snapshot) ||
      !all(c("representation", "eigen_prop") %in% names(formals(api))))
    stop("Pinned sblr 0.2.0 local-source contract failed.", call. = FALSE)
  invisible(list(path = snapshot, version = config$required_sblr_version,
    sha = config$required_sblr_sha))
}
