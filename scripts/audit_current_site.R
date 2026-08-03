html <- list.files("_site", pattern = "[.]html$", recursive = TRUE,
  full.names = TRUE)
missing <- character()
absolute_windows <- character()
for (file in html) {
  text <- paste(readLines(file, warn = FALSE), collapse = "\n")
  links <- regmatches(text,
    gregexpr('href="[^"]+"', text, perl = TRUE))[[1L]]
  if (!length(links) || identical(links[[1L]], -1L)) next
  links <- substring(links, 7L, nchar(links) - 1L)
  absolute_windows <- c(absolute_windows,
    links[grepl("^[A-Za-z]:[/\\\\]", links)])
  links <- links[!grepl("^(https?:|mailto:|javascript:|#|/)", links)]
  links <- sub("[?#].*$", "", links)
  links <- links[nzchar(links)]
  for (link in links) {
    target <- file.path(dirname(file), utils::URLdecode(link))
    if (dir.exists(target)) target <- file.path(target, "index.html")
    if (!file.exists(target)) missing <- c(missing,
      paste(file, link, sep = " -> "))
  }
}
missing <- unique(missing)
absolute_windows <- unique(absolute_windows)
cat("HTML", length(html), "MISSING", length(missing), "ABS_WINDOWS",
  length(absolute_windows), "\n")
if (length(missing)) cat(paste(missing, collapse = "\n"), "\n")
if (length(absolute_windows)) cat(paste(absolute_windows, collapse = "\n"),
  "\n")
if (length(missing) || length(absolute_windows)) quit(status = 1L)
