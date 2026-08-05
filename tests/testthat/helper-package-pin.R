current_sblr_description_sha <- function() {
  description <- read.dcf(testthat::test_path("..", "..", "DESCRIPTION"),
    fields = "Remotes")[[1L]]
  remotes <- trimws(strsplit(description, ",", fixed = TRUE)[[1L]])
  remote <- remotes[grepl("^psoerensen/sblr@[0-9a-f]{40}$", remotes)]
  if (length(remote) != 1L)
    stop("DESCRIPTION must contain exactly one pinned psoerensen/sblr remote.",
      call. = FALSE)
  sub("^psoerensen/sblr@", "", remote)
}

with_current_sblr_description_pin <- function(spec) {
  spec$packages$sblr$sha <- current_sblr_description_sha()
  spec
}
