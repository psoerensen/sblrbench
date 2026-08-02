#!/usr/bin/env Rscript
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config <- source(file.path(root, "studies", "06_ld_operator", "v2",
  "config.R"), local = TRUE)$value
source(file.path(root, "studies", "06_ld_operator", "v2", "guard.R"),
  local = TRUE)

grid <- .study06v2_validate_grid(config)
stopifnot(nrow(grid) == 60L)
for (id in grep("^low_rank_", config$configurations, value = TRUE)) {
  controls <- .study06v2_low_rank_configuration(id, config)
  stopifnot(identical(controls$representation, "low_rank"),
    is.finite(controls$eigen_prop), controls$eigen_prop < 1,
    !any(c("eigen_filter", "eigen_tau", "eigen_eta") %in% names(controls)))
}
cat("Study 06 v2 retained-low-rank contract smoke test: PASSED\n")
