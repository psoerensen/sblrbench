.study04_specs <- function(config) {
  n <- if (identical(config$profile, "five_replicate_validation")) 5L else 1L
  out <- list()
  for (r in seq_len(n)) for (i in seq_len(nrow(config$matched_grid))) {
    z <- as.list(config$matched_grid[i, ])
    z$replicate <- r
    out[[length(out) + 1L]] <- z
  }
  out
}

.study04_sim_specs <- function(config) {
  architectures <- names(config$simulation$architectures)
  n <- if (identical(config$profile, "five_replicate_validation")) 5L else 1L
  out <- list()
  for (i in seq_along(architectures)) for (r in seq_len(n)) out[[length(out) + 1L]] <- list(
    architecture = architectures[i], replicate = r,
    simulation_seed = as.integer(config$simulation$base_seed + i * 1000L + r))
  out
}

.study04_diagnostic_grid <- function(draws, config) do.call(rbind,
  lapply(config$burnin_candidates, function(b) do.call(rbind,
    lapply(config$retained_draw_candidates[config$retained_draw_candidates + b <= 3000L],
      function(n) .study04_diagnostics(draws, b, n, .study04_registry(), config$thresholds)))))

.study04_write <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); write.csv(x, path, row.names = FALSE); path }
