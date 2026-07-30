.study04_specs <- function(config) lapply(seq_len(nrow(config$matched_grid)), function(i)
  as.list(config$matched_grid[i, ]))

.study04_sim_specs <- function(config) {
  architectures <- names(config$simulation$architectures)
  lapply(seq_along(architectures), function(i) list(
    architecture = architectures[i], replicate = 1L,
    simulation_seed = as.integer(config$simulation$base_seed + i * 1000L + 1L)
  ))
}

.study04_diagnostic_grid <- function(draws, config) do.call(rbind,
  lapply(config$burnin_candidates, function(b) do.call(rbind,
    lapply(config$retained_draw_candidates[config$retained_draw_candidates + b <= 3000L],
      function(n) .study04_diagnostics(draws, b, n, .study04_registry(), config$thresholds)))))

.study04_write <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); write.csv(x, path, row.names = FALSE); path }
