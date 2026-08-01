run_study07_contract_smoke_test <- function() {
  root <- if (file.exists(file.path("studies", "07_mtblr_validation",
    "config.R"))) "." else file.path("..", "..")
  config <- source(file.path(root, "studies", "07_mtblr_validation",
    "config.R"), local = TRUE)$value
  for (f in c("state_contract.R", "simulation.R", "alignment.R"))
    source(file.path(root, "studies", "07_mtblr_validation", f),
      local = TRUE)
  .study07_validate_state_contract(config)
  set.seed(7); Z <- scale(matrix(rnorm(300 * 80), 300, 80))
  rownames(Z) <- paste0("i", seq_len(nrow(Z)))
  colnames(Z) <- paste0("m", seq_len(ncol(Z)))
  for (a in config$contract_architectures) {
    old <- config$simulation$state_counts[[a]]
    config$simulation$state_counts[[a]] <- pmin(old, c(5L, 5L, 10L))
    x <- .study07_simulate(Z, a, 1L, config)
    .study07_validate_simulation(x, Z, config)
  }
  TRUE
}

if (sys.nframe() == 0L) {
  stopifnot(isTRUE(run_study07_contract_smoke_test()))
  cat("Study 07 sampler-free contract smoke test passed.\n")
}
