# Study 03 analysis entry point. Run from the repository root.

load_study03_sources <- function(environment = new.env(parent = globalenv())) {
  files <- c("config.R", "simulation.R", "estimands.R", "methods.R", "metrics.R", "promotion.R")
  for (file in files) sys.source(file.path("studies", "03_parameter_estimation", file), envir = environment)
  invisible(environment)
}

run_study03 <- function() {
  Sys.setenv(SBLR_BENCH_STUDY = "03_parameter_estimation", SBLR_BENCH_PROFILE = "development")
  targets::tar_make()
}

message("Study 03 helpers loaded. run_study03() wraps the current targets pipeline; worked_parameter_estimation_example.R is the smaller interface example.")
