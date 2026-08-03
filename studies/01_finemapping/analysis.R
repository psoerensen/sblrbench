# Study 01 analysis entry point. Run from the repository root.
# This file defines helpers only: sourcing it does not simulate or fit anything.

load_study01_sources <- function(environment = new.env(parent = globalenv())) {
  for (file in c("config.R", "pilot.R", "promotion.R"))
    sys.source(file.path("studies", "01_finemapping", file), envir = environment)
  invisible(environment)
}

run_study01 <- function() {
  Sys.setenv(SBLR_BENCH_STUDY = "01_finemapping")
  targets::tar_make()
}

message("Study 01 helpers loaded. Call load_study01_sources() to inspect the implementation or run_study01() to run the targets pipeline.")
