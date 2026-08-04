# Study 06 analysis entry point. The current scientific outcome is a stop.
# Do not start a full benchmark unless a new convergence design is prespecified.

load_study06_sources <- function(environment = new.env(parent = globalenv())) {
  files <- c("config.R", "annotation_design.R", "simulation.R", "chain_extraction.R",
    "diagnostics.R", "methods.R", "metrics.R", "promotion.R")
  for (file in files) sys.source(file.path("studies", "06_annotation_models", file), envir = environment)
  invisible(environment)
}

run_study06_convergence <- function() {
  Sys.setenv(SBLR_BENCH_STUDY = "06_annotation_models")
  targets::tar_make()
}

message("Study 06 helpers loaded. run_study06_convergence() reruns only when called; the published stop capsule remains authoritative.")
