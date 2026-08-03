# Study 05 analysis entry point. The current scientific outcome is a stop.
# Do not start a full benchmark unless a new convergence design is prespecified.

load_study05_sources <- function(environment = new.env(parent = globalenv())) {
  files <- c("config.R", "annotation_design.R", "simulation.R", "chain_extraction.R",
    "diagnostics.R", "methods.R", "metrics.R", "promotion.R")
  for (file in files) sys.source(file.path("studies", "05_annotation_models", file), envir = environment)
  invisible(environment)
}

run_study05_convergence <- function() {
  Sys.setenv(SBLR_BENCH_STUDY = "05_annotation_models")
  targets::tar_make()
}

message("Study 05 helpers loaded. run_study05_convergence() reruns only when called; the published stop capsule remains authoritative.")
