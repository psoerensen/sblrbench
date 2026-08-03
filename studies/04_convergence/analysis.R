# Study 04 analysis entry point. Selection and validation remain separate stages.

load_study04_sources <- function(environment = new.env(parent = globalenv())) {
  files <- c("config.R", "chain_extraction.R", "diagnostic_registry.R", "diagnostics.R",
    "methods.R", "recommendations.R", "promotion.R")
  for (file in files) sys.source(file.path("studies", "04_convergence", file), envir = environment)
  invisible(environment)
}

run_study04_selection <- function() {
  Sys.setenv(SBLR_BENCH_STUDY = "04_convergence")
  targets::tar_make()
}

run_study04_validation <- function() {
  targets::tar_make(script = "studies/04_convergence/validation_targets.R")
}

message("Study 04 helpers loaded. Run selection and fixed-setting validation explicitly with the two run_study04_* functions.")
