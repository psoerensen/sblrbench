# Study 06 v2 analysis entry point. Run from the repository root.
# The report consumes the current capsule and never sources this script.

load_study06_sources <- function(environment = new.env(parent = globalenv())) {
  base <- file.path("studies", "06_ld_operator", "v2")
  files <- c("config.R", "source_loader.R", "guard.R", "design_crosswalk.R",
    "operator_validation.R", "deterministic_run.R", "methods.R", "runtime_data.R",
    "final_validation.R", "phases.R", "provenance.R", "promotion.R", "verification.R")
  for (file in files) sys.source(file.path(base, file), envir = environment)
  invisible(environment)
}

run_study06_validation <- function(extra_args = character()) {
  script <- file.path("results", "reference", "06_ld_operator", "current", "reproduce.R")
  system2(file.path(R.home("bin"), "Rscript"), c(script, "--validate-only", extra_args))
}

message("Study 06 helpers loaded. Inspect with load_study06_sources(); run_study06_validation() invokes the capsule's validation-only entry point.")
