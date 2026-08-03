# Study 02 analysis entry point. Run from the repository root.
# Sourcing this file is side-effect free apart from defining functions/messages.

load_study02_sources <- function(environment = new.env(parent = globalenv())) {
  for (file in c("config.R", "pilot.R", "promotion.R"))
    sys.source(file.path("studies", "02_prediction", file), envir = environment)
  invisible(environment)
}

run_study02 <- function() {
  Sys.setenv(SBLR_BENCH_STUDY = "02_prediction")
  targets::tar_make()
}

message("Study 02 helpers loaded. Call load_study02_sources() or run_study02(); the website never calls either function.")
