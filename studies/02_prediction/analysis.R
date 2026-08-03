# Study 02 exact analysis workflow. Run from the repository root.

find_repository_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "studies", "02_prediction", "spec.R")))
      return(path)
    parent <- dirname(path)
    if (identical(parent, path))
      stop("Could not locate the sblrbench repository root.", call. = FALSE)
    path <- parent
  }
}

root <- find_repository_root()
old_working_directory <- setwd(root)
on.exit(setwd(old_working_directory), add = TRUE)

# Setup and provenance -----------------------------------------------------
library(sblr)
library(sblrbench)
spec <- read_benchmark_spec(file.path("studies", "02_prediction", "spec.R"))

profile <- Sys.getenv("SBLR_BENCH_PROFILE", "benchmark")
output_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path("results", "local", "02_prediction"))

# Data ---------------------------------------------------------------------
# The shared execution path loads the pinned qgdata panel, preserves canonical
# sample/marker order, creates the fixed 70/30 split, and learns scaling, LD,
# and summary statistics from training individuals only.

# Scenarios ----------------------------------------------------------------
# `spec$scenarios` exposes the homogeneous sparse and variance-mixture designs.
# Both retain 50 causal markers, target h2 = 0.30, and the committed seed rules.

# Methods ------------------------------------------------------------------
# `spec$methods` exposes BED BayesC, BED BayesR, CSR SBayesC, and CSR SBayesR.
# Benchmark controls are the frozen Study 04 recommendations.

# Execution ----------------------------------------------------------------
results <- run_benchmark(spec = spec, output_dir = output_dir,
  profile = profile, resume = TRUE, validate_only = FALSE)

# Tables -------------------------------------------------------------------
fit_status <- results$status
prediction_metrics <- results$metrics
runtime <- results$runtime

# Plots --------------------------------------------------------------------
# Plotting belongs downstream of the tidy tables. The website report reads the
# frozen reference capsule and never executes this analysis.

# Extension points ---------------------------------------------------------
# Add supported scenarios, methods, or metrics by editing spec.R and the small
# task-specific shared functions. Do not encode new scientific choices inside
# run_benchmark().
