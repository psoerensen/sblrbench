# Study 03 exact parameter-estimation workflow. Run from repository root.
find_repository_root <- function(path=getwd()) {
  path <- normalizePath(path,winslash="/",mustWork=TRUE)
  repeat {
    if(file.exists(file.path(path,"DESCRIPTION")) &&
       file.exists(file.path(path,"studies","03_parameter_estimation","spec.R")))
      return(path)
    parent <- dirname(path)
    if(identical(parent,path)) stop("Could not locate sblrbench root.",call.=FALSE)
    path <- parent
  }
}
root <- find_repository_root(); old <- setwd(root); on.exit(setwd(old),add=TRUE)

# Setup and provenance -----------------------------------------------------
library(sblr)
library(sblrbench)
spec <- read_benchmark_spec(file.path("studies","03_parameter_estimation",
  "spec.R"))
profile <- Sys.getenv("SBLR_BENCH_PROFILE","benchmark")
output_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path("results","local","03_parameter_estimation"))

# Data ---------------------------------------------------------------------
# The shared path loads the pinned qgdata panel, preserves all 5,000 samples
# and canonical marker order, and constructs full-sample sparse LD.

# Scenarios ----------------------------------------------------------------
# sparse_homogeneous and sparse_mixture each retain 50 causal markers and
# realized h2=0.30 under the exact committed random-number ordering.

# Methods ------------------------------------------------------------------
# BED BayesC/BayesR and CSR SBayesC/SBayesR use Study 04's frozen controls.

# Estimands ----------------------------------------------------------------
# Six realized quantities are audited: causal proportion, nonzero-effect
# variance, total marker-effect variance, genetic variance, residual variance,
# and h2 = vgs/(vgs+ves). Nonlinear quantities are transformed draw by draw.

# Execution ----------------------------------------------------------------
results <- run_benchmark(spec,output_dir,profile=profile,resume=TRUE)

# Summaries ----------------------------------------------------------------
parameter_estimates <- results$estimates
parameter_metrics <- results$metrics
runtime <- results$runtime

# Plots --------------------------------------------------------------------
# Plot tidy result tables here. The website report reads only the frozen capsule.

# Extension points ---------------------------------------------------------
# Add supported scenarios, methods, estimands, or metrics in spec.R and the
# corresponding small shared task functions; keep scientific choices explicit.
