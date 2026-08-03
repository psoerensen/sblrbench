# Reusable scalar parameter-estimation template. Workshop runs check mechanics
# only and are unsuitable for method-performance or interval-coverage claims.
library(sblr)
library(sblrbench)

profile <- "workshop" # Change to "benchmark" for the audited full profile.
output_dir <- file.path("results","local","my_parameter_estimation")
spec <- read_benchmark_spec(file.path("studies","03_parameter_estimation",
  "spec.R"))

# Replace example data through spec$data and a compatible SBLR_BENCH_GLIST.
# Extend spec$scenarios or supported spec$methods explicitly.
# Select existing scalar parameters in spec$estimands and metrics in spec$metrics.
# Change output_dir above; do not hide scientific settings inside the runner.
results <- run_benchmark(spec,output_dir,profile=profile,resume=TRUE)
print(results$status)
print(results$estimates)
