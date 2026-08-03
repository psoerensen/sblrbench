# Reusable single-trait prediction analysis template.
# Workshop output checks mechanics only and is unsuitable for method-performance
# claims. Use the benchmark profile only with an audited scientific spec.

library(sblr)
library(sblrbench)

profile <- "workshop"             # Change to "benchmark" for the full profile.
output_dir <- file.path("results", "local", "my_prediction_analysis")

# Start from the committed prediction contract, then review every scientific
# field before changing it.
spec <- read_benchmark_spec(file.path("studies", "02_prediction", "spec.R"))

# Replace example data: edit spec$data and supply a compatible Glist through
# SBLR_BENCH_GLIST. Arbitrary data formats are intentionally not inferred.
# Add scenarios: extend spec$scenarios with a supported effect distribution.
# Add methods: extend the supported plain method lists and dispatch explicitly.
# Choose metrics: edit spec$metrics using supported prediction metric names.
# Change output paths: edit output_dir above.
print(benchmark_scenario_table(spec, profile))
print(benchmark_method_table(spec, profile))

results <- run_benchmark(spec = spec, output_dir = output_dir,
  profile = profile, resume = TRUE, validate_only = FALSE)

print(results$status)
print(results$metrics)
prediction_plot <- plot_prediction_metrics(results$metrics)
runtime_plot <- plot_benchmark_runtime(results$runtime)
print(prediction_plot)
print(runtime_plot)
