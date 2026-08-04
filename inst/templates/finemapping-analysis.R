# Fine-mapping analysis template. Workshop settings are for learning and code
# validation only; they are unsuitable for method-performance claims.
library(sblr)
library(sblrbench)

profile <- "workshop" # change to "benchmark" for a prespecified analysis
spec <- read_benchmark_spec("studies/01_finemapping/spec.R")

# Extension points: replace spec$data; define the locus and simulated causal
# variants in causal_design/locus_design; select supported methods and metrics.
spec$supported_profiles$workshop$replicate_count <- 1L
output_dir <- "results/local/my-finemapping-analysis"

print(benchmark_scenario_table(spec,profile))
print(benchmark_method_table(spec,profile))
results <- run_benchmark(spec,output_dir,profile,resume=TRUE)

marker_results <- results$marker_results
credible_sets <- results$credible_sets
causal_pip_plot <- plot_causal_marker_pip(marker_results)
credible_set_plot <- plot_credible_set_size(credible_sets)
print(causal_pip_plot)
print(credible_set_plot)
print(benchmark_output_inventory(results))
