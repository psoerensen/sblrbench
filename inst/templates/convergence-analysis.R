# Convergence-analysis template. Workshop settings exercise mechanics only and
# are unsuitable for method-performance or universal-convergence claims.

library(sblr)
library(sblrbench)

profile <- "workshop" # change to "benchmark" for the complete validated design
output_dir <- file.path("results", "local", "my-convergence-analysis")
spec <- read_benchmark_spec(file.path("studies", "04_convergence", "spec.R"))

# Extension points: choose from registered true-chain quantities, revise an
# explicitly justified threshold, or replace the matched grid in a copied spec.
# Do not infer traces from posterior means or final states.
spec$diagnostics$registry$required <-
  spec$diagnostics$registry$quantity %in% c("effect_variance",
    "genetic_variance", "residual_variance", "heritability")
spec$diagnostics$thresholds$rhat <- 1.01
spec$diagnostics$thresholds$ess_bulk <- 400
spec$diagnostics$thresholds$ess_tail <- 400
spec$diagnostics$thresholds$relative_mcse <- 0.05
validate_benchmark_spec(spec)

print(benchmark_convergence_design(spec))
print(benchmark_convergence_seeds(spec, profile))

results <- run_benchmark(spec, output_dir = output_dir, profile = profile,
  resume = TRUE)
convergence <- results$convergence
recommendations <- results$recommendations
print(recommendations)

selection <- convergence[convergence$stage == "selection", , drop = FALSE]
rhat_plot <- plot_convergence_rhat(selection,
  spec$diagnostics$thresholds$rhat)
ess_plot <- plot_convergence_ess(selection,
  spec$diagnostics$thresholds$ess_bulk)
print(rhat_plot)
print(ess_plot)

# Further extensions belong in a copied specification or focused table-only
# summaries, not in sampler, checkpoint, or trace-extraction code here.
