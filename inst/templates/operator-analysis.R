# LD-operator analysis template
# Workshop output is for learning and structural validation only; it is not
# suitable for method-performance or operator-equivalence claims.

profile <- "workshop"                 # or "benchmark"
output_dir <- "results/local/my-operator-analysis"

library(sblr)
library(sblrbench)
spec <- read_benchmark_spec("studies/05_ld_operator/spec.R")
source("studies/05_ld_operator/operator-design.R")
options(sblrbench.ld_operator_runner = study05_reference_operator_runner)

# Extension points: replace the example data description, define reference
# and approximate operators, edit block sizes/eigen retention, and choose the
# existing comparison metrics. Preserve explicit marker order and provenance.
spec$operators$block$size <- 1000L
spec$operators$eigen$proportions["low_rank_0995"] <- 0.995
validate_benchmark_spec(spec)

print(benchmark_scenario_table(spec, profile))
print(study05_operator_table(spec))

results <- run_benchmark(spec, output_dir, profile = profile, resume = TRUE)

operator_error_plot <- plot_operator_errors(results$operator_comparisons)
retained_rank_plot <- plot_operator_retained_rank(results$eigenvalue_summary)
print(operator_error_plot)
print(retained_rank_plot)

# Further extensions: add an approximate operator to spec$operators, retain
# the reference operator unchanged, and add only scientifically justified
# metrics to the tidy extracted result tables.
