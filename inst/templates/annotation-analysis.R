# Annotation-informed scalar-model analysis template
library(sblr)
library(sblrbench)

# Safe default. Use "qualification" first; use "final" only after it passes.
mode <- Sys.getenv("SBLR_BENCH_MODE", "validate_only")
profile <- "workshop"                 # or "benchmark"
output_dir <- "results/local/my-annotation-analysis"

source("studies/06_annotation_models/spec.R", local = TRUE)
source("studies/06_annotation_models/annotation-design.R", local = TRUE)

# Extension points ---------------------------------------------------------
# Replace example data and edit the spec's data provenance deliberately.
# Define deterministic informative and uninformative annotations here.
# Keep marker rows aligned and the intercept first.
# Add only methods exposing comparable retained annotation traces.
# Select annotation-prior and marker-recovery metrics in spec$metrics.

validate_benchmark_spec(spec)
marker_ids <- sprintf("marker_%05d", seq_len(spec$data$expected_marker_count))
annotations <- construct_annotation_design(marker_ids, spec)
annotation_truth <- construct_annotation_truth(annotations, spec)
print(annotation_design_summary(annotations, annotation_truth, spec))
print(benchmark_annotation_seeds(spec, profile, mode = "qualification"))

results <- run_benchmark(spec, output_dir = output_dir, profile = profile,
  resume = TRUE, validate_only = identical(mode, "validate_only"), mode = mode)

annotation_prior_plot <- causal_marker_plot <- NULL
if (is.data.frame(results$metrics) && nrow(results$metrics)) {
  annotation_prior_plot <- plot_annotation_prior_recovery(results$metrics)
  causal_marker_plot <- plot_annotation_marker_recovery(results$metrics)
  print(annotation_prior_plot)
  print(causal_marker_plot)
}

print(benchmark_output_inventory(results), row.names = FALSE)

# Workshop settings are for learning and interface checks, not performance
# claims. Qualification must pass under the prespecified benchmark design before
# final fitting; final execution is always a separate explicit command.
