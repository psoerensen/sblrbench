# Exact workflow entry point for the separated development benchmark.
# Genotype input is not redistributed. This reproduces the workflow only when
# the expected input is available; a compatible Glist may be substituted for testing.
# The frozen reference files reproduce the report without genotype data.
Sys.setenv(
  SBLR_BENCH_STUDY = "01_finemapping",
  SBLR_BENCH_REPLICATES = "10"
)

targets::tar_make()

targets::tar_manifest()
targets::tar_read(marker_metrics)
targets::tar_read(credible_set_summary)
targets::tar_read(computational_summary)
targets::tar_read(replicate_status)
