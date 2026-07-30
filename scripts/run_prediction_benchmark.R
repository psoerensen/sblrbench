# Study 02: single-trait BayesC/BayesR prediction benchmark.
#
# Start with one replicate. Only ST-BED BayesC/BayesR and ST-CSR
# SBayesC/SBayesR are active. Public qgdata files are pinned, checksum
# validated and cached; scaling, sparse LD and summary statistics use training
# individuals only. Test phenotypes are evaluation-only.
Sys.setenv(
  SBLR_BENCH_STUDY = "02_prediction",
  SBLR_BENCH_REPLICATES = "1"
)

targets::tar_outdated()
targets::tar_make()
targets::tar_read(prediction_replicate_status)
targets::tar_read(prediction_metrics)
targets::tar_read(prediction_paired_method_differences)
targets::tar_read(prediction_simulation_summary)
targets::tar_meta(fields = c(name, seconds, error))

# After the one-replicate gate passes, expand without deleting _targets/:
# Sys.setenv(SBLR_BENCH_REPLICATES = "5")
# targets::tar_make()
# Sys.setenv(SBLR_BENCH_REPLICATES = "10")
# targets::tar_make()
#
# Changing the count adds missing branches and reuses genotype, split, scaling,
# and training-LD targets. Exact numbers can vary slightly across platforms and
# numerical libraries.
