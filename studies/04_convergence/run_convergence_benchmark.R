# Complete Study 04 development benchmark
# Four matched methods x four identifiable native chains = 16 chains.
# Public genotype inputs are pinned, checksum-validated, and cached. Existing
# valid targets and Study 03 genotype/LD preparation are reused.

Sys.setenv(
  SBLR_BENCH_STUDY = "04_convergence",
  OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"
)

targets::tar_outdated()
targets::tar_make()
targets::tar_read(convergence_recommendations)
targets::tar_read(convergence_checkpoint_summary)
targets::tar_read(convergence_chain_status)
targets::tar_meta(fields = c(name, seconds, error))
