# Complete Study 02 development benchmark.
# Public qgdata files are downloaded from the pinned revision, checksum
# validated, and cached under results/local/. Genotype scaling, LD, and summary
# statistics are learned from training individuals only; test phenotypes never
# enter fitting. targets reuses every valid upstream object, so do not delete
# _targets/ during routine reproduction.
Sys.setenv(
  SBLR_BENCH_STUDY = "02_prediction",
  SBLR_BENCH_REPLICATES = "10"
)

targets::tar_outdated()
targets::tar_make()

# Compact aggregate evidence.
targets::tar_read(prediction_method_summary)
targets::tar_read(prediction_paired_summary)
targets::tar_read(prediction_replicate_status)
targets::tar_meta(fields = c(name, status, seconds, error))

# Exact numerical equality additionally depends on recorded package versions,
# seeds, platform, BLAS, and other numerical libraries.
