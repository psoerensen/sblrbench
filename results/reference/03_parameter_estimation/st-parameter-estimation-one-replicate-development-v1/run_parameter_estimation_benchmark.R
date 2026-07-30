# Complete one-replicate development benchmark. The pinned public genotype
# files are downloaded and checksum-validated automatically and cached locally.
Sys.setenv(SBLR_BENCH_STUDY = "03_parameter_estimation", SBLR_BENCH_PROFILE = "development")
targets::tar_outdated()
targets::tar_make()
targets::tar_read(parameter_recovery_summary)
targets::tar_read(parameter_paired_summary)
targets::tar_read(parameter_computational_summary)
