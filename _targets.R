study <- Sys.getenv("SBLR_BENCH_STUDY", "00_contract_smoke")
if (study %in% c("01_finemapping", "02_prediction", "03_parameter_estimation",
    "04_convergence", "05_ld_operator"))
  stop("Migrated studies use scripts/run_benchmark.R, not a per-study targets graph.")
path <- file.path("studies", study, "targets.R")
if (!file.exists(path)) stop("Unknown study: ", study)
source(path, local = TRUE)$value
