study <- Sys.getenv("SBLR_BENCH_STUDY", "00_contract_smoke")
path <- file.path("studies", study, "targets.R")
if (!file.exists(path)) stop("Unknown study: ", study)
source(path, local = TRUE)$value
