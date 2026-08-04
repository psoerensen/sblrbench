study <- Sys.getenv("SBLR_BENCH_STUDY", "")
development_studies <- c("06_annotation_models", "07_mt_validation")

if (!nzchar(study))
  stop("Set SBLR_BENCH_STUDY to an active development study: ",
    paste(development_studies, collapse = " or "), ".", call. = FALSE)
if (!study %in% development_studies)
  stop("No targets graph is available for `", study,
    "`. Completed Studies 01-05 use scripts/run_benchmark.R; ",
    "Studies 06-07 remain development-only.", call. = FALSE)

source(file.path("studies", study, "targets.R"), local = TRUE)$value
