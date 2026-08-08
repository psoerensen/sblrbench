study <- Sys.getenv("SBLR_BENCH_STUDY", "")
development_studies <- "08_mt_validation"

if (!nzchar(study))
  stop("Set SBLR_BENCH_STUDY to an active development study: ",
    paste(development_studies, collapse = " or "), ".", call. = FALSE)
if (identical(study, "06_annotation_models"))
  stop("Study 06 uses scripts/run_benchmark.R in validate_only, qualification, or final mode; no targets graph remains.",
    call. = FALSE)
if (!study %in% development_studies)
  stop("No targets graph is available for `", study,
    "`. Completed Studies 01-05 use scripts/run_benchmark.R; ",
    "Study 06 uses the shared runner; Study 08 remains development-only.",
    call. = FALSE)

source(file.path("studies", study, "targets.R"), local = TRUE)$value
