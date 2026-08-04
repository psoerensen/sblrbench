#!/usr/bin/env Rscript

# Optional, deliberately tiny installed-package interface smoke test.
# This is never called by validation-only mode or by report rendering.
run_study06_interface_fit_smoke_test <- function(
    store = file.path("results", "local", "study06_annotation_models",
      "_targets_convergence")) {
  pkgload::load_all(".", quiet = TRUE)
  source("studies/06_annotation_models/methods.R", local = TRUE)
  Sys.setenv(SBLR_BENCH_STUDY06_PHASE = "convergence",
    OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
    OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1")
  on.exit(Sys.unsetenv("SBLR_BENCH_STUDY06_PHASE"), add = TRUE)
  targets::tar_make(
    names = c(study06_convergence_stats, study06_alpha,
      study06_split, study06_annotation, study06_ld_bundle,
      study06_convergence_simulation, study06_method_specs),
    script = "studies/06_annotation_models/targets.R", store = store,
    callr_function = NULL, reporter = "verbose")
  read <- function(name) targets::tar_read_raw(name, store = store)
  cfg <- read("study06_config")
  cfg$convergence_pilot$nit <- 2L
  cfg$convergence_pilot$nburn <- 0L
  cfg$convergence_pilot$max_history <- 2L
  methods <- Filter(function(x) isTRUE(x$annotation),
    read("study06_method_specs"))
  simulation <- read("study06_convergence_simulation")
  stats <- read("study06_convergence_stats")$stats
  Glist <- read("study06_ld_bundle")$Glist
  split <- read("study06_split")
  A <- read("study06_annotation")
  alpha <- read("study06_alpha")
  runs <- lapply(methods, .study06_fit, simulation = simulation,
    stats = stats, Glist = Glist, split = split, A = A,
    alpha_bundle = alpha, config = cfg, phase = "convergence")
  failed <- vapply(runs, function(x) !identical(x$status, "ok"), logical(1))
  if (any(failed)) stop(paste(vapply(runs[failed], `[[`, "", "reason"),
    collapse = " | "), call. = FALSE)
  for (run in runs) {
    if (length(run$fit$chains) != 4L ||
        dim(run$fit$convergence_traces$values)[2L] != 4L ||
        !all(c("alpha", "sigmaSqAlpha") %in%
          run$fit$convergence_traces$quantities$parameter_name))
      stop("Tiny installed-package annotation output contract failed.",
        call. = FALSE)
  }
  data.frame(method = vapply(runs, function(x) x$method$id, ""),
    status = "passed", chain_count = 4L, retained_draws = 2L,
    runtime_seconds = vapply(runs, `[[`, numeric(1), "runtime"),
    stringsAsFactors = FALSE)
}

if (sys.nframe() == 0L) {
  result <- run_study06_interface_fit_smoke_test()
  dir.create("results/local/study06_annotation_models", recursive = TRUE,
    showWarnings = FALSE)
  write.csv(result,
    "results/local/study06_annotation_models/interface_fit_smoke_test.csv",
    row.names = FALSE)
  print(result)
}
