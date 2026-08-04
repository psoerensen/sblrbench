#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
resume <- "--resume" %in% args
phase_arg <- grep("^--phase=", args, value = TRUE)
if (!length(phase_arg)) {
  i <- match("--phase", args)
  phase <- if (!is.na(i) && length(args) >= i + 1L) args[i + 1L] else "all"
} else phase <- sub("^--phase=", "", phase_arg[1L])
if (!phase %in% c("convergence", "benchmark", "all"))
  stop("--phase must be convergence, benchmark, or all.", call. = FALSE)

find_root <- function(path = getwd()) {
  path <- normalizePath(path, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "sblrbench.Rproj"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path))
      stop("Cannot locate the sblrbench repository root.", call. = FALSE)
    path <- parent
  }
}

root <- find_root()
setwd(root)
Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1")

local_dir <- file.path(root, "results", "local", "study06_annotation_models")
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(local_dir, "study06_status.csv")
summary_path <- file.path(local_dir, "study06_summary.txt")
preflight_path <- file.path(local_dir, "preflight.txt")
initial_head <- "765346d0d5fd5d63c77f9b537badd415a7a894bc"

empty_status <- function() data.frame(
  phase = character(), scenario = character(), replicate = integer(),
  method = character(), state = character(), start_time = character(),
  finish_time = character(), elapsed_seconds = numeric(),
  output_path = character(), validation_status = character(),
  error_message = character(), stringsAsFactors = FALSE)
status <- empty_status()
if (resume && file.exists(status_path)) {
  prior <- try(utils::read.csv(status_path, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(prior, "try-error") &&
      identical(names(prior), names(status))) status <- prior
}

atomic_csv <- function(x, path) {
  tmp <- tempfile(".study06-status-", dirname(path), ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  con <- file(tmp, "ab"); flush(con); close(con)
  if (!file.exists(tmp) || file.info(tmp)$size < 2L)
    stop("Atomic status temporary file was not written.", call. = FALSE)
  old <- paste0(path, ".previous")
  if (file.exists(old)) unlink(old)
  if (file.exists(path) && !file.rename(path, old))
    stop("Could not preserve the previous status file.", call. = FALSE)
  if (!file.rename(tmp, path)) {
    if (file.exists(old)) file.rename(old, path)
    stop("Atomic status replacement failed.", call. = FALSE)
  }
  if (file.exists(old)) unlink(old)
}

write_summary <- function() {
  lines <- c("Study 06 annotation-model run",
    paste("Updated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("PID:", Sys.getpid()), paste("Requested phase:", phase),
    paste("Validation only:", validate_only), paste("Resume:", resume), "",
    if (nrow(status)) apply(status, 1L, function(x)
      paste(x[c("phase", "scenario", "replicate", "method", "state",
        "validation_status", "error_message")], collapse = " | "))
    else "No phases have completed.")
  writeLines(lines, summary_path, useBytes = TRUE)
}

record <- function(phase_name, state, started, output_path = "",
                   validation_status = "", error = "",
                   scenario = "", replicate = NA_integer_, method = "") {
  finished <- Sys.time()
  status <<- rbind(status, data.frame(phase = phase_name,
    scenario = scenario, replicate = replicate, method = method, state = state,
    start_time = format(started, tz = "UTC", usetz = TRUE),
    finish_time = format(finished, tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
    output_path = output_path, validation_status = validation_status,
    error_message = error, stringsAsFactors = FALSE))
  atomic_csv(status, status_path)
  write_summary()
}

with_log <- function(phase_name, log_name, expr, output_path = "") {
  started <- Sys.time()
  path <- file.path(local_dir, log_name)
  con <- file(path, open = "at", encoding = "UTF-8")
  sink(con, type = "output", split = TRUE)
  sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(type = "output"); close(con) }, add = TRUE)
  cat("\n[", format(started, tz = "UTC", usetz = TRUE), "] ",
    phase_name, " started\n", sep = "")
  flush.console()
  ans <- tryCatch({ force(expr); list(ok = TRUE, error = "") },
    error = function(e) list(ok = FALSE, error = paste0(conditionMessage(e),
      " [", paste(deparse(conditionCall(e)), collapse = " "), "]")))
  cat("[", format(Sys.time(), tz = "UTC", usetz = TRUE), "] ",
    phase_name, if (ans$ok) " completed\n" else
      paste0(" failed: ", ans$error, "\n"), sep = "")
  flush.console()
  record(phase_name, if (ans$ok) "completed" else "failed", started,
    output_path, if (ans$ok) "passed" else "failed", ans$error)
  ans$ok
}

cfg <- source(file.path(root, "studies", "06_annotation_models", "config.R"),
  local = TRUE)$value
for (f in c("annotation_design.R", "simulation.R", "methods.R",
            "chain_extraction.R", "diagnostics.R", "metrics.R", "pilot.R",
            "promotion.R"))
  source(file.path(root, "studies", "06_annotation_models", f), local = TRUE)

active_pid <- function(pid) {
  if (is.na(pid) || pid == Sys.getpid()) return(FALSE)
  code <- suppressWarnings(system2("powershell", c("-NoProfile", "-Command",
    sprintf("if (Get-Process -Id %d -ErrorAction SilentlyContinue) { exit 0 } else { exit 1 }",
      pid)), stdout = FALSE, stderr = FALSE))
  identical(code, 0L)
}

preflight <- function() {
  required <- c("sblr", "sblrbench", "targets", "posterior", "jsonlite",
    "ggplot2", "knitr", "rmarkdown", "testthat", "pkgload", "devtools")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) stop("Missing required installed packages: ",
    paste(names(available)[!available], collapse = ", "), call. = FALSE)
  pkgload::load_all(root, quiet = TRUE)
  if (!identical(normalizePath(root), normalizePath(getwd())))
    stop("Working directory is outside the repository root.", call. = FALSE)
  head <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  branch <- system2("git", c("branch", "--show-current"), stdout = TRUE)
  staged <- system2("git", c("diff", "--cached", "--name-only"), stdout = TRUE)
  if (!identical(head, initial_head))
    stop("HEAD changed since the clean initial gate.", call. = FALSE)
  if (length(staged)) stop("Staged changes are present.", call. = FALSE)
  if (!identical(as.character(packageVersion("sblr")), "0.1.2"))
    stop("Unexpected installed sblr version.", call. = FALSE)
  data_dirs <- c(file.path(root, "results/local/02_prediction/data"),
    file.path(root, "results/local/03_parameter_estimation/data"))
  for (dir in data_dirs) for (nm in cfg$example_data$files) {
    path <- file.path(dir, nm)
    if (!file.exists(path) ||
        unname(file.info(path)$size) != unname(cfg$example_data$size_bytes[nm]) ||
        unname(tools::md5sum(path)) != unname(cfg$example_data$md5[nm]))
      stop("Cached qgdata validation failed: ", path, call. = FALSE)
  }
  pid_file <- file.path(local_dir, "process_id.txt")
  if (file.exists(pid_file)) {
    pid <- suppressWarnings(as.integer(readLines(pid_file, n = 1L, warn = FALSE)))
    if (active_pid(pid)) stop("Another Study 06 process is active: ", pid,
      call. = FALSE)
  }
  paths <- .study06_paths()
  example_files <- sblrbench:::benchmark_example_files(paths$data_dir,
    list(example_data=cfg$example_data))
  base_glist <- sblrbench:::benchmark_load_glist(paths, example_files)
  marker_ids <- sblrbench:::benchmark_filter_markers(base_glist,cfg$chr,
    cfg$qc,cfg$sparse_ld)$marker_ids
  A <- .study06_annotation_design(marker_ids, cfg)
  alpha <- .study06_true_alpha(A, cfg)
  .study06_validate_annotation(A, marker_ids, cfg)
  for (scenario in cfg$scenarios)
    .study06_marker_probabilities(A, alpha[[scenario]], cfg$mixture_var)
  stopifnot(isTRUE(run_study06_contract_smoke_test()))
  for (p in c(file.path(local_dir, "_targets_convergence"),
              file.path(local_dir, "_targets_benchmark"),
              file.path(local_dir, "convergence_output"),
              file.path(local_dir, "benchmark_output"))) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    probe <- file.path(p, ".write-probe")
    writeLines("ok", probe); if (!file.exists(probe)) stop("Unwritable path: ", p)
    unlink(probe)
  }
  for (which_phase in c("convergence", "benchmark")) {
    Sys.setenv(SBLR_BENCH_STUDY06_PHASE = which_phase)
    manifest <- targets::tar_manifest(
      script = "studies/06_annotation_models/targets.R")
    if (!nrow(manifest)) stop("Empty target graph for ", which_phase)
  }
  Sys.unsetenv("SBLR_BENCH_STUDY06_PHASE")
  lines <- c("Study 06 preflight", "Result: PASSED",
    paste("Time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Repository:", normalizePath(root, winslash = "/")),
    paste("Branch:", branch), paste("HEAD:", head),
    "Initial working tree: clean",
    paste("Current unstaged Study 06 implementation files:",
      length(system2("git", c("status", "--porcelain"), stdout = TRUE))),
    paste("sblr:", packageVersion("sblr")),
    paste("sblr source commit:", packageDescription("sblr")$RemoteSha),
    paste("sblrbench:", packageVersion("sblrbench")),
    paste("qgdata commit:", cfg$example_data$commit),
    "qgdata: both cached copies match pinned sizes and checksums",
    paste("annotation dimensions:", nrow(A), "x", ncol(A)),
    paste("annotation rank:", qr(A)$rank),
    paste("expected non-null:", alpha$expected_nonnull),
    paste("expected enriched non-null share:",
      alpha$enriched_expected_nonnull_share),
    "tarchetypes: not required or checked; Study 06 has no tarchetypes calls",
    "thread controls: OMP=4, OMP limit=4, BLAS/MKL/vecLib=1",
    "target stores and output paths: writable",
    "model interface and sampler-free output contracts: passed",
    "network access: not used", "package installation: not attempted")
  cat(paste(lines, collapse = "\n"), "\n")
  invisible(TRUE)
}

# Source after helper definitions so validation-only can call it without sampling.
source(file.path(root, "studies", "06_annotation_models",
  "contract_smoke_test.R"), local = TRUE)

ok <- with_log("preflight", "preflight.txt", preflight(), preflight_path)
if (validate_only)
  quit(status = if (ok) 0L else 1L, save = "no")
if (!ok) quit(status = 1L, save = "no")

run_targets_phase <- function(which_phase) {
  capsule <- if (which_phase == "convergence")
    cfg$convergence_capsule else cfg$benchmark_capsule
  validator <- if (which_phase == "convergence")
    .study06_validate_convergence_capsule else .study06_validate_benchmark_capsule
  if (dir.exists(capsule)) {
    validator(capsule)
    message("Validated existing frozen capsule; target execution skipped: ",
      capsule)
    return(invisible(capsule))
  }
  output <- file.path(local_dir, paste0(which_phase, "_output"))
  store <- file.path(local_dir, paste0("_targets_", which_phase))
  Sys.setenv(SBLR_BENCH_STUDY06_PHASE = which_phase,
    SBLR_BENCH_OUTPUT_DIR = output)
  on.exit(Sys.unsetenv(c("SBLR_BENCH_STUDY06_PHASE",
    "SBLR_BENCH_OUTPUT_DIR")), add = TRUE)
  targets::tar_make(script = "studies/06_annotation_models/targets.R",
    store = store, callr_function = NULL, reporter = "verbose",
    seconds_meta_append = 1)
  if (which_phase == "convergence")
    .study06_promote_convergence(output, cfg)
  else .study06_promote_benchmark(output, cfg)
  validator(capsule)
  capsule
}

phase_ok <- logical()
if (phase %in% c("convergence", "all")) {
  phase_ok["convergence"] <- with_log("convergence", "convergence.log",
    run_targets_phase("convergence"), cfg$convergence_capsule)
  if (!phase_ok["convergence"]) quit(status = 1L, save = "no")
}
if (phase %in% c("benchmark", "all")) {
  if (!dir.exists(cfg$convergence_capsule))
    stop("Validated convergence capsule is required before the benchmark.",
      call. = FALSE)
  .study06_validate_convergence_capsule(cfg$convergence_capsule)
  phase_ok["benchmark"] <- with_log("benchmark", "benchmark.log",
    run_targets_phase("benchmark"), cfg$benchmark_capsule)
  if (!phase_ok["benchmark"]) quit(status = 1L, save = "no")
}

if (phase %in% c("benchmark", "all")) {
  phase_ok["render_and_tests"] <- with_log("render_and_tests",
    "render_and_tests.log", {
      quarto <- Sys.which("quarto")
      if (!nzchar(quarto)) {
        candidate <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
        if (!file.exists(candidate)) stop("Quarto executable not found.")
        quarto <- candidate
      }
      code <- system2(quarto, "render")
      if (!identical(code, 0L)) stop("quarto render failed: ", code)
      testthat::test_file("tests/testthat/test-study06-annotation-models.R",
        reporter = "summary", stop_on_failure = TRUE)
      testthat::test_local(root, reporter = "summary", stop_on_failure = TRUE)
      check <- devtools::check(root, document = FALSE, error_on = "warning")
      if (length(check$errors) || length(check$warnings))
        stop("Package check reported errors or warnings.")
      diff_check <- system2("git", c("diff", "--check"))
      if (!identical(diff_check, 0L)) stop("git diff --check failed.")
    }, root)
}

quit(status = if (all(phase_ok)) 0L else 1L, save = "no")
