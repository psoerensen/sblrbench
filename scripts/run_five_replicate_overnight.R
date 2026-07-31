#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
skip_check <- "--skip-check" %in% args

find_root <- function(path = getwd()) {
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "_targets.R"))) return(normalizePath(path))
    parent <- dirname(path)
    if (identical(parent, path)) stop("Cannot locate sblrbench repository root.", call. = FALSE)
    path <- parent
  }
}

root <- find_root()
setwd(root)
Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1")

base_dir <- file.path(root, "results", "local", "five_replicate_overnight")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
run_dir <- file.path(base_dir, paste0("run-", stamp, "-", Sys.getpid()))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(base_dir, "overnight_status.csv")
summary_path <- file.path(base_dir, "overnight_summary.txt")
preflight_path <- file.path(base_dir, "preflight.txt")

status <- data.frame(phase = character(), status = character(), started_at = character(),
  ended_at = character(), elapsed_seconds = numeric(), message = character(),
  log = character(), stringsAsFactors = FALSE)

write_atomic_csv <- function(x, path) {
  tmp <- tempfile("status-", tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.exists(tmp) || file.info(tmp)$size <= 0L ||
      inherits(try(utils::read.csv(tmp), silent = TRUE), "try-error"))
    stop("Status update was not fully written; preserving the prior status file.", call. = FALSE)
  backup <- paste0(path, ".previous")
  if (file.exists(backup)) unlink(backup)
  if (file.exists(path) && !file.rename(path, backup))
    stop("Could not preserve prior status file.", call. = FALSE)
  if (!file.rename(tmp, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("Atomic status update failed; prior status was restored.", call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup)
}

write_summary <- function() {
  lines <- c("Five-replicate overnight run", paste("Run directory:", run_dir),
    paste("Updated:", format(Sys.time(), tz = "UTC", usetz = TRUE)), "",
    if (nrow(status)) apply(status, 1L, function(z)
      paste(z[["phase"]], z[["status"]], z[["message"]], sep = " | ")) else "No phases recorded.")
  writeLines(lines, summary_path)
}

record <- function(phase, phase_status, started, message = "", log = "") {
  ended <- Sys.time()
  status <<- rbind(status, data.frame(phase = phase, status = phase_status,
    started_at = format(started, tz = "UTC", usetz = TRUE),
    ended_at = format(ended, tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
    message = message, log = log, stringsAsFactors = FALSE))
  write_atomic_csv(status, status_path)
  write_summary()
}

run_phase <- function(phase, log_name, code) {
  started <- Sys.time(); log_path <- file.path(base_dir, log_name)
  con <- file(log_path, open = "at", encoding = "UTF-8")
  sink(con, type = "output", split = TRUE); sink(con, type = "message")
  on.exit({ sink(type = "message"); sink(type = "output"); close(con) }, add = TRUE)
  cat("\n[", format(started, tz = "UTC", usetz = TRUE), "] phase ", phase, " started\n", sep = "")
  flush.console()
  result <- tryCatch({ force(code); list(ok = TRUE, message = "completed") },
    error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  cat("[", format(Sys.time(), tz = "UTC", usetz = TRUE), "] phase ", phase, " ",
    if (result$ok) "completed" else paste("failed:", result$message), "\n", sep = "")
  flush.console()
  record(phase, if (result$ok) "success" else "failed", started,
    result$message, normalizePath(log_path, winslash = "/", mustWork = FALSE))
  result$ok
}

preflight <- function() {
  required <- c("sblr", "targets", "posterior", "jsonlite", "ggplot2", "knitr",
    "rmarkdown", "testthat", "pkgload", "devtools")
  available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  if (!all(available)) stop("Missing required offline packages: ",
    paste(names(available)[!available], collapse = ", "), call. = FALSE)
  active_rel <- c("_targets.R", "studies/02_prediction/targets.R",
    "studies/03_parameter_estimation/targets.R", "studies/04_convergence/targets.R",
    "studies/04_convergence/validation_targets.R",
    "scripts/run_five_replicate_overnight.R",
    "scripts/run_five_replicate_overnight.ps1")
  source_files <- file.path(root, active_rel)
  rel <- substring(normalizePath(source_files, winslash = "/"), nchar(normalizePath(root, winslash = "/")) + 2L)
  keep <- file.exists(source_files)
  call_patterns <- c(paste0("tar", "chetypes::"),
    paste0("library(", "tar", "chetypes)"),
    paste0("require(", "tar", "chetypes)"))
  hits <- unlist(lapply(source_files[keep], function(f) {
    x <- readLines(f, warn = FALSE)
    i <- which(Reduce(`|`, lapply(call_patterns, function(pattern)
      grepl(pattern, x, fixed = TRUE))))
    if (length(i)) paste(rel[match(f, source_files)], i, x[i], sep = ":") else character()
  }), use.names = FALSE)
  if (length(hits)) stop("Active tarchetypes calls found: ", paste(hits, collapse = " | "), call. = FALSE)
  pkgload::load_all(root, quiet = TRUE)
  if (!identical(as.character(packageVersion("sblr")), "0.1.2"))
    stop("Unexpected installed sblr version.", call. = FALSE)
  manifests <- list(
    study02 = source(file.path(root, "studies", "02_prediction", "config.R"), local = TRUE)$value$example_data,
    study03 = source(file.path(root, "studies", "03_parameter_estimation", "config.R"), local = TRUE)$value$example_data)
  data_dirs <- c(study02 = file.path(root, "results", "local", "02_prediction", "data"),
    study03 = file.path(root, "results", "local", "03_parameter_estimation", "data"))
  for (study in names(data_dirs)) for (name in manifests[[study]]$files) {
    path <- file.path(data_dirs[[study]], name)
    if (!file.exists(path) || file.info(path)$size != manifests[[study]]$size_bytes[[name]] ||
        unname(tools::md5sum(path)) != manifests[[study]]$md5[[name]])
      stop("Cached qgdata validation failed: ", study, "/", name, call. = FALSE)
  }
  fam_n <- length(readLines(file.path(data_dirs[["study02"]], "human.fam")))
  bim_n <- length(readLines(file.path(data_dirs[["study02"]], "human.bim")))
  if (fam_n < 5000L || bim_n < 37991L) stop("Insufficient cached sample/marker metadata.", call. = FALSE)
  source(file.path(root, "studies", "02_prediction", "promotion.R"), local = TRUE)
  source(file.path(root, "studies", "03_parameter_estimation", "promotion.R"), local = TRUE)
  source(file.path(root, "studies", "04_convergence", "promotion.R"), local = TRUE)
  .study02_validate_capsule(file.path(root, "results/reference/02_prediction/st-bayesc-bayesr-one-replicate-development-v1"))
  .study03_validate_capsule(file.path(root, "results/reference/03_parameter_estimation/st-parameter-estimation-one-replicate-development-v1"))
  .study04_validate_capsule(file.path(root, "results/reference/04_convergence/st-multichain-convergence-development-v1"))
  source(file.path(root, "studies", "five_replicate_helpers.R"), local = TRUE)
  .five_replicate_recommendations()
  .five_replicate_thread_settings()
  stores <- file.path(base_dir, c("_targets_02", "_targets_03", "_targets_04"))
  for (store in stores) {
    dir.create(store, recursive = TRUE, showWarnings = FALSE)
    probe <- file.path(store, ".write_probe")
    writeLines("ok", probe); if (!file.exists(probe)) stop("Target store is not writable.")
    unlink(probe)
  }
  pid_path <- file.path(base_dir, "process_id.txt")
  if (file.exists(pid_path)) {
    old_pid <- suppressWarnings(as.integer(readLines(pid_path, n = 1L, warn = FALSE)))
    if (!is.na(old_pid) && old_pid != Sys.getpid()) {
      active <- system2("powershell", c("-NoProfile", "-Command",
        sprintf("if (Get-Process -Id %d -ErrorAction SilentlyContinue) { exit 0 } else { exit 1 }", old_pid)))
      if (identical(active, 0L)) stop("An overnight process is already running with PID ", old_pid, call. = FALSE)
    }
  }
  lines <- c("Five-replicate overnight preflight", "Result: PASSED",
    paste("Time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Repository:", normalizePath(root, winslash = "/")),
    paste("Commit:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    paste("Git dirty at orchestrator start:", length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0L),
    paste("sblr:", packageVersion("sblr")), paste("sblrbench:", packageVersion("sblrbench")),
    "tarchetypes: not required by active pipelines; no calls found and no installation attempted",
    "qgdata: cached Study 02 and shared Study 03/04 inputs match pinned sizes and checksums",
    paste("sample metadata rows:", fam_n), paste("marker metadata rows:", bim_n),
    "existing Study 02-04 capsules: valid", "target stores: writable",
    "network access: not used", "../sblr: not accessed")
  writeLines(lines, preflight_path)
  writeLines(lines, file.path(run_dir, "preflight.txt"))
  invisible(TRUE)
}

initial <- list(commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  branch = system2("git", c("branch", "--show-current"), stdout = TRUE),
  dirty = length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0L,
  started_at = format(Sys.time(), tz = "UTC", usetz = TRUE), pid = Sys.getpid(),
  validate_only = validate_only, thread_settings = as.list(Sys.getenv(c("OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))))
jsonlite::write_json(initial, file.path(run_dir, "run_metadata.json"), pretty = TRUE, auto_unbox = TRUE)

ok <- run_phase("preflight", "preflight.log", preflight())
if (validate_only) quit(status = if (ok) 0L else 1L, save = "no")

study_phase <- function(study, profile, script, store, output, promote) {
  Sys.setenv(SBLR_BENCH_PROFILE = profile, SBLR_BENCH_OUTPUT_DIR = output)
  if (study == "02_prediction") Sys.setenv(SBLR_BENCH_REPLICATES = "5") else Sys.unsetenv("SBLR_BENCH_REPLICATES")
  targets::tar_make(script = script, store = store, callr_function = NULL,
    reporter = "verbose", seconds_meta_append = 1)
  promote(output)
}

source(file.path(root, "studies", "five_replicate_promotion.R"), local = TRUE)
phase_ok <- c(preflight = ok)
phase_ok["study02"] <- run_phase("study02", "study02.log", study_phase(
  "02_prediction", "five_replicate_development", "studies/02_prediction/targets.R",
  file.path(base_dir, "_targets_02"), file.path(base_dir, "02_prediction"),
  .five_promote_study02))
phase_ok["study03"] <- run_phase("study03", "study03.log", study_phase(
  "03_parameter_estimation", "five_replicate_development",
  "studies/03_parameter_estimation/targets.R", file.path(base_dir, "_targets_03"),
  file.path(base_dir, "03_parameter_estimation"), .five_promote_study03))
phase_ok["study04"] <- run_phase("study04", "study04.log", study_phase(
  "04_convergence", "five_replicate_validation",
  "studies/04_convergence/validation_targets.R", file.path(base_dir, "_targets_04"),
  file.path(base_dir, "04_convergence"), .five_promote_study04))

phase_ok["quarto"] <- run_phase("quarto", "quarto.log", {
  quarto <- Sys.which("quarto")
  if (!nzchar(quarto)) {
    candidate <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
    if (file.exists(candidate)) quarto <- candidate else stop("Quarto executable not found.")
  }
  quarto_status <- system2(quarto, c("render"), stdout = "", stderr = "")
  if (quarto_status != 0L) stop("quarto render failed with status ", quarto_status)
})

Sys.setenv(SBLR_BENCH_PROFILE = "development")
Sys.unsetenv(c("SBLR_BENCH_REPLICATES", "SBLR_BENCH_OUTPUT_DIR"))
phase_ok["focused_tests"] <- run_phase("focused_tests", "tests.log", {
  files <- c("tests/testthat/test-five-replicate-overnight.R",
    "tests/testthat/test-prediction-helpers.R",
    "tests/testthat/test-study03-parameter-estimation.R",
    "tests/testthat/test-study04-convergence.R")
  for (f in files) testthat::test_file(f, reporter = "summary", stop_on_failure = TRUE)
})
phase_ok["full_tests"] <- run_phase("full_tests", "tests.log",
  testthat::test_local(root, reporter = "summary", stop_on_failure = TRUE))
phase_ok["check"] <- if (skip_check) {
  started <- Sys.time(); record("check", "skipped", started, "--skip-check requested", file.path(base_dir, "check.log")); TRUE
} else run_phase("check", "check.log", {
  result <- devtools::check(root, document = FALSE, error_on = "warning")
  if (length(result$errors) || length(result$warnings)) stop("Package check reported errors or warnings.")
})

final_ok <- all(phase_ok)
final_started <- Sys.time()
record("final", if (final_ok) "success" else "failed", final_started,
  if (final_ok) "all required phases completed" else
    paste("failed phases:", paste(names(phase_ok)[!phase_ok], collapse = ", ")), summary_path)
quit(status = if (final_ok) 0L else 1L, save = "no")
