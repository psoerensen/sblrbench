#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
validate_only <- "--validate-only" %in% args
resume <- "--resume" %in% args
phase_arg <- grep("^--phase=", args, value = TRUE)
if (length(phase_arg)) {
  phase <- sub("^--phase=", "", phase_arg[1L])
} else {
  i <- match("--phase", args)
  phase <- if (!is.na(i) && length(args) >= i + 1L)
    args[i + 1L] else "all"
}
valid_phases <- c("audit", "operator", "filter", "pilot", "convergence",
  "benchmark", "aggregate", "verify", "all")
if (!phase %in% valid_phases)
  stop("--phase must be audit, operator, filter, pilot, convergence, benchmark, aggregate, verify, or all.",
    call. = FALSE)

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
initial_head <- "0d5b7d854e655c88aac69cef59279be513f4b37d"
local_dir <- file.path(root, "results", "local",
  "study06_ld_operator")
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(local_dir, "status.csv")
summary_path <- file.path(local_dir, "summary.txt")
preflight_path <- file.path(local_dir, "preflight.txt")

empty_status <- function() data.frame(
  phase = character(), architecture = character(),
  replicate = integer(), configuration = character(),
  method = character(), operator = character(), filter = character(),
  state = character(), start_time = character(),
  finish_time = character(), elapsed_seconds = numeric(),
  output_path = character(), validation_status = character(),
  error_message = character(), stringsAsFactors = FALSE)
status <- empty_status()
if (resume && file.exists(status_path)) {
  prior <- try(utils::read.csv(status_path,
    stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(prior, "try-error") &&
      identical(names(prior), names(status))) status <- prior
}

atomic_csv <- function(x, path) {
  tmp <- tempfile(".study06-status-", dirname(path), ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  con <- file(tmp, "ab")
  flush(con)
  close(con)
  if (!file.exists(tmp) || file.info(tmp)$size < 2L)
    stop("Atomic status temporary file was not written.",
      call. = FALSE)
  backup <- paste0(path, ".previous")
  if (file.exists(backup)) unlink(backup)
  if (file.exists(path) && !file.rename(path, backup))
    stop("Could not preserve previous Study 06 status.",
      call. = FALSE)
  if (!file.rename(tmp, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    stop("Atomic Study 06 status replacement failed.",
      call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup)
}

write_summary <- function() {
  lines <- c("Study 06 single-trait LD operator benchmark",
    paste("Updated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("PID:", Sys.getpid()), paste("Requested phase:", phase),
    paste("Validation only:", validate_only),
    paste("Resume:", resume), "",
    if (nrow(status)) apply(status, 1L, function(x)
      paste(x[c("phase", "architecture", "replicate",
        "configuration", "state", "validation_status",
        "error_message")], collapse = " | "))
    else "No phases have completed.")
  writeLines(lines, summary_path, useBytes = TRUE)
}

record <- function(phase_name, state, started, output_path = "",
                   validation_status = "", error = "") {
  finished <- Sys.time()
  status <<- rbind(status, data.frame(phase = phase_name,
    architecture = "", replicate = NA_integer_,
    configuration = "", method = "", operator = "", filter = "",
    state = state,
    start_time = format(started, tz = "UTC", usetz = TRUE),
    finish_time = format(finished, tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(finished, started,
      units = "secs")), output_path = output_path,
    validation_status = validation_status,
    error_message = error, stringsAsFactors = FALSE))
  atomic_csv(status, status_path)
  write_summary()
}

with_log <- function(phase_name, log_name, expr,
                     output_path = "") {
  started <- Sys.time()
  path <- file.path(local_dir, log_name)
  con <- file(path, open = "at", encoding = "UTF-8")
  sink(con, type = "output", split = TRUE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(con)
  }, add = TRUE)
  cat("\n[", format(started, tz = "UTC", usetz = TRUE),
    "] ", phase_name, " started\n", sep = "")
  flush.console()
  ans <- tryCatch({
    force(expr)
    list(ok = TRUE, error = "")
  }, error = function(e) list(ok = FALSE,
    error = paste0(conditionMessage(e), " [",
      paste(deparse(conditionCall(e)), collapse = " "), "]")))
  cat("[", format(Sys.time(), tz = "UTC", usetz = TRUE),
    "] ", phase_name, if (ans$ok) " completed\n" else
      paste0(" failed: ", ans$error, "\n"), sep = "")
  flush.console()
  record(phase_name, if (ans$ok) "completed" else "failed",
    started, output_path, if (ans$ok) "passed" else "failed",
    ans$error)
  ans$ok
}

cfg <- source(file.path("studies", "06_ld_operator", "config.R"),
  local = TRUE)$value
source(file.path("studies", "01_finemapping",
  "setup_example_data.R"), local = TRUE)
for (f in c("blocks.R", "operators.R", "operator_validation.R",
            "simulation.R", "methods.R", "chain_extraction.R",
            "diagnostics.R", "metrics.R", "pilot.R",
            "promotion.R", "contract_smoke_test.R"))
  source(file.path("studies", "06_ld_operator", f),
    local = TRUE)

active_pid <- function(pid) {
  if (is.na(pid) || pid == Sys.getpid()) return(FALSE)
  code <- suppressWarnings(system2("powershell",
    c("-NoProfile", "-Command",
      sprintf("if (Get-Process -Id %d -ErrorAction SilentlyContinue) { exit 0 } else { exit 1 }",
        pid)), stdout = FALSE, stderr = FALSE))
  identical(code, 0L)
}

preflight <- function() {
  required <- c("sblr", "sblrbench", "targets", "posterior",
    "jsonlite", "ggplot2", "knitr", "rmarkdown", "testthat",
    "pkgload", "devtools", "qgg", "Matrix")
  available <- vapply(required, requireNamespace, logical(1),
    quietly = TRUE)
  if (!all(available))
    stop("Missing required installed packages: ",
      paste(names(available)[!available], collapse = ", "),
      call. = FALSE)
  pkgload::load_all(root, quiet = TRUE)
  head <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  branch <- system2("git", c("branch", "--show-current"),
    stdout = TRUE)
  staged <- system2("git", c("diff", "--cached", "--name-only"),
    stdout = TRUE)
  if (!identical(head, initial_head))
    stop("HEAD changed since the clean initial Study 06 gate.",
      call. = FALSE)
  if (length(staged))
    stop("Staged changes are present.", call. = FALSE)
  study05_diff <- system2("git", c("diff", "--name-only", "--",
    "studies/05_annotation_models",
    "scripts/run_study05_annotation_models.R",
    "scripts/run_study05_annotation_models.ps1",
    "docs/dev/study05_annotation_models_run.md"),
    stdout = TRUE)
  if (length(study05_diff))
    stop("Study 05 changed during Study 06 implementation: ",
      paste(study05_diff, collapse = ", "), call. = FALSE)
  if (!identical(as.character(packageVersion("sblr")), "0.1.2") ||
      !identical(packageDescription("sblr")$RemoteSha,
        "92ff3f6e7a0b1228f9f04b693d91a36d86934b0f"))
    stop("Installed sblr version or source commit changed.",
      call. = FALSE)
  exports <- getNamespaceExports("sblr")
  if (!all(c("stblr_bed", "stblr_csr",
      "stblr_block_eigen") %in% exports))
    stop("A required public sblr interface is unavailable.",
      call. = FALSE)
  filters <- eval(formals(sblr::stblr_block_eigen)$eigen_filter)
  if (!all(c("hard_truncate", "ridge_fixed",
      "ridge_lw") %in% filters))
    stop("The public Ledoit-Wolf or required block filter is unavailable.",
      call. = FALSE)
  data_dirs <- c(file.path(root,
    "results/local/02_prediction/data"),
    file.path(root,
      "results/local/03_parameter_estimation/data"))
  for (dir in data_dirs) for (name in cfg$example_data$files) {
    path <- file.path(dir, name)
    if (!file.exists(path) ||
        unname(file.info(path)$size) !=
          unname(cfg$example_data$size_bytes[name]) ||
        unname(tools::md5sum(path)) !=
          unname(cfg$example_data$md5[name]))
      stop("Cached qgdata validation failed: ", path,
        call. = FALSE)
  }
  pid_file <- file.path(local_dir, "run.pid")
  if (file.exists(pid_file)) {
    pid <- suppressWarnings(as.integer(readLines(pid_file,
      n = 1L, warn = FALSE)))
    if (active_pid(pid))
      stop("Another Study 06 process is active: ", pid,
        call. = FALSE)
  }
  stopifnot(isTRUE(run_study06_contract_smoke_test()))
  rec <- .study06_baseline_recommendations(cfg)
  if (nrow(rec) != 4L)
    stop("Study 04 recommendations are incomplete.",
      call. = FALSE)
  for (p in c(file.path(local_dir, "_targets"),
              file.path(local_dir, "operator_output"),
              file.path(local_dir, "convergence_output"),
              file.path(local_dir, "benchmark_output"))) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    probe <- file.path(p, ".write-probe")
    writeLines("ok", probe)
    if (!file.exists(probe)) stop("Unwritable path: ", p)
    unlink(probe)
  }
  manifest <- targets::tar_manifest(
    script = "studies/06_ld_operator/targets.R")
  if (!nrow(manifest) ||
      !all(c("study06_operator_files",
        "study06_convergence_files",
        "study06_benchmark_files") %in% manifest$name))
    stop("Study 06 target graph is incomplete.", call. = FALSE)
  lines <- c("Study 06 preflight", "Result: PASSED",
    paste("Time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Repository:", normalizePath(root, winslash = "/")),
    paste("Branch:", branch), paste("HEAD:", head),
    "Initial working tree gate: clean",
    "Current changes: Study 06 implementation only; unstaged",
    "Study 05 tracked files: unchanged",
    paste("sblr:", packageVersion("sblr")),
    paste("sblr source commit:",
      packageDescription("sblr")$RemoteSha),
    paste("sblrbench:", packageVersion("sblrbench")),
    paste("qgdata commit:", cfg$example_data$commit),
    "qgdata: both cached copies match pinned sizes and checksums",
    paste("public block filters:", paste(filters, collapse = ", ")),
    "unfiltered policy: ridge_fixed with eigen_eta = 0",
    "hard threshold zero: not unfiltered (native 0.01 floor)",
    "hard threshold 0.10: retained even when effective no-op",
    "Ledoit-Wolf ridge: public sensitivity route; coefficient 1 is complete diagonal shrinkage",
    "main ridge policy: ridge_fixed with frozen 1% off-diagonal shrinkage",
    "revision: the original compression/rank-reduction hard gate was superseded by the operator-agreement design",
    "tarchetypes: not required or checked; no Study 06 use",
    "thread controls: OMP=4, OMP limit=4, BLAS/MKL/vecLib=1",
    "sampler-free contract smoke test: passed",
    "target stores and output paths: writable",
    "network access: not used",
    "package installation: not attempted")
  writeLines(lines, preflight_path, useBytes = TRUE)
  cat(paste(lines, collapse = "\n"), "\n")
  invisible(TRUE)
}

run_targets <- function(names) {
  known <- targets::tar_meta(store = file.path(local_dir, "_targets"),
    fields = name)$name
  use_shortcut <- resume && all(names %in% known)
  targets::tar_make(names = tidyselect::any_of(names),
    script = "studies/06_ld_operator/targets.R",
    store = file.path(local_dir, "_targets"),
    callr_function = NULL, reporter = "verbose",
    seconds_meta_append = 1, shortcut = use_shortcut)
}

ok <- with_log("audit", "preflight.txt", preflight(),
  preflight_path)
if (!ok) quit(status = 1L, save = "no")

# Validation-only includes the full deterministic operator gate but no sampler.
if (validate_only) {
  ok <- with_log("operator_validation",
    "operator_validation.log",
    run_targets("study06_operator_files"),
    file.path(local_dir, "operator_output"))
  quit(status = if (ok) 0L else 1L, save = "no")
}

run_operator <- function() {
  run_targets("study06_operator_files")
  gate <- read.csv(file.path(local_dir, "operator_output",
    "operator_equivalence_summary.csv"))
  if (!isTRUE(gate$pass[1L]))
    stop("Study 06 deterministic operator gate did not pass.")
}

run_filter <- function() {
  run_targets(c("study06_hard_filter_pilot",
    "study06_selected_hard_filter",
    "study06_ridge_lw_operator", "study06_fixed_ridge_pilot",
    "study06_selected_ridge_filter", "study06_operator_files"))
  selected <- read.csv(file.path(local_dir, "operator_output",
    "selected_hard_filter.csv"))
  if (nrow(selected) != 1L || !isTRUE(selected$pass[1L]))
    stop("Study 06 public hard-filter route is invalid.")
  ridge <- read.csv(file.path(local_dir, "operator_output",
    "selected_ridge_filter.csv"))
  if (nrow(ridge) != 1L || !isTRUE(ridge$pass[1L]))
    stop("Study 06 fixed-ridge policy is invalid.")
}

run_pilot <- function() {
  run_targets(c("study06_operator_files", "study06_operator_pilot_file"))
  selected <- read.csv(file.path(local_dir, "operator_output",
    "selected_hard_filter.csv"))
  ridge <- read.csv(file.path(local_dir, "operator_output",
    "selected_ridge_filter.csv"))
  if (nrow(selected) != 1L || !isTRUE(selected$pass[1L]) ||
      nrow(ridge) != 1L || !isTRUE(ridge$pass[1L]))
    stop("Study 06 frozen filter policies are invalid.")
  x <- targets::tar_read(study06_operator_pilot_summary,
    store = file.path(local_dir, "_targets"))
  for (architecture in cfg$architectures) {
    z <- x[x$architecture == architecture, , drop = FALSE]
    csr <- z[z$configuration == "block_csr", , drop = FALSE]
    eig <- z[z$configuration == "block_eigen_unfiltered", , drop = FALSE]
    if (nrow(csr) != 1L || nrow(eig) != 1L ||
        abs(eig$posterior_heritability - csr$posterior_heritability) >
          cfg$pilot_gate$maximum_heritability_difference ||
        abs(eig$prediction_correlation - csr$prediction_correlation) >
          cfg$pilot_gate$maximum_prediction_correlation_difference)
      stop("One-replicate operator pilot found a large unexplained runtime-matched discrepancy for ",
        architecture, call. = FALSE)
  }
}

run_convergence <- function() {
  if (dir.exists(cfg$convergence_capsule)) {
    .study06_validate_convergence_capsule(
      cfg$convergence_capsule)
    return(invisible(cfg$convergence_capsule))
  }
  run_targets("study06_convergence_files")
  .study06_promote_convergence(
    file.path(local_dir, "operator_output"),
    file.path(local_dir, "convergence_output"), cfg)
  .study06_validate_convergence_capsule(
    cfg$convergence_capsule)
}

run_benchmark <- function() {
  .study06_validate_convergence_capsule(
    cfg$convergence_capsule)
  run_targets("study06_benchmark_files")
  status <- read.csv(file.path(local_dir, "benchmark_output",
    "fit_status.csv"))
  if (nrow(status) != 60L || any(status$status != "ok"))
    stop("Study 06 main 60-fit grid is incomplete.")
}

run_aggregate <- function() {
  run_benchmark()
  if (!dir.exists(cfg$benchmark_capsule))
    .study06_promote_benchmark(
      file.path(local_dir, "benchmark_output"), cfg)
  .study06_validate_benchmark_capsule(cfg$benchmark_capsule)
}

find_quarto <- function() {
  candidates <- c(Sys.which("quarto"),
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe")
  candidates <- candidates[nzchar(candidates) &
    file.exists(candidates)]
  if (!length(candidates)) stop("Quarto executable not found.")
  candidates[1L]
}

run_verify <- function() {
  .study06_validate_convergence_capsule(cfg$convergence_capsule)
  .study06_validate_benchmark_capsule(cfg$benchmark_capsule)
  code <- system2(find_quarto(), "render")
  if (!identical(code, 0L)) stop("quarto render failed.")
  testthat::test_file(
    "tests/testthat/test-study06-ld-operator.R",
    stop_on_failure = TRUE, stop_on_warning = TRUE)
  testthat::test_local(stop_on_failure = TRUE)
  # pkgbuild first copies the whole working directory before R applies
  # .Rbuildignore. The ignored benchmark stores are intentionally large, so
  # stage only the package paths that R CMD build would retain. This preserves
  # every cache while checking the exact working-tree R code and tests.
  tracked <- system2("git", c("ls-files", "--cached", "--others",
    "--exclude-standard"), stdout = TRUE)
  package_roots <- c("R", "man", "tests")
  package_files <- c(".Rbuildignore", ".gitattributes", ".gitignore",
    "DESCRIPTION", "NAMESPACE", "README.md")
  keep <- tracked %in% package_files |
    sub("[/\\\\].*$", "", tracked) %in% package_roots
  tracked <- sort(unique(tracked[keep & file.exists(tracked)]))
  check_source <- file.path(local_dir, "package_check_source",
    format(Sys.time(), "%Y%m%dT%H%M%S"))
  dir.create(check_source, recursive = TRUE, showWarnings = FALSE)
  for (path in tracked) {
    destination <- file.path(check_source, path)
    dir.create(dirname(destination), recursive = TRUE,
      showWarnings = FALSE)
    if (!file.copy(path, destination, overwrite = FALSE,
        copy.date = TRUE))
      stop("Could not stage package-check source: ", path)
  }
  check <- devtools::check(pkg = check_source, document = FALSE,
    error_on = "warning", quiet = FALSE)
  if (length(check$errors) || length(check$warnings))
    stop("devtools::check reported errors or warnings.")
  code <- system2("git", c("diff", "--check"))
  if (!identical(code, 0L)) stop("git diff --check failed.")
  invisible(TRUE)
}

requested <- if (phase == "all")
  c("pilot", "convergence", "benchmark",
    "aggregate", "verify") else phase
phase_ok <- logical()
for (item in requested) {
  if (item == "audit") next
  fun <- switch(item, operator = run_operator,
    filter = run_filter, pilot = run_pilot, convergence = run_convergence,
    benchmark = run_benchmark, aggregate = run_aggregate,
    verify = run_verify)
  log <- switch(item,
    operator = "operator_validation.log",
    filter = "pilot.log", pilot = "pilot.log", convergence = "convergence.log",
    benchmark = "benchmark.log", aggregate = "aggregate.log",
    verify = "render_and_tests.log")
  output <- switch(item,
    operator = file.path(local_dir, "operator_output"),
    filter = file.path(local_dir, "operator_output"),
    pilot = file.path(local_dir, "operator_output"),
    convergence = cfg$convergence_capsule,
    benchmark = file.path(local_dir, "benchmark_output"),
    aggregate = cfg$benchmark_capsule,
    verify = file.path(root, "_site"))
  phase_ok[item] <- with_log(item, log, fun(), output)
  if (!phase_ok[item]) quit(status = 1L, save = "no")
}

quit(status = 0L, save = "no")
