#!/usr/bin/env Rscript
options(warn = 1)
Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
value_arg <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, paste0(prefix, "="))]
  if (length(hit)) sub(paste0("^", prefix, "="), "", hit[[1L]]) else default
}
phase <- value_arg("--phase", "all")
if ("--phase" %in% args) {
  i <- match("--phase", args)
  if (i == length(args)) stop("--phase requires a value.")
  phase <- args[[i + 1L]]
}
validate_only <- "--validate-only" %in% args
resume <- "--resume" %in% args
valid_phases <- c("audit", "contract", "runtime", "convergence",
  "benchmark", "stress", "aggregate", "verify", "all")
if (!phase %in% valid_phases) stop("Unknown Study 07 phase: ", phase)

locate_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(here, "DESCRIPTION")) &&
        file.exists(file.path(here, "studies", "07_mtblr_validation",
          "config.R"))) return(here)
    parent <- dirname(here)
    if (identical(parent, here)) stop("Cannot locate sblrbench root.")
    here <- parent
  }
}
root <- locate_root(); setwd(root)
cfg <- source("studies/07_mtblr_validation/config.R", local = TRUE)$value
source("studies/07_mtblr_validation/state_contract.R", local = TRUE)
source("studies/07_mtblr_validation/promotion.R", local = TRUE)
source("studies/07_mtblr_validation/contract_smoke_test.R", local = TRUE)
local_dir <- cfg$local_dir
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
status_path <- file.path(local_dir, "status.csv")
summary_path <- file.path(local_dir, "summary.txt")
preflight_path <- file.path(local_dir, "preflight.txt")
store <- file.path(local_dir, "_targets")

atomic_csv <- function(x, path) {
  tmp <- tempfile(".status-", dirname(path), ".csv")
  write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) { unlink(tmp); stop("Atomic status update failed.") }
}
read_status <- function() if (file.exists(status_path))
  read.csv(status_path, stringsAsFactors = FALSE) else data.frame()
record <- function(phase, state, started, output = "", validation = "",
                   error = "") {
  row <- data.frame(phase = phase, architecture = "", replicate = NA_integer_,
    implementation = "", marker_count = NA_integer_, chain_count = NA_integer_,
    state = state, start_time = format(started, tz = "UTC", usetz = TRUE),
    finish_time = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    memory_estimate = NA_real_, output_path = output,
    validation_status = validation, error_message = error,
    stringsAsFactors = FALSE)
  old <- read_status(); atomic_csv(if (nrow(old)) rbind(old, row) else row,
    status_path)
  cat(sprintf("[%s] %s %s\n", row$finish_time, phase, state),
    file = summary_path, append = TRUE)
  flush.console()
}

active_pid <- function(pid) {
  if (!is.finite(pid) || pid <= 0 || pid == Sys.getpid()) return(FALSE)
  identical(system2("powershell", c("-NoProfile", "-Command",
    sprintf("if (Get-Process -Id %d -ErrorAction SilentlyContinue) { exit 0 } else { exit 1 }", pid)),
    stdout = FALSE, stderr = FALSE), 0L)
}

preflight <- function() {
  started <- Sys.time()
  head <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
  branch <- trimws(system2("git", c("branch", "--show-current"), stdout = TRUE))
  status <- system2("git", c("status", "--short"), stdout = TRUE)
  allowed <- c("DESCRIPTION", "_quarto.yml", "index.qmd", "studies/index.qmd",
    "docs/dev/sblr_implementation_map.md", "docs/dev/study07_mtblr_validation_run.md",
    "scripts/run_study07_mtblr_validation.R",
    "scripts/run_study07_mtblr_validation.ps1",
    "studies/07_mtblr_validation", "tests/testthat/test-study07-mtblr-validation.R",
    "results/reference/07_mtblr_validation")
  changed <- trimws(sub("^..", "", status))
  unrelated <- changed[!vapply(changed, function(x)
    any(startsWith(x, allowed)), logical(1))]
  if (length(unrelated)) stop("Unrelated working-tree changes: ",
    paste(unrelated, collapse = ", "))
  required_packages <- c("sblr", "sblrbench", "targets", "posterior",
    "jsonlite", "Matrix", "qgg", "testthat", "devtools")
  missing <- required_packages[!vapply(required_packages,
    requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ",
    paste(missing, collapse = ", "))
  if (!identical(as.character(packageVersion("sblr")), "0.1.2") ||
      !identical(packageDescription("sblr")$RemoteSha,
        "92ff3f6e7a0b1228f9f04b693d91a36d86934b0f"))
    stop("Installed sblr provenance changed.")
  if (!all(c("mtblr_bed", "mtblr_csr", "mtblr_block_eigen") %in%
      getNamespaceExports("sblr"))) stop("Required MTBLR interface missing.")
  data_dirs <- c("results/local/02_prediction/data",
    "results/local/03_parameter_estimation/data")
  for (dir in data_dirs) for (name in cfg$example_data$files) {
    path <- file.path(dir, name)
    if (!file.exists(path) || unname(file.info(path)$size) !=
        unname(cfg$example_data$size_bytes[name]) ||
        unname(tools::md5sum(path)) != unname(cfg$example_data$md5[name]))
      stop("Cached qgdata validation failed: ", path)
  }
  pid_file <- file.path(local_dir, "run.pid")
  if (file.exists(pid_file)) {
    pid <- suppressWarnings(as.integer(readLines(pid_file, n = 1L)))
    if (active_pid(pid)) stop("Another Study 07 process is active: ", pid)
  }
  stopifnot(isTRUE(run_study07_contract_smoke_test()))
  manifest <- targets::tar_manifest(script =
    "studies/07_mtblr_validation/targets.R")
  expected <- c("study07_contract_output_files",
    "study07_convergence_output_files", "study07_main_output_files")
  if (!all(expected %in% manifest$name)) stop("Study 07 target graph incomplete.")
  lines <- c("Study 07 preflight", "Result: PASSED",
    paste("Time:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Repository:", normalizePath(root, winslash = "/")),
    paste("Branch:", branch), paste("HEAD:", head),
    paste("sblr:", packageVersion("sblr")),
    paste("sblr source commit:", packageDescription("sblr")$RemoteSha),
    paste("sblrbench:", packageVersion("sblrbench")),
    paste("qgdata commit:", cfg$example_data$commit),
    "qgdata: both cached copies match pinned sizes and checksums",
    "interfaces: mtblr_bed, mtblr_csr, mtblr_block_eigen",
    "joint states: 0_0, 1_0, 0_1, 1_1",
    "sampler-free contracts: passed", "active duplicate process: absent",
    "thread controls: OMP=4, OMP limit=4, BLAS/MKL/vecLib=1",
    "network access: not used", "package installation: not attempted")
  writeLines(lines, preflight_path); record("audit", "completed", started,
    preflight_path, "passed"); invisible(TRUE)
}

run_target <- function(name, log_name) {
  log <- file.path(local_dir, log_name)
  con <- file(log, open = "at"); sink(con, type = "output"); sink(con, type = "message")
  on.exit({sink(type = "message"); sink(type = "output"); close(con)}, add = TRUE)
  targets::tar_make(names = name, script =
    "studies/07_mtblr_validation/targets.R", store = store,
    reporter = "timestamp")
}
run_phase <- function(id, target, log, after = NULL, output = "") {
  started <- Sys.time()
  tryCatch({run_target(target, log); if (!is.null(after)) after();
    record(id, "completed", started, output, "passed")}, error = function(e) {
      record(id, "failed", started, output, "failed", conditionMessage(e)); stop(e)
  })
}

preflight()
if (validate_only) quit(save = "no", status = 0L)
paused_phases <- c("runtime", "convergence", "benchmark", "stress",
  "aggregate", "all")
if (phase %in% paused_phases) {
  message <- paste(
    "Study 07 MT block-eigen execution is paused.",
    "The current mtblr_block_eigen() backend is reconstructed dense.",
    "Resume this phase only after the retained low-rank MT operator",
    "has been implemented and validated in sblr."
  )
  record(phase, "paused", Sys.time(), "", "paused", message)
  stop(message, call. = FALSE)
}
sequence <- if (phase == "all") c("contract", "runtime", "convergence",
  "benchmark", "aggregate", "verify") else phase
for (p in sequence) switch(p,
  audit = preflight(),
  contract = run_phase("contract", "study07_tiny_fit_status", "contract.log"),
  runtime = run_phase("runtime", "study07_contract_output_files", "runtime.log",
    function() .study07_promote_contract(cfg), cfg$contract_capsule),
  convergence = run_phase("convergence", "study07_convergence_output_files",
    "convergence.log", function() .study07_promote_convergence(cfg),
    cfg$convergence_capsule),
  benchmark = run_phase("benchmark", "study07_main_output_files",
    "benchmark.log"),
  aggregate = run_phase("aggregate", "study07_main_output_files",
    "aggregate.log", function() .study07_promote_benchmark(cfg),
    cfg$benchmark_capsule),
  stress = record("stress", "skipped", Sys.time(), "",
    "optional_not_requested"),
  verify = {
    started <- Sys.time(); log <- file.path(local_dir, "render_and_tests.log")
    status <- tryCatch({
      q <- Sys.which("quarto"); if (!nzchar(q))
        q <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
      if (system2(q, "render", stdout = log, stderr = log) != 0L)
        stop("Quarto render failed.")
      testthat::test_file("tests/testthat/test-study07-mtblr-validation.R",
        stop_on_failure = TRUE)
      testthat::test_local(stop_on_failure = TRUE)
      result <- devtools::check(error_on = "warning")
      if (length(result$errors) || length(result$warnings))
        stop("Package check reported errors or warnings.")
      if (system2("git", c("diff", "--check")) != 0L)
        stop("git diff --check failed.")
      "passed"
    }, error = function(e) {record("verify", "failed", started, "_site",
      "failed", conditionMessage(e)); stop(e)})
    record("verify", "completed", started, "_site", status)
  })
cat("Study 07 requested phases completed.\n")
