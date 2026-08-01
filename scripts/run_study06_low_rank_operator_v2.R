#!/usr/bin/env Rscript
options(warn = 1)
Sys.setenv(OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) return(args[[hit + 1L]])
  inline <- args[startsWith(args, paste0(flag, "="))]
  if (length(inline)) sub(paste0("^", flag, "="), "", inline[[1L]]) else default
}
phase <- arg_value("--phase", "all")
resume <- "--resume" %in% args
validate_only <- "--validate-only" %in% args
valid <- c("audit", "deterministic", "operator-pilot", "convergence",
  "benchmark", "aggregate", "report", "verify", "all")
if (!phase %in% valid) stop("Unknown Study 06 v2 phase: ", phase)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION")) ||
    !file.exists(file.path(root, "studies", "06_ld_operator", "v2", "config.R")))
  stop("Run Study 06 v2 from the sblrbench repository root.", call. = FALSE)
config <- source("studies/06_ld_operator/v2/config.R", local = TRUE)$value
pkgload::load_all(root, recompile = FALSE, quiet = TRUE,
  export_all = FALSE, helpers = FALSE)
for (file in c("source_loader.R", "guard.R", "design_crosswalk.R",
    "operator_validation.R", "deterministic_run.R", "methods.R",
    "runtime_data.R", "phases.R", "provenance.R", "promotion.R",
    "verification.R"))
  source(file.path("studies", "06_ld_operator", "v2", file), local = TRUE)

local_dir <- config$local_dir
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
for (directory in c("preflight", "historical_v1_inventory", "deterministic",
    "operator_pilot", "convergence", "fit_checkpoints", "benchmark",
    "aggregate", "verification", "logs", "_targets"))
  dir.create(file.path(local_dir, directory), recursive = TRUE,
    showWarnings = FALSE)
status_path <- file.path(local_dir, "status.csv")
summary_path <- file.path(local_dir, "summary.txt")

atomic_status <- function(row) {
  old <- if (file.exists(status_path)) read.csv(status_path,
    stringsAsFactors = FALSE) else data.frame()
  out <- if (nrow(old)) rbind(old, row) else row
  tmp <- tempfile(".status-", dirname(status_path), ".csv")
  write.csv(out, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, status_path)) {
    unlink(tmp); stop("Atomic Study 06 v2 status write failed.")
  }
}
record <- function(id, state, started, output = "", error = "") {
  row <- data.frame(phase = id, state = state,
    start_time = format(started, tz = "UTC", usetz = TRUE),
    finish_time = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    output_path = output, error_message = error, stringsAsFactors = FALSE)
  atomic_status(row)
  cat(sprintf("[%s] %s %s\n", row$finish_time, id, state),
    file = summary_path, append = TRUE)
}

run_audit <- function() {
  started <- Sys.time()
  tryCatch({
    bench_root <- trimws(system2("git", c("rev-parse", "--show-toplevel"),
      stdout = TRUE))
    sibling <- normalizePath(file.path(root, "..", "sblr"), winslash = "/",
      mustWork = TRUE)
    sibling_head <- trimws(system2("git", c("-C", sibling, "rev-parse", "HEAD"),
      stdout = TRUE))
    sibling_status <- system2("git", c("-C", sibling, "status", "--short"),
      stdout = TRUE)
    if (!identical(sibling_head, config$required_sblr_sha) ||
        length(sibling_status))
      stop("Pinned sibling sblr revision is unavailable or dirty.")
    source_info <- .study06v2_load_pinned_sblr(config, recompile = FALSE)
    .study06v2_validate_grid(config)
    for (capsule in config$historical_capsules)
      if (!dir.exists(capsule)) stop("Historical Study 06 v1 capsule missing: ", capsule)
    for (data_dir in c("results/local/02_prediction/data",
        "results/local/03_parameter_estimation/data"))
      for (name in config$example_data$files) {
        path <- file.path(data_dir, name)
        if (!file.exists(path) || file.info(path)$size !=
            unname(config$example_data$size_bytes[name]) ||
            unname(tools::md5sum(path)) != unname(config$example_data$md5[name]))
          stop("Pinned qgdata validation failed: ", path)
      }
    lines <- c("Study 06 v2 audit: PASSED", paste("repository:", bench_root),
      paste("branch:", trimws(system2("git", c("branch", "--show-current"),
        stdout = TRUE))), paste("HEAD:", trimws(system2("git",
        c("rev-parse", "HEAD"), stdout = TRUE))),
      paste("sblr_source_sha:", source_info$sha),
      paste("sblr_version:", source_info$version),
      paste("sblr_source_path:", source_info$path),
      "operator_contract: block_low_rank_v1",
      "representation: low_rank", "legacy_sampler_allowed: FALSE",
      "qgdata: both pinned cache copies valid",
      "historical_v1_capsules: present and read-only")
    path <- file.path(local_dir, "preflight", "audit.txt")
    writeLines(lines, path)
    record("audit", "completed", started, path)
    TRUE
  }, error = function(e) {
    record("audit", "failed", started, "", conditionMessage(e)); stop(e)
  })
}

run_deterministic <- function() {
  started <- Sys.time()
  tryCatch({
    .study06v2_load_pinned_sblr(config, recompile = FALSE)
    .study06v2_run_deterministic(config)
    output <- file.path(local_dir, "deterministic",
      "deterministic_identity_summary.csv")
    record("deterministic", "completed", started, output)
  }, error = function(e) {
    record("deterministic", "failed", started, "", conditionMessage(e)); stop(e)
  })
}

run_phase <- function(id, fun, output) {
  started <- Sys.time()
  tryCatch({
    .study06v2_load_pinned_sblr(config, recompile = FALSE)
    fun(config)
    record(id, "completed", started, output)
    TRUE
  }, error = function(e) {
    record(id, "failed", started, "", conditionMessage(e)); stop(e)
  })
}

already_complete <- function(id, output) {
  resume && file.exists(output) && nrow(read.csv(status_path,
    stringsAsFactors = FALSE)) && any(with(read.csv(status_path,
      stringsAsFactors = FALSE), phase == id & state == "completed"))
}

dispatch <- function(id) switch(id,
  deterministic = {
    output <- file.path(local_dir, "deterministic",
      "deterministic_identity_summary.csv")
    if (!already_complete(id, output)) run_deterministic()
  },
  `operator-pilot` = {
    output <- file.path(local_dir, "operator_pilot", "pilot_gate.csv")
    if (!already_complete(id, output))
      run_phase(id, .study06v2_pilot, output)
  },
  convergence = {
    output <- file.path(local_dir, "convergence",
      "method_recommendations.csv")
    if (!already_complete(id, output))
      run_phase(id, .study06v2_convergence, output)
  },
  benchmark = {
    output <- file.path(local_dir, "fit_checkpoints", "benchmark")
    run_phase(id, .study06v2_benchmark, output)
  },
  aggregate = {
    output <- file.path(local_dir, "aggregate", "fit_status.csv")
    if (!already_complete(id, output)) {
      run_phase(id, .study06v2_aggregate, output)
      .study06v2_promote("convergence", config)
      readiness <- .study06v2_readiness_decision(config)
      .study06v2_promote("benchmark", config, readiness)
    }
  },
  report = {
    output <- file.path(local_dir, "logs", "report.log")
    if (!already_complete(id, output))
      run_phase(id, .study06v2_render_report, output)
  },
  verify = {
    output <- file.path(local_dir, "verification", "verification_summary.txt")
    if (!already_complete(id, output))
      run_phase(id, .study06v2_verify, output)
  },
  audit = NULL)

run_audit()
if (validate_only) quit(save = "no", status = 0L)
sequence <- if (phase == "all") valid[valid != "all"] else phase
for (item in sequence) dispatch(item)
