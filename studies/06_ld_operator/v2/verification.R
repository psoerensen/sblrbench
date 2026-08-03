.study06v2_quarto_executable <- function() {
  candidates <- c(Sys.which("quarto"),
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe")
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (!length(candidates)) stop("Quarto executable is unavailable.", call. = FALSE)
  candidates[1L]
}

.study06v2_render_report <- function(config) {
  source(file.path("studies", "06_ld_operator", "compatibility.R"),
    local = TRUE)
  .study06_validate_current_capsule(config$benchmark_capsule)
  report <- file.path("studies", "06_ld_operator",
    "low-rank-operator.qmd")
  log <- file.path(config$local_dir, "logs", "report.log")
  status <- system2(.study06v2_quarto_executable(), c("render", report),
    stdout = log, stderr = log)
  if (!identical(status, 0L))
    stop("Study 06 v2 Quarto report render failed; see ", log,
      call. = FALSE)
  invisible(report)
}

.study06v2_package_files <- function() {
  tracked <- system2("git", "ls-files", stdout = TRUE)
  untracked <- system2("git", c("ls-files", "--others",
    "--exclude-standard"), stdout = TRUE)
  files <- unique(c(tracked, untracked))
  files <- files[file.exists(files) & !dir.exists(files)]
  files[!grepl(paste0(
    "^(results/local|_targets|_site|fit_checkpoints|\\.git)(/|$)|",
    "(^|/)low-rank-operator(_files/|[.]html$)"),
    files)]
}

.study06v2_stage_package_source <- function(config) {
  destination <- file.path(config$local_dir, "verification", "staged_source")
  resolved <- normalizePath(file.path(config$local_dir, "verification"),
    winslash = "/", mustWork = FALSE)
  if (!startsWith(normalizePath(destination, winslash = "/", mustWork = FALSE),
      paste0(resolved, "/")))
    stop("Unsafe Study 06 v2 package staging path.", call. = FALSE)
  if (dir.exists(destination)) unlink(destination, recursive = TRUE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  files <- .study06v2_package_files()
  for (path in files) {
    target <- file.path(destination, path)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(path, target, overwrite = FALSE, copy.date = TRUE))
      stop("Could not stage package source file: ", path, call. = FALSE)
  }
  write.csv(data.frame(file = sort(files), size_bytes = file.info(files)$size),
    file.path(config$local_dir, "verification", "staged_source_files.csv"),
    row.names = FALSE)
  destination
}

.study06v2_free_disk_bytes <- function(path) {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  value <- suppressWarnings(ps::ps_disk_usage(
    normalizePath(path, winslash = "/", mustWork = TRUE))$available)
  if (length(value) == 1L && is.finite(value)) value else NA_real_
}

.study06v2_verify <- function(config) {
  output <- file.path(config$local_dir, "verification")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  before <- .study06v2_free_disk_bytes(getwd())
  focused <- testthat::test_file(
    "tests/testthat/test-study06-low-rank-operator-v2.R",
    reporter = testthat::SummaryReporter$new())
  historical <- testthat::test_file(
    "tests/testthat/test-study06-ld-operator.R",
    reporter = testthat::SummaryReporter$new())
  full <- testthat::test_local(reporter = testthat::SummaryReporter$new())
  failures <- function(x) sum(vapply(x, function(y)
    length(y$results) && any(vapply(y$results, inherits, logical(1),
      what = c("expectation_failure", "expectation_error"))), logical(1)))
  if (failures(focused) || failures(historical) || failures(full))
    stop("Study 06 v2 test validation failed.", call. = FALSE)
  .study06v2_render_report(config)
  staged <- .study06v2_stage_package_source(config)
  staged <- normalizePath(staged, winslash = "/", mustWork = TRUE)
  build_log <- file.path(output, "r_cmd_build.log")
  old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
  setwd(output)
  build <- system2(file.path(R.home("bin"), "R.exe"),
    c("CMD", "build", shQuote(staged), "--no-manual", "--no-build-vignettes"),
    stdout = basename(build_log), stderr = basename(build_log))
  if (!identical(build, 0L))
    stop("Disk-safe staged R CMD build failed; see ", build_log,
      call. = FALSE)
  package <- sort(list.files(".", pattern = "^sblrbench_.*[.]tar[.]gz$",
    full.names = TRUE), decreasing = TRUE)[1L]
  check_log <- "r_cmd_check.log"
  check <- system2(file.path(R.home("bin"), "R.exe"),
    c("CMD", "check", shQuote(normalizePath(package, winslash = "/")),
      "--no-manual", "--no-build-vignettes"), stdout = check_log,
    stderr = check_log)
  if (!identical(check, 0L))
    stop("Disk-safe staged R CMD check failed; see ", check_log,
      call. = FALSE)
  devtools_log <- "devtools_check.log"
  sink(devtools_log, split = TRUE); on.exit(sink(), add = TRUE)
  result <- devtools::check(staged, document = FALSE, manual = FALSE,
    vignettes = FALSE, cran = FALSE, quiet = FALSE)
  sink(); on.exit(NULL, add = FALSE)
  if (length(result$errors) || length(result$warnings))
    stop("Disk-safe staged devtools check returned errors or warnings.",
      call. = FALSE)
  setwd(old_wd); on.exit(NULL, add = FALSE)
  diff_check <- system2("git", c("diff", "--check"), stdout = TRUE,
    stderr = TRUE)
  if (!is.null(attr(diff_check, "status")) && attr(diff_check, "status") != 0L)
    stop("git diff --check failed.", call. = FALSE)
  after <- .study06v2_free_disk_bytes(getwd())
  write.csv(data.frame(free_bytes_before = before, free_bytes_after = after,
    staged_source_bytes = sum(file.info(.study06v2_package_files())$size),
    stringsAsFactors = FALSE), file.path(output, "disk_space.csv"),
    row.names = FALSE)
  writeLines(c("Study 06 v2 verification: PASSED",
    "focused tests: passed", "historical Study 06 tests: passed",
    "full testthat::test_local(): passed", "Quarto report: passed",
    "staged R CMD check: 0 errors, 0 warnings",
    "staged devtools::check(): 0 errors, 0 warnings",
    "git diff --check: passed"), file.path(output, "verification_summary.txt"))
  invisible(TRUE)
}
