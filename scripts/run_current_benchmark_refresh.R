args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[[1L]] == length(args)) stop(flag, " requires a value.", call. = FALSE)
  args[[hit[[1L]] + 1L]]
}

phase <- arg_value("--phase", "audit")
resume <- "--resume" %in% args
allowed <- c("preflight", "audit", "study04-selection", "study04-validation",
  "study03", "study02", "study01", "study05", "study06-compatibility",
  "cleanup", "verify", "all")
if (!phase %in% allowed) stop("Unsupported refresh phase: ", phase, call. = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION")) ||
    !identical(read.dcf(file.path(root, "DESCRIPTION"), fields = "Package")[[1L]],
      "sblrbench")) stop("Run this script from the sblrbench repository root.", call. = FALSE)

refresh_sha <- "02e8c74baa906e83c4a08d42a9cc6339b4e81072"
refresh_root <- file.path(root, "results", "local", "current_benchmark_refresh")
rlib <- file.path(refresh_root, "rlib")
recommendations <- file.path(root, "results", "reference", "04_convergence",
  "current-selection", "method_recommendations.csv")
directories <- file.path(refresh_root, c("audit", "study04", "study03", "study02",
  "study01", "study05", "study06_compatibility", "rendering", "package_checks",
  "logs", "targets"))
invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
.libPaths(c(rlib, .libPaths()))

atomic_write <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  con <- file(temporary, open = "wb")
  tryCatch({ writeLines(lines, con, useBytes = TRUE); flush(con) },
    finally = close(con))
  if (file.exists(path) && !file.remove(path))
    stop("Cannot replace status file: ", path, call. = FALSE)
  if (!file.rename(temporary, path)) stop("Atomic rename failed: ", path, call. = FALSE)
  invisible(path)
}

git <- function(..., directory = root) {
  out <- system2("git", c("-C", shQuote(directory), ...), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status") %||% 0L
  if (status) stop(paste(out, collapse = "\n"), call. = FALSE)
  trimws(out)
}
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

verify_package <- function() {
  if (!dir.exists(rlib)) stop("The isolated refresh library is missing.", call. = FALSE)
  suppressPackageStartupMessages(library("sblr", lib.loc = rlib, character.only = TRUE))
  desc <- utils::packageDescription("sblr")
  observed <- desc$RemoteSha %||% desc$GithubSHA1 %||% ""
  if (!identical(as.character(utils::packageVersion("sblr")), "0.2.0") ||
      !identical(observed, refresh_sha) ||
      !startsWith(normalizePath(find.package("sblr"), winslash = "/"),
        normalizePath(rlib, winslash = "/")))
    stop("The active sblr package does not match the frozen refresh provenance.", call. = FALSE)
  required <- c("stblr_bed", "stblr_csr", "stblr_csr_annot", "stblr_block_eigen")
  if (!all(vapply(required, exists, logical(1), envir = asNamespace("sblr"),
      inherits = FALSE))) stop("A required sblr API is unavailable.", call. = FALSE)
  list(version = "0.2.0", sha = observed, library = find.package("sblr"))
}

preflight <- function() {
  if (!identical(git("branch", "--show-current"), "master"))
    stop("The refresh must run on master.", call. = FALSE)
  sibling <- normalizePath(file.path(root, "..", "sblr"), winslash = "/", mustWork = TRUE)
  if (!identical(git("branch", "--show-current", directory = sibling), "master") ||
      length(git("status", "--porcelain", directory = sibling)) ||
      !identical(git("rev-parse", "HEAD", directory = sibling), refresh_sha) ||
      !identical(git("rev-parse", "origin/master", directory = sibling), refresh_sha))
    stop("The sibling sblr repository is not the frozen clean synchronized master.", call. = FALSE)
  pkg <- verify_package()
  record <- c(
    paste("recorded_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("sblrbench_head:", git("rev-parse", "HEAD")),
    paste("sblrbench_branch:", git("branch", "--show-current")),
    paste("sblr_refresh_sha:", refresh_sha),
    paste("sblr_subject:", git("show", "-s", "--format=%s", refresh_sha,
      directory = sibling)),
    paste("sblr_version:", pkg$version), paste("sblr_library:", pkg$library),
    paste("R:", R.version.string), paste("platform:", R.version$platform),
    paste("BLAS:", extSoftVersion()[["BLAS"]] %||% "unavailable"))
  atomic_write(record, file.path(refresh_root, "audit", "preflight.txt"))
  invisible(pkg)
}

classify_inventory <- function(path) {
  if (grepl("^results/reference/", path)) {
    if (grepl("/(current|current-selection|current-validation|current-stop|current-convergence)(/|$)", path))
      return("required current capsule")
    return("superseded and removable")
  }
  if (grepl("development|overnight_run|study0[5-7]_.*_run", path, ignore.case = TRUE))
    return("stale prose")
  if (grepl("(^|/)tests/", path)) return("required source code")
  if (grepl("(^|/)studies/|(^|/)scripts/", path)) return("required source code")
  if (grepl("\\.qmd$|_quarto\\.yml$", path)) return("stale prose")
  "unresolved dependency"
}

write_inventory <- function() {
  tracked <- git("ls-files")
  keep <- grepl("^(studies/|results/reference/|docs/dev/|scripts/|tests/|_quarto\\.yml$|index\\.qmd$|framework\\.qmd$|reproducibility\\.qmd$)", tracked)
  paths <- tracked[keep]
  info <- file.info(file.path(root, paths))
  inventory <- data.frame(path = paths, size_bytes = as.numeric(info$size),
    classification = vapply(paths, classify_inventory, character(1)),
    refresh_action = ifelse(grepl("^results/reference/", paths) &
      !grepl("/(current|current-selection|current-validation|current-stop|current-convergence)(/|$)", paths),
      "remove after replacement validates", "review and retain or consolidate"),
    stringsAsFactors = FALSE)
  inventory <- inventory[order(inventory$path), ]
  path <- file.path(refresh_root, "audit", "benchmark_inventory.csv")
  utils::write.csv(inventory, path, row.names = FALSE)
  path
}

write_prior_impact <- function() {
  rows <- data.frame(
    study = c("01", "01", "02", "02", "03", "03", "04", "04", "05", "05",
      "06", "06"),
    model = rep(c("BayesC/SBayesC", "BayesR/SBayesR"), 6L),
    old_policy = c(rep(c("m * pi_nonnull", "m * pi_nonnull"), 4L),
      "m * marker-averaged annotation nonnull probability",
      "m * marker-averaged annotation nonnull probability",
      "resolved global expected multiplier", "resolved global expected multiplier"),
    current_policy = rep("resolved expected genetic variance calibration", 12L),
    B_formula = rep("B = vy * h2 / weight_initial", 12L),
    ssb_prior_formula = rep(
      "ssb_prior = ((nub - 2) / nub) * vy * h2 / weight_prior", 12L),
    initial_component_weight = rep(c("pi_nonnull", "sum(pi * mixture_multiplier)"), 6L),
    prior_mean_component_weight = rep(c("pi_prior_mean",
      "sum(normalized_alpha * mixture_multiplier)"), 6L),
    mixture_multiplier_weight = rep(c("1", "c(0, 0.01, 0.1, 1)"), 6L),
    marker_or_annotation_weight = c(rep("none", 8L), rep("annotation dependent", 2L),
      rep("none", 2L)),
    numerical_change = c(rep(c(FALSE, TRUE), 4L), rep(TRUE, 2L), rep(FALSE, 2L)),
    decision = c(rep("fresh complete rerun", 8L), rep("fresh convergence pilot", 2L),
      rep("compatibility validation; retain fits only if sampler hashes and resolved priors match", 2L)),
    stringsAsFactors = FALSE)
  path <- file.path(refresh_root, "audit", "prior_impact.csv")
  utils::write.csv(rows, path, row.names = FALSE)
  path
}

run_targets <- function(study, profile, output_dir, store_name, script = "_targets.R") {
  preflight()
  # The benchmark package under development is loaded from this committed source tree so
  # targets that call sblrbench:: APIs see the matching namespace. The scientific backend
  # remains the ordinary installed sblr package from the frozen isolated library.
  pkgload::load_all(root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  verify_package()
  vars <- c(SBLR_BENCH_STUDY = study, SBLR_BENCH_PROFILE = profile,
    SBLR_BENCH_OUTPUT_DIR = output_dir,
    SBLR_BENCH_RECOMMENDATIONS = recommendations,
    SBLR_BENCH_REFRESH_ROOT = refresh_root,
    OMP_NUM_THREADS = "4", OMP_THREAD_LIMIT = "4", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  old <- Sys.getenv(names(vars), unset = NA_character_)
  do.call(Sys.setenv, as.list(vars))
  on.exit({
    restore <- !is.na(old); if (any(restore)) do.call(Sys.setenv, as.list(old[restore]))
    if (any(!restore)) Sys.unsetenv(names(old)[!restore])
  }, add = TRUE)
  store <- file.path(refresh_root, "targets", store_name)
  targets::tar_make(script = script, store = store, callr_function = NULL)
}

run_phase <- function(x) switch(x,
  preflight = preflight(),
  audit = { preflight(); write_inventory(); write_prior_impact() },
  `study04-selection` = run_targets("04_convergence", "development",
    file.path(refresh_root, "study04", "selection"), "study04-selection"),
  `study04-validation` = run_targets("04_convergence", "five_replicate_validation",
    file.path(refresh_root, "study04", "validation"), "study04-validation",
    script = file.path("studies", "04_convergence", "validation_targets.R")),
  `study03` = run_targets("03_parameter_estimation", "five_replicate_development",
    file.path(refresh_root, "study03"), "study03"),
  `study02` = run_targets("02_prediction", "five_replicate_development",
    file.path(refresh_root, "study02"), "study02"),
  `study01` = run_targets("01_finemapping", "pilot",
    file.path(refresh_root, "study01"), "study01"),
  stop("Phase is declared but not implemented yet: ", x, call. = FALSE))

if (phase == "all") {
  for (x in c("audit", "study04-selection", "study04-validation", "study03",
      "study02", "study01", "study05", "study06-compatibility", "cleanup", "verify"))
    run_phase(x)
} else run_phase(phase)

atomic_write(c(paste("phase:", phase), "status: complete",
  paste("finished_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE))),
  file.path(refresh_root, "logs", paste0(phase, ".complete")))
