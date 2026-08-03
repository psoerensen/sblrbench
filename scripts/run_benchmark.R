args <- commandArgs(trailingOnly = TRUE)

find_repository_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(path, "studies", "02_prediction", "spec.R")))
      return(path)
    parent <- dirname(path)
    if (identical(parent, path))
      stop("Could not locate the sblrbench repository root.", call. = FALSE)
    path <- parent
  }
}

root <- find_repository_root()
isolated <- file.path(root, "results", "local", "current_benchmark_refresh",
  "rlib")
if (!dir.exists(isolated))
  stop("Validated isolated benchmark library is unavailable: ", isolated,
    call. = FALSE)
.libPaths(c(normalizePath(isolated, winslash = "/"), .libPaths()))
if (!requireNamespace("devtools", quietly = TRUE))
  stop("The CLI requires devtools to load the current sblrbench source tree.",
    call. = FALSE)
devtools::load_all(root, quiet = TRUE)

options <- sblrbench:::parse_benchmark_cli_arguments(args)
spec <- read_benchmark_spec(file.path(root, "studies", "02_prediction",
  "spec.R"))
result <- run_benchmark(spec = spec, output_dir = options$output_dir,
  profile = options$profile, resume = options$resume,
  validate_only = options$validate_only)
cat("Study:", spec$study, "\n")
cat("Profile:", options$profile, "\n")
cat("Validate only:", options$validate_only, "\n")
cat("Coordinates:", nrow(result$status), "\n")
cat("Manifest:", result$paths$manifest, "\n")
