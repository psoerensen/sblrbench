# Compatibility loader for reports in their pre-migration locations.
# Authoritative implementation: R/benchmark-reporting.R
# Remove after all report callers load shared package helpers directly.
.sblrbench_root <- getwd()
while (!file.exists(file.path(.sblrbench_root, "DESCRIPTION"))) {
  .sblrbench_parent <- dirname(.sblrbench_root)
  if (identical(.sblrbench_parent, .sblrbench_root)) stop("Cannot locate sblrbench root.")
  .sblrbench_root <- .sblrbench_parent
}
source(file.path(.sblrbench_root, "R", "benchmark-reporting.R"), local = TRUE)
