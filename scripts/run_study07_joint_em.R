args <- commandArgs(trailingOnly = TRUE)
gate_arg <- grep("^--gate=", args, value = TRUE)
if (!length(gate_arg)) {
  index <- match("--gate", args)
  gate <- if (!is.na(index) && index < length(args)) args[[index + 1L]] else
    "validate"
} else gate <- sub("^--gate=", "", gate_arg[[1L]])
if (!gate %in% c("validate", "baseline", "joint", "em", "selection",
    "summarize", "all"))
  stop("Unknown Study 07 gate: ", gate, call. = FALSE)

pkgload::load_all(".", quiet = TRUE)
study_dir <- file.path("studies", "07_joint_em_sbayesrc")
spec <- source(file.path(study_dir, "spec.R"), local = TRUE)$value
source(file.path(study_dir, "data.R"), local = FALSE)
source(file.path(study_dir, "methods.R"), local = FALSE)
source(file.path(study_dir, "extraction.R"), local = FALSE)
source(file.path(study_dir, "figures.R"), local = FALSE)
inputs <- study07_load_inputs(spec)
validation <- study07_validate_design(spec, inputs)
cat("Gate 0 PASS: design/input frozen\n")

run <- function(method, start = "primary") {
  cat("Running ", method, " [", start, "]\n", sep = "")
  value <- study07_run_fit(method, spec, inputs, validation, start = start)
  cat("Completed ", method, " [", start, "] in ",
    format(value$elapsed_seconds, digits = 5), " seconds; reused=",
    value$reused, "\n", sep = "")
  invisible(value)
}

if (gate %in% c("baseline", "all")) run("SBayesR")
if (gate %in% c("joint", "all")) run("SBayesRC")
if (gate %in% c("em", "all"))
  for (start in spec$em$starts) run("SBayesRC-EM", start)
if (gate %in% c("selection", "all")) {
  run("SBayesRC-S")
  for (start in spec$em$starts) run("SBayesRC-S-EM", start)
}
if (gate %in% c("summarize", "all")) {
  output <- study07_summarize(spec, inputs)
  output$figures <- study07_make_figures(output, inputs, spec)
  cat("Wrote Study 07 compact tables:\n",
    paste(unname(output$files), collapse = "\n"), "\n")
  cat("Wrote Study 07 figures:\n",
    paste(unname(output$figures), collapse = "\n"), "\n")
}
