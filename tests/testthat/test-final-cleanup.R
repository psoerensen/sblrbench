root <- testthat::test_path("..", "..")

test_that("completed studies expose clean source contracts", {
  completed <- sprintf("%02d_%s", 1:5, c("finemapping", "prediction",
    "parameter_estimation", "convergence", "ld_operator"))
  for (study in completed) {
    path <- file.path(root, "studies", study)
    expect_true(all(file.exists(file.path(path,
      c("spec.R", "analysis.R", "report.qmd")))))
    expect_false(file.exists(file.path(path, "targets.R")))
  }
})

test_that("obsolete smoke and compatibility infrastructure is absent", {
  obsolete <- c(file.path("studies", "00_contract_smoke"),
    file.path("studies", "reporting_helpers.R"),
    file.path("studies", "five_replicate_helpers.R"),
    file.path("scripts", "contract_smoke_test.R"),
    file.path("scripts", "prediction_contract_smoke_test.R"),
    file.path("scripts", "run_current_benchmark_refresh.R"),
    file.path("scripts", "run_study06_annotation_models.R"),
    file.path("scripts", "run_study08_mtblr_validation.R"))
  expect_false(any(file.exists(file.path(root, obsolete))))
})

test_that("root targets dispatch is development-only", {
  text <- paste(readLines(file.path(root, "_targets.R"), warn = FALSE),
    collapse = "\n")
  expect_match(text, "06_annotation_models", fixed = TRUE)
  expect_match(text, "08_mt_validation", fixed = TRUE)
  expect_false(grepl("00_contract_smoke", text, fixed = TRUE))
  expect_false(grepl("05_ld_operator.*targets.R", text))
})

test_that("current studies publish honest status", {
  study06 <- paste(readLines(file.path(root, "studies", "06_annotation_models",
    "README.md"), warn = FALSE), collapse = "\n")
  expect_match(study06, "CLOSED — EST-R2", fixed = TRUE)
  expect_match(study06, "Official qualifier: EST-R5", fixed = TRUE)

  study08 <- paste(readLines(file.path(root, "studies", "08_mt_validation",
    "README.md"), warn = FALSE), collapse = "\n")
  expect_match(study08, "Status: In development", fixed = TRUE)
  expect_match(study08, "not a completed benchmark|no completed benchmark",
    ignore.case = TRUE)
})

test_that("five reusable templates use the shared runner", {
  templates <- c("finemapping-analysis.R", "prediction-analysis.R",
    "parameter-estimation-analysis.R", "convergence-analysis.R",
    "operator-analysis.R")
  for (template in templates) {
    text <- paste(readLines(file.path(root, "inst", "templates", template),
      warn = FALSE), collapse = "\n")
    expect_match(text, "run_benchmark", fixed = TRUE)
    expect_match(text, "workshop", fixed = TRUE)
    expect_match(text, "unsuitable|not\\s+(#\\s*)?suitable",
      ignore.case = TRUE)
    expect_match(text, "_plot <-")
  }
})

test_that("Study 06 validates while Study 08 remains unavailable", {
  args <- c("--profile", "benchmark", "--output-dir", "unused",
    "--resume", "true", "--validate-only", "true")
  parsed <- parse_benchmark_cli_arguments(c("--study",
    "06_annotation_models", args))
  expect_identical(parsed$study, "06_annotation_models")
  expect_true(parsed$validate_only)
  expect_error(parse_benchmark_cli_arguments(c("--study",
    "08_mt_validation", args)), "in development")
  expect_error(parse_benchmark_cli_arguments(c("--study",
    "08_mtblr_validation", args)), "retired.*08_mt_validation")
})
