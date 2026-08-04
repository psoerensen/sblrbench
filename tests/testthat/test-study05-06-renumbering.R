test_that("final study numbering and integrated paths are authoritative", {
  root <- testthat::test_path("..", "..")
  expect_true(dir.exists(file.path(root, "studies", "05_ld_operator")))
  expect_true(dir.exists(file.path(root, "studies", "06_annotation_models")))
  expect_false(dir.exists(file.path(root, "studies", "06_ld_operator")))
  expect_false(dir.exists(file.path(root, "studies", "05_annotation_models")))
  expect_true(file.exists(file.path(root, "results", "reference",
    "05_ld_operator", "current", "sbayesr_variance_summary.csv")))
  expect_false(dir.exists(file.path(root, "results", "reference",
    "05_ld_operator", "sbayesr_ld_robustness")))
})

test_that("retired CLI IDs are rejected with their replacements", {
  expect_error(parse_benchmark_cli_arguments(c("--study", "06_ld_operator",
    "--profile", "benchmark", "--output-dir", "out")),
    "05_ld_operator", fixed = TRUE)
  expect_error(parse_benchmark_cli_arguments(c("--study",
    "05_annotation_models", "--profile", "benchmark", "--output-dir",
    "out")), "06_annotation_models", fixed = TRUE)
})

test_that("website has one Study 05 report and no additional validation", {
  root <- testthat::test_path("..", "..")
  quarto <- paste(readLines(file.path(root, "_quarto.yml"), warn = FALSE),
    collapse = "\n")
  expect_match(quarto, "studies/05_ld_operator/report.qmd", fixed = TRUE)
  expect_false(grepl("Additional validation", quarto, fixed = TRUE))
  expect_false(grepl("sbayesr_ld_robustness/report.qmd", quarto,
    fixed = TRUE))
  annotation <- paste(readLines(file.path(root, "studies",
    "06_annotation_models", "annotation-convergence.qmd"), warn = FALSE),
    collapse = "\n")
  expect_match(annotation, "In development", fixed = TRUE)
})
