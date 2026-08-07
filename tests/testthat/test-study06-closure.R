test_that("Study 06 final decision is closed at EST-R2", {
  path <- test_path("..", "..", "results", "reference", "06_annotation_models", "final_decision.json")
  decision <- jsonlite::read_json(path, simplifyVector = TRUE)

  expect_identical(decision$status, "closed")
  expect_identical(decision$primary_decision, "EST-R2")
  expect_identical(decision$official_qualifier, "EST-R5")
  expect_identical(decision$specification_hash,
    "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56")
  expect_identical(decision$truth_hash,
    "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb")
  expect_true(decision$study06_closed)
  expect_false(decision$additional_final_benchmark_authorized)
  expect_false(decision$additional_scientific_fits_required)
})

test_that("Study 06 compact closure evidence has required comparisons", {
  base <- test_path("..", "..", "results", "reference", "06_annotation_models")
  comparison <- read.csv(file.path(base, "final_cross_implementation_comparison.csv"),
    check.names = FALSE)
  hierarchy <- read.csv(file.path(base, "final_hierarchy_of_evidence.csv"),
    check.names = FALSE)

  expect_setequal(
    comparison$quantity[comparison$section == "raw_alpha_stick1"],
    c("intercept", "enriched_binary", "continuous_signal", "null_annotation")
  )
  expect_setequal(
    comparison$quantity[comparison$section == "active_contrast"],
    c("enriched_binary", "continuous_signal", "null_annotation")
  )
  raw <- comparison[comparison$section == "raw_alpha_stick1", ]
  expect_equal(as.numeric(raw$truth[raw$quantity == "enriched_binary"]), 1.6)
  expect_equal(as.numeric(raw$official_sbayesrc_D1[raw$quantity == "continuous_signal"]),
    0.585737, tolerance = 1e-6)
  expect_equal(as.numeric(raw$sblr_learned_BED[raw$quantity == "null_annotation"]),
    -0.042171, tolerance = 1e-6)
  expect_equal(as.numeric(raw$sblr_learned_block[raw$quantity == "enriched_binary"]),
    1.642282, tolerance = 1e-6)
  expect_true(all(c("raw alpha", "q / pi probabilities", "SNP PIP") %in% hierarchy$quantity))
})

test_that("authoritative Study 06 entry point names the final decision", {
  path <- test_path("..", "..", "studies", "06_annotation_models", "README.md")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "CLOSED — EST-R2", fixed = TRUE)
  expect_match(text, "EST-R5", fixed = TRUE)
})
