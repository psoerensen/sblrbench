test_that("current refresh paths and package pin are explicit", {
  runner <- paste(readLines(test_path("..", "..", "scripts",
    "run_current_benchmark_refresh.R"), warn = FALSE), collapse = "\n")
  expect_match(runner, "02e8c74baa906e83c4a08d42a9cc6339b4e81072", fixed = TRUE)
  expect_match(runner, "current_benchmark_refresh", fixed = TRUE)
  expect_match(runner, "current-selection", fixed = TRUE)
  expect_false(grepl("pkgload::load_all", runner, fixed = TRUE))
})

test_that("five-replicate recommendations resolve to the current capsule", {
  env <- new.env(parent = baseenv())
  sys.source(test_path("..", "..", "studies", "five_replicate_helpers.R"), env)
  old <- Sys.getenv("SBLR_BENCH_RECOMMENDATIONS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("SBLR_BENCH_RECOMMENDATIONS") else
    Sys.setenv(SBLR_BENCH_RECOMMENDATIONS = old), add = TRUE)
  Sys.unsetenv("SBLR_BENCH_RECOMMENDATIONS")
  expect_match(env$.five_replicate_recommendation_path(), "current-selection")
  Sys.setenv(SBLR_BENCH_RECOMMENDATIONS = "custom/recommendations.csv")
  expect_identical(env$.five_replicate_recommendation_path(),
    "custom/recommendations.csv")
})
