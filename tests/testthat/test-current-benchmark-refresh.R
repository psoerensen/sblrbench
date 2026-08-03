test_that("current refresh paths and package pin are explicit", {
  runner <- paste(readLines(test_path("..", "..", "scripts",
    "run_current_benchmark_refresh.R"), warn = FALSE), collapse = "\n")
  expect_match(runner, "02e8c74baa906e83c4a08d42a9cc6339b4e81072", fixed = TRUE)
  expect_match(runner, "current_benchmark_refresh", fixed = TRUE)
  expect_match(runner, "current-selection", fixed = TRUE)
  expect_match(runner, "SBLR_BENCH_REPLICATES = \"10\"", fixed = TRUE)
  expect_match(runner, "pkgload::load_all(root", fixed = TRUE)
  expect_match(runner, "library(\"sblr\", lib.loc = rlib", fixed = TRUE)
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

test_that("Study 02 compacts simulation branches before aggregation", {
  source <- paste(readLines(test_path("..", "..", "studies", "02_prediction",
    "targets.R"), warn = FALSE), collapse = "\n")
  expect_match(source, "prediction_simulation_summary_branch")
  expect_match(source, "pattern = map(prediction_simulation_bundle)", fixed = TRUE)
  expect_false(grepl("lapply(prediction_simulation_bundle, function", source,
    fixed = TRUE))
})

test_that("Study 01 current capsule contract is exact and current", {
  promotion <- paste(readLines(test_path("..", "..", "studies", "01_finemapping",
    "promotion.R"), warn = FALSE), collapse = "\n")
  targets <- paste(readLines(test_path("..", "..", "studies", "01_finemapping",
    "targets.R"), warn = FALSE), collapse = "\n")
  expect_match(promotion, "results.*reference.*01_finemapping.*current")
  expect_match(promotion, "exact successful 40-fit grid", fixed = TRUE)
  expect_match(promotion, "target_warnings.csv", fixed = TRUE)
  expect_match(promotion, "warning_target_count", fixed = TRUE)
  expect_match(promotion, "ld_warning_validation.csv", fixed = TRUE)
  expect_match(promotion, "threshold_consistent", fixed = TRUE)
  expect_match(promotion, "02e8c74baa906e83c4a08d42a9cc6339b4e81072", fixed = TRUE)
  expect_match(targets, "sblr_source_commit")
  expect_match(targets, "simulation_summary")
  expect_match(targets, "seed_registry")
  expect_match(targets, "benchmark_summary")
  expect_match(targets, "computational_summary_branch")
  expect_match(targets, "marker_metrics_branch")
  expect_match(targets, "credible_set_metrics_branch")
  expect_match(targets, "credible_set_summary_branch")
  expect_false(grepl("lapply(method_run, `\\[\\[`, \"marker_metrics\")", targets))
})

test_that("Study 05 refresh promotes an explicit convergence decision only", {
  runner <- paste(readLines(test_path("..", "..", "scripts",
    "run_current_benchmark_refresh.R"), warn = FALSE), collapse = "\n")
  config <- paste(readLines(test_path("..", "..", "studies",
    "05_annotation_models", "config.R"), warn = FALSE), collapse = "\n")
  promotion <- paste(readLines(test_path("..", "..", "studies",
    "05_annotation_models", "promotion.R"), warn = FALSE), collapse = "\n")
  expect_match(runner, "SBLR_BENCH_STUDY05_PHASE = \"convergence\"", fixed = TRUE)
  expect_match(runner, "full_benchmark_started: false", fixed = TRUE)
  expect_match(config, "current-selection", fixed = TRUE)
  expect_match(promotion, "current-convergence", fixed = TRUE)
  expect_match(promotion, "current-stop", fixed = TRUE)
  expect_match(promotion, "prespecified_convergence_stop_triggered", fixed = TRUE)
  expect_match(promotion, "full_benchmark_started")
  diagnostics <- paste(readLines(test_path("..", "..", "studies",
    "05_annotation_models", "diagnostics.R"), warn = FALSE), collapse = "\n")
  expect_match(diagnostics, "recommendation_status = if (supported)", fixed = TRUE)
  expect_match(diagnostics, "\"unsupported\"", fixed = TRUE)
  expect_false(grepl("No supported annotation-model MCMC setting", diagnostics,
    fixed = TRUE))
})
