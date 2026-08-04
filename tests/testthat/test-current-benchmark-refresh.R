refresh_repository_only <- file.exists(test_path("..", "..", "scripts",
  "run_current_benchmark_refresh.R"))

test_that("current refresh paths and package pin are explicit", {
  skip_if_not(refresh_repository_only,
    "repository-only refresh sources are excluded from package builds")
  runner <- paste(readLines(test_path("..", "..", "scripts",
    "run_current_benchmark_refresh.R"), warn = FALSE), collapse = "\n")
  expect_match(runner, "02e8c74baa906e83c4a08d42a9cc6339b4e81072", fixed = TRUE)
  expect_match(runner, "current_benchmark_refresh", fixed = TRUE)
  expect_match(runner, "current-selection", fixed = TRUE)
  expect_match(runner, "run_common_benchmark(\"01_finemapping\"", fixed = TRUE)
  expect_match(runner, "pkgload::load_all(root", fixed = TRUE)
  expect_match(runner, "library(\"sblr\", lib.loc = rlib", fixed = TRUE)
})

test_that("five-replicate recommendations resolve to the current capsule", {
  skip_if_not(refresh_repository_only,
    "repository-only refresh sources are excluded from package builds")
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

test_that("Study 02 refresh uses the shared execution entry point", {
  skip_if_not(refresh_repository_only,
    "repository-only refresh sources are excluded from package builds")
  runner <- paste(readLines(test_path("..", "..", "scripts",
    "run_current_benchmark_refresh.R"), warn = FALSE), collapse = "\n")
  expect_match(runner, "read_benchmark_spec")
  expect_match(runner, "run_benchmark")
  expect_false(grepl("SBLR_BENCH_STUDY='02_prediction'", runner, fixed = TRUE))
})

test_that("Study 01 current capsule contract is exact and current", {
  skip_if_not(refresh_repository_only,
    "repository-only refresh sources are excluded from package builds")
  spec <- read_benchmark_spec(test_path("..","..","studies",
    "01_finemapping","spec.R"))
  report <- paste(readLines(test_path("..","..","studies",
    "01_finemapping","report.qmd"),warn=FALSE),collapse="\n")
  expect_equal(nrow(benchmark_coordinates(spec,"benchmark")),40L)
  expect_identical(spec$frozen_capsule,
    "results/reference/01_finemapping/current")
  expect_identical(spec$packages$sblr$sha,
    "02e8c74baa906e83c4a08d42a9cc6339b4e81072")
  expect_match(report,"target_warnings.csv",fixed=TRUE)
  expect_match(report,"ld_warning_validation.csv",fixed=TRUE)
  expect_false(grepl("run_benchmark\\(|readRDS\\(",report))
})

test_that("Study 05 refresh promotes an explicit convergence decision only", {
  skip_if_not(refresh_repository_only,
    "repository-only refresh sources are excluded from package builds")
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
