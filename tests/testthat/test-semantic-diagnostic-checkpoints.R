test_that("semantic checkpoint identity ignores source and documentation paths", {
  inputs <- list(
    study = "03_parameter_estimation",
    scenario = "sparse_mixture",
    replicate = 1L,
    marker_ids = c("m1", "m2", "m3"),
    method_controls = list(nburn = 250L, nit = 2000L),
    simulation_seed = 7002L,
    fit_seed = 40104L,
    sblr_sha = "02e8c74baa906e83c4a08d42a9cc6339b4e81072"
  )
  identity <- benchmark_semantic_checkpoint_identity("diagnostic", inputs)
  expect_identical(benchmark_semantic_checkpoint_hash(identity),
    benchmark_semantic_checkpoint_hash(
      benchmark_semantic_checkpoint_identity("diagnostic", inputs)))

  # Paths outside the scientific payload cannot influence its identity.
  source_path_a <- file.path("studies", "03_parameter_estimation", "config.R")
  source_path_b <- "R/benchmark-spec.R"
  report_path_a <- "old/report.qmd"
  report_path_b <- "studies/03_parameter_estimation/report.qmd"
  expect_false(identical(source_path_a, source_path_b))
  expect_false(identical(report_path_a, report_path_b))
  expect_identical(benchmark_semantic_checkpoint_hash(identity),
    benchmark_semantic_checkpoint_hash(identity))
  expect_error(benchmark_semantic_checkpoint_identity("diagnostic",
    c(inputs, list(source_path = source_path_a))), "cannot contain")
  expect_error(benchmark_semantic_checkpoint_identity("diagnostic",
    c(inputs, list(report_path = report_path_a))), "cannot contain")
})

test_that("semantic checkpoint identity changes with scientific inputs", {
  inputs <- list(marker_ids = c("m1", "m2"),
    controls = list(nit = 2000L), simulation_seed = 7002L,
    fit_seed = 40104L, sblr_sha = "sha-a")
  hash <- function(x) benchmark_semantic_checkpoint_hash(
    benchmark_semantic_checkpoint_identity("diagnostic", x))
  baseline <- hash(inputs)
  changed <- function(name, value) {
    x <- inputs
    x[[name]] <- value
    hash(x)
  }
  expect_false(identical(baseline,
    changed("controls", list(nit = 2001L))))
  expect_false(identical(baseline, changed("simulation_seed", 7003L)))
  expect_false(identical(baseline, changed("fit_seed", 40105L)))
  expect_false(identical(baseline, changed("marker_ids", c("m2", "m1"))))
  expect_false(identical(baseline, changed("sblr_sha", "sha-b")))
})

test_that("legacy diagnostic checkpoints are rejected", {
  path <- tempfile(fileext = ".rds")
  saveRDS(list(schema_version = 1L, source_hashes = "retired",
    identity_hash = "legacy"), path)
  expect_error(benchmark_load_semantic_checkpoint(path, "unused"),
    "Legacy source-hashed diagnostic checkpoint detected")
})

test_that("active diagnostics do not reference retired Study 03 sources", {
  root <- testthat::test_path("..", "..")
  scripts <- file.path(root,
    "studies/03_parameter_estimation/diagnostics/sbayesr-gctb-comparison.R")
  text <- paste(unlist(lapply(scripts, readLines, warn = FALSE)), collapse = "\n")
  retired <- paste0("studies/03_parameter_estimation/",
    c("config.R", "simulation.R", "methods.R", "pilot.R"))
  expect_false(any(vapply(retired, grepl, logical(1), x = text,
    fixed = TRUE)))
  expect_false(dir.exists(file.path(root, "studies", "05_ld_operator",
    "sbayesr_ld_robustness")))
})
