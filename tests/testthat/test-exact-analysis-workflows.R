.workflow_spec <- function(study) read_benchmark_spec(testthat::test_path(
  "..", "..", "studies", study, "spec.R"))

test_that("shared design tables expose the committed scientific contracts", {
  prediction <- .workflow_spec("02_prediction")
  parameter <- .workflow_spec("03_parameter_estimation")

  scenarios <- benchmark_scenario_table(prediction, "benchmark")
  expect_identical(scenarios$scenario,
    c("sparse_homogeneous", "sparse_mixture"))
  expect_identical(scenarios$causal_markers, c(50L, 50L))
  expect_equal(scenarios$target_heritability, c(.3, .3))
  expect_identical(scenarios$replicate_count, c(5L, 5L))

  coordinates <- benchmark_coordinate_table(prediction, "benchmark")
  expect_equal(nrow(coordinates), 40L)
  expect_true(all(c("simulation_seed", "method_seed", "chain_seeds") %in%
    names(coordinates)))
  expect_identical(coordinates$simulation_seed[1L], 5002L)

  methods <- benchmark_method_table(prediction, "benchmark")
  expect_identical(methods$method, names(prediction$methods))
  expect_identical(methods$nburn, rep(250L, 4L))
  expect_identical(methods$nit, c(250L, 2000L, 250L, 2000L))
  expect_true(all(nzchar(methods$inclusion_prior)))
  expect_true(all(nzchar(methods$ld_policy)))

  data <- benchmark_data_summary(prediction, "benchmark")
  expect_identical(data$training_count, 3500L)
  expect_identical(data$test_count, 1500L)
  parameter_data <- benchmark_data_summary(parameter, "benchmark")
  expect_identical(parameter_data$sample_count, 5000L)
  expect_true(is.na(parameter_data$test_count))

  estimands <- benchmark_estimand_table(parameter)
  expect_identical(estimands$estimand, parameter$estimands$estimand_id)
  expect_identical(estimands$formula_or_fit_field,
    parameter$estimands$posterior_source)
  expect_identical(estimands$truth_definition,
    parameter$estimands$truth_source)
})

test_that("output inventory describes returned paths", {
  root <- tempfile("workflow-output-")
  dir.create(root)
  manifest <- file.path(root, "manifest.json")
  writeLines("{}", manifest)
  inventory <- benchmark_output_inventory(list(paths = list(root = root,
    manifest = manifest, absent = file.path(root, "absent.csv"))))
  expect_identical(inventory$output, c("root", "manifest", "absent"))
  expect_identical(inventory$exists, c(TRUE, TRUE, FALSE))
})

test_that("workflow plotting helpers return named ggplot objects", {
  metrics <- expand.grid(scenario = c("sparse_homogeneous",
    "sparse_mixture"), replicate = 1:2,
    method = sblrbench_method_levels(), stringsAsFactors = FALSE)
  metrics <- metrics[rep(seq_len(nrow(metrics)), each = 5L), ]
  metrics$metric <- rep(c("prediction_correlation", "prediction_nmse",
    "prediction_calibration_intercept", "prediction_calibration_slope",
    "effect_rmse"), times = nrow(metrics) / 5L)
  metrics$value <- seq_len(nrow(metrics)) / 100
  metrics$status <- "ok"
  expect_s3_class(plot_prediction_metrics(metrics), "ggplot")
  expect_s3_class(plot_prediction_calibration(metrics), "ggplot")
  expect_s3_class(plot_effect_recovery(metrics), "ggplot")

  estimates <- expand.grid(architecture = c("sparse_homogeneous",
    "sparse_mixture"), replicate = 1:2, method = sblrbench_method_levels(),
    estimand_id = c("heritability", "genetic_variance",
      "causal_proportion"), stringsAsFactors = FALSE)
  estimates$truth <- .3
  estimates$posterior_mean <- .3 + seq_len(nrow(estimates)) / 1000
  estimates$bias <- estimates$posterior_mean - estimates$truth
  estimates$status <- "ok"
  expect_s3_class(plot_parameter_recovery(estimates, "heritability"),
    "ggplot")
  expect_s3_class(plot_parameter_bias(estimates), "ggplot")
  expect_s3_class(plot_component_probabilities(estimates), "ggplot")

  runtime <- unique(metrics[c("scenario", "replicate", "method", "status")])
  runtime$elapsed_seconds <- seq_len(nrow(runtime))
  expect_s3_class(plot_benchmark_runtime(runtime), "ggplot")
})

test_that("exact analyses execute summaries and assign principal plots", {
  root <- testthat::test_path("..", "..")
  files <- file.path(root, "studies",
    c("02_prediction", "03_parameter_estimation"), "analysis.R")
  text <- lapply(files, readLines, warn = FALSE)
  executable <- lapply(text, function(x) x[!grepl("^\\s*#", x)])
  expect_true(all(vapply(executable, function(x)
    any(grepl("results <- run_benchmark", x, fixed = TRUE)), logical(1))))
  expect_true(all(vapply(executable, function(x)
    any(grepl("benchmark_scenario_table", x, fixed = TRUE)), logical(1))))
  expect_true(all(vapply(executable, function(x)
    any(grepl("benchmark_method_table", x, fixed = TRUE)), logical(1))))
  expect_true(all(vapply(executable, function(x)
    any(grepl("benchmark_output_inventory", x, fixed = TRUE)), logical(1))))

  prediction_names <- c("prediction_plot", "calibration_plot",
    "effect_recovery_plot", "runtime_plot")
  parameter_names <- c("heritability_recovery_plot",
    "genetic_variance_plot", "residual_variance_plot",
    "parameter_bias_plot", "component_probability_plot", "runtime_plot")
  expect_true(all(vapply(prediction_names, function(name)
    any(grepl(paste0(name, " <- plot_"), executable[[1L]], fixed = TRUE)),
    logical(1))))
  expect_true(all(vapply(parameter_names, function(name)
    any(grepl(paste0(name, " <- plot_"), executable[[2L]], fixed = TRUE)),
    logical(1))))
})

test_that("migration is complete, templates stay focused, and reports stay frozen", {
  root <- testthat::test_path("..", "..")
  migration <- paste(readLines(file.path(root, "docs", "dev",
    "study03_migration.md"), warn = FALSE), collapse = "\n")
  expect_match(migration, "now uses the common ordinary-R runner")
  expect_false(grepl("remain temporarily|pending an explicit", migration))

  templates <- file.path(root, "inst", "templates",
    c("prediction-analysis.R", "parameter-estimation-analysis.R"))
  template_text <- lapply(templates, readLines, warn = FALSE)
  expect_true(all(vapply(template_text, function(x)
    any(grepl("unsuitable", x, fixed = TRUE)), logical(1))))
  expect_true(all(vapply(template_text, function(x)
    any(grepl("results <- run_benchmark", x, fixed = TRUE)), logical(1))))
  expect_true(all(vapply(template_text, function(x)
    any(grepl("_plot <- plot_", x)), logical(1))))

  reports <- file.path(root, "studies",
    c("02_prediction", "03_parameter_estimation"), "report.qmd")
  for (report in reports) {
    lines <- readLines(report, warn = FALSE)
    fences <- cumsum(grepl("^```", lines))
    code <- paste(lines[fences %% 2L == 1L & !grepl("^```", lines)],
      collapse = "\n")
    expect_false(grepl("results/local|run_benchmark|readRDS|targets::|tar_read",
      code))
    expect_true(any(grepl("analysis.R", lines, fixed = TRUE)))
    expect_true(any(grepl("spec.R", lines, fixed = TRUE)))
    expect_true(any(grepl("templates/", lines, fixed = TRUE)))
  }
})
