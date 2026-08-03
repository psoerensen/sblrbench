test_that("prediction split is deterministic, complete, and ordered", {
  ids <- paste0("id", 1:10)
  a <- make_prediction_split(ids, 0.7, 31)
  b <- make_prediction_split(ids, 0.7, 31)
  expect_identical(a, b)
  expect_length(a$train_ids, 7L)
  expect_length(intersect(a$train_ids, a$test_ids), 0L)
  expect_setequal(c(a$train_ids, a$test_ids), ids)
  expect_true(all(diff(a$train_rows) > 0))
  expect_true(all(diff(a$test_rows) > 0))
  expect_error(make_prediction_split(ids, 0), "strictly")
  expect_error(make_prediction_split(ids, 1), "strictly")
})

test_that("genotype scaling is learned only from training rows", {
  x <- matrix(c(0, 0, 1, 1, 2, 2, 0, 1, 2, 1, 0, 2), 4, 3,
    dimnames = list(paste0("s", 1:4), paste0("m", 1:3)))
  a <- training_scaled_genotypes(x, 1:3)
  changed <- x
  changed[4, ] <- c(2, 0, 1)
  b <- training_scaled_genotypes(changed, 1:3)
  expect_equal(a$allele_frequency, colMeans(x[1:3, , drop = FALSE]) / 2)
  expect_identical(a$allele_frequency, b$allele_frequency)
  expect_identical(a$train, b$train)
  expect_identical(colnames(a$train), colnames(a$test))
  expect_equal(a$all, sweep(sweep(x, 2, a$center, "-"), 2, a$scale, "/"))
})

test_that("simulation views and prediction attachment preserve contracts", {
  s <- bench_fixture()
  view <- subset_sblrbench_simulation_samples(s, rev(s$data$sample_ids[1:2]))
  expect_identical(view$data$sample_ids, rev(s$data$sample_ids[1:2]))
  expect_identical(view$truth$effects, s$truth$effects)
  expect_error(subset_sblrbench_simulation_samples(s, c("absent")), "missing")
  expect_error(subset_sblrbench_simulation_samples(s,
    rep(s$data$sample_ids[[1]], 2)), "unique")
  base <- new_sblrbench_result("method", effects = s$truth$effects,
    provenance = list(note = "kept"))
  p <- s$truth$genetic_values[view$data$sample_ids, , drop = FALSE]
  out <- add_sblrbench_predictions(base, p, view)
  expect_identical(out$predictions$genetic_value, p)
  expect_identical(out$provenance, base$provenance)
})

test_that("prediction metrics are hand calculable and reject zero variance", {
  s <- bench_fixture()
  delta <- matrix(seq_len(length(s$truth$genetic_values)) / 10,
    nrow(s$truth$genetic_values), ncol(s$truth$genetic_values))
  predicted <- s$truth$genetic_values + delta
  r <- new_sblrbench_result("method", genetic_value = predicted)
  mse <- metric_prediction_mse(s, r)
  nmse <- metric_prediction_nmse(s, r)
  phen <- metric_phenotype_prediction_correlation(s, r)
  cal <- metric_prediction_calibration(s, r)
  expect_equal(mse$value,
    unname(colMeans((predicted - s$truth$genetic_values)^2)))
  expect_equal(nmse$value,
    unname(mse$value / apply(s$truth$genetic_values, 2, var)))
  expect_equal(phen$value, vapply(seq_len(ncol(predicted)), function(j)
    cor(predicted[, j], s$truth$phenotypes[, j]), numeric(1)))
  for (j in seq_along(s$data$trait_names)) {
    expected <- coef(lm(s$truth$genetic_values[, j] ~ predicted[, j]))
    z <- cal[cal$trait == s$data$trait_names[[j]], ]
    expect_equal(z$value, unname(expected))
  }
  flat <- predicted
  flat[, 1] <- 1
  failed <- metric_prediction_calibration(s,
    new_sblrbench_result("method", genetic_value = flat))
  expect_true(all(failed$status[failed$trait ==
    s$data$trait_names[[1]]] == "failed"))
  zero <- s
  zero$truth$genetic_values[, 1] <- 1
  expect_identical(metric_prediction_nmse(zero, r)$status[[1]], "failed")
})

test_that("generic paired advantages orient benefit and retain missing pairs", {
  metrics <- data.frame(architecture = "a", replicate = 1L,
    method = c("bayesc", "bayesr", "bayesc", "bayesr", "bayesc"), trait = "t",
    metric = c("prediction_correlation", "prediction_correlation",
      "prediction_mse", "prediction_mse", "effect_rmse"),
    value = c(.4, .6, 2, 1.5, .2), status = "ok")
  comparisons <- data.frame(comparison_id = "bayesr_vs_bayesc",
    focal_method = "bayesr", comparison_method = "bayesc")
  z <- paired_method_advantages(metrics, comparisons)
  expect_equal(z$advantage[z$paired_metric == "prediction_correlation"], .2)
  expect_equal(z$advantage[z$paired_metric == "prediction_mse"], .5)
  expect_false(z$complete_pair[z$paired_metric == "effect_rmse"])
  expect_true(is.na(z$advantage[z$paired_metric == "effect_rmse"]))
})

test_that("calibration advantages use absolute error", {
  metrics <- data.frame(architecture = "a", replicate = 1L,
    method = rep(c("bed", "csr"), each = 2), trait = "t",
    metric = rep(c("prediction_calibration_intercept",
      "prediction_calibration_slope"), 2), value = c(.2, 1.3, .1, .9),
    status = "ok")
  cmp <- data.frame(comparison_id = "csr_vs_bed", focal_method = "csr",
    comparison_method = "bed")
  z <- paired_method_advantages(metrics, cmp)
  expect_equal(z$advantage[
    z$paired_metric == "absolute_calibration_intercept_error"], .1)
  expect_equal(z$advantage[
    z$paired_metric == "absolute_calibration_slope_error"], .2)
})

.study02_spec_path <- function() testthat::test_path("..", "..", "studies",
  "02_prediction", "spec.R")

test_that("Study 02 spec and profiles validate with actionable failures", {
  path <- .study02_spec_path()
  skip_if_not(file.exists(path), "repository Study 02 spec is unavailable")
  spec <- read_benchmark_spec(path)
  expect_invisible(validate_benchmark_spec(spec))
  expect_identical(resolve_benchmark_profile(spec, "workshop")$replicate_count,
    1L)
  expect_identical(resolve_benchmark_profile(spec, "benchmark")$replicate_count,
    5L)
  incomplete <- spec
  incomplete$seeds <- NULL
  expect_error(validate_benchmark_spec(incomplete), "seeds")
  unsupported <- spec
  unsupported$task <- "parameter_estimation"
  expect_error(validate_benchmark_spec(unsupported), "Unsupported benchmark task")
})

test_that("Study 02 coordinates and all seed mappings are preserved", {
  spec <- read_benchmark_spec(.study02_spec_path())
  coordinates <- benchmark_seeds(spec, "benchmark")
  expect_equal(nrow(coordinates), 40L)
  expect_identical(unique(coordinates$scenario),
    c("sparse_homogeneous", "sparse_mixture"))
  expect_identical(unique(coordinates$replicate), 1:5)
  expect_identical(unique(coordinates$method), names(spec$methods))
  first <- coordinates[1, ]
  expect_identical(first$architecture_seed, 5001L)
  expect_identical(first$simulation_seed, 5002L)
  expect_identical(first$fit_seed, 20101L)
  expect_identical(first$chain_seeds[[1]],
    c(120101L, 220101L, 320101L, 420101L))
  last <- coordinates[nrow(coordinates), ]
  expect_identical(last$simulation_seed, 6006L)
  expect_identical(last$fit_seed, 30504L)
  expect_equal(nrow(benchmark_coordinates(spec, "workshop")), 8L)
})

test_that("migrated simulation kernel matches the old deterministic fixture", {
  spec <- read_benchmark_spec(.study02_spec_path())
  spec$controls$simulation$n_causal <- 10L
  set.seed(9182L)
  z <- matrix(rnorm(3000), 100, 30,
    dimnames = list(paste0("s", 1:100), paste0("m", 1:30)))
  expected <- c(
    sparse_homogeneous = "b88f71bf594303ad79d7b7f8924d5116ba3ceeb304bd836a300c1b3a1f873401",
    sparse_mixture = "0cfef6247c447a95ba0ee94c1728d1197965db74d278296a393c861e07150671")
  for (scenario in names(spec$scenarios)) {
    coordinate <- list(scenario = scenario, replicate = 1L,
      simulation_seed = 5001L + match(scenario, names(spec$scenarios)))
    simulation <- simulate_prediction_architecture(coordinate, z, spec)
    value <- digest::digest(list(effects = simulation$truth$effects,
      phenotypes = simulation$truth$phenotypes,
      causal = simulation$truth$causal, extras = simulation$extras),
      algo = "sha256")
    expect_identical(value, unname(expected[[scenario]]))
    expect_equal(simulation$truth$parameters$h2_observed, .30,
      tolerance = 1e-12)
    expect_true(check_oracle_genetic_values(simulation)$ok)
  }
})

test_that("Study 02 methods, priors, and controls are exact", {
  spec <- read_benchmark_spec(.study02_spec_path())
  methods <- resolve_prediction_methods(spec)
  expect_identical(vapply(methods, `[[`, character(1), "id"),
    names(spec$methods))
  coordinates <- benchmark_seeds(spec, "benchmark")
  for (method_id in names(spec$methods)) {
    row <- coordinates[coordinates$scenario == "sparse_homogeneous" &
      coordinates$replicate == 1L & coordinates$method == method_id, ]
    controls <- prediction_method_controls(spec, method_id, "benchmark",
      row$fit_seed, row$chain_seeds[[1]])
    expect_identical(controls$nburn, 250L)
    expect_identical(controls$nchains, 4L)
    expect_identical(controls$ncores, 4L)
    expect_identical(controls$nit,
      if (grepl("bayesr", method_id)) 2000L else 250L)
    if (grepl("bayesr", method_id)) {
      expect_equal(controls$pi, c(.99, rep(.01 / 3, 3)))
      expect_identical(controls$mixture_var, c(0, .01, .1, 1))
    } else expect_identical(controls$pi_init, .01)
  }
})

test_that("prediction summaries preserve one-replicate schemas", {
  metrics <- data.frame(architecture = "sparse_homogeneous", replicate = 1L,
    method = "st_bed_bayesc", metric = "prediction_correlation",
    value = .8, status = "ok")
  summary <- prediction_metric_summary(metrics)
  expect_equal(summary$value, .8)
  expect_true(is.na(summary$sd))
  paired <- data.frame(architecture = "sparse_homogeneous", replicate = 1L,
    comparison_id = "bayesr_vs_bayesc_bed",
    focal_method = "st_bed_bayesr", comparison_method = "st_bed_bayesc",
    paired_metric = "prediction_correlation", advantage = .1,
    complete_pair = TRUE)
  paired_summary <- prediction_paired_summary(paired)
  expect_true(is.na(paired_summary$sd_advantage))
  expect_equal(paired_summary$complete_pairs, 1L)
})

test_that("validate-only execution cannot call the fit dispatch", {
  spec <- read_benchmark_spec(.study02_spec_path())
  output <- tempfile("study02-validation-")
  withr::local_options(list(sblrbench.fit_dispatch = function(...)
    stop("FIT DISPATCH MUST NOT RUN")))
  result <- run_benchmark(spec, output, profile = "benchmark",
    resume = TRUE, validate_only = TRUE)
  expect_equal(nrow(result$status), 40L)
  expect_true(all(result$status$status == "not_run_validate_only"))
  expect_true(file.exists(file.path(output, "manifest.json")))
  expect_true(file.exists(file.path(output, "tables", "fit_status.csv")))
  expect_length(list.files(file.path(output, "checkpoints"), recursive = TRUE),
    0L)
})

test_that("CLI arguments are strict and default only booleans", {
  args <- c("--study", "02_prediction", "--profile", "benchmark",
    "--output-dir", "out")
  parsed <- sblrbench:::parse_benchmark_cli_arguments(args)
  expect_true(parsed$resume)
  expect_false(parsed$validate_only)
  expect_error(sblrbench:::parse_benchmark_cli_arguments(c(args,
    "--unknown", "x")), "Unknown")
  expect_error(sblrbench:::parse_benchmark_cli_arguments(c("--study",
    "03_parameter_estimation", "--profile", "benchmark", "--output-dir",
    "out")), "Unsupported --study")
  expect_error(sblrbench:::parse_benchmark_cli_arguments(c(args,
    "--resume", "yes")), "true.*false")
})

test_that("frozen Study 02 capsule validates against the migrated spec", {
  path <- testthat::test_path("..", "..", "results", "reference",
    "02_prediction", "current")
  skip_if_not(dir.exists(path), "frozen capsule is unavailable")
  spec <- read_benchmark_spec(.study02_spec_path())
  expect_invisible(sblrbench:::validate_prediction_capsule(path, spec))
})

test_that("removed Study 02 internals have no executable callers", {
  root <- testthat::test_path("..", "..")
  removed <- file.path(root, "studies", "02_prediction",
    c("config.R", "pilot.R", "targets.R", "promotion.R"))
  expect_false(any(file.exists(removed)))
  files <- list.files(root, pattern = "\\.(R|qmd)$", recursive = TRUE,
    full.names = TRUE)
  files <- files[!grepl("results[/\\\\](reference|local)|_freeze|_targets",
    files)]
  files <- files[basename(files) != "test-prediction-helpers.R"]
  text <- paste(unlist(lapply(files, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("studies[/\\\\]02_prediction[/\\\\](config|pilot|targets|promotion)\\.R",
    text))
  report <- readLines(file.path(root, "studies", "02_prediction", "report.qmd"),
    warn = FALSE)
  executable <- paste(report[!grepl("^```|^#|^The |^This |^A ", report)],
    collapse = "\n")
  expect_false(grepl("results/local|readRDS|run_benchmark|targets::|tar_read|stblr_(bed|csr)\\(",
    executable))
})
