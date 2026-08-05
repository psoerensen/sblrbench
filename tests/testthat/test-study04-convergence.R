study04_root <- testthat::test_path("..", "..")
study04_spec <- read_benchmark_spec(file.path(study04_root, "studies",
  "04_convergence", "spec.R"))

.study04_fixture <- function(n = 3000L, shifted = FALSE) {
  set.seed(22)
  values <- replicate(4L, stats::arima.sim(list(ar = .25), n = n))
  if (shifted) values[, 4L] <- values[, 4L] + 2
  wide <- data.frame(chain = rep(1:4, each = n),
    iteration = rep(seq_len(n), 4L),
    effect_variance = as.vector(exp(values) / 10),
    genetic_variance = as.vector(exp(values) / 3),
    residual_variance = as.vector(exp(-values)), stringsAsFactors = FALSE)
  wide$heritability <- wide$genetic_variance /
    (wide$genetic_variance + wide$residual_variance)
  long <- reshape(wide, varying = c("effect_variance", "genetic_variance",
    "residual_variance", "heritability"), v.names = "value",
    timevar = "quantity", times = c("effect_variance", "genetic_variance",
      "residual_variance", "heritability"), direction = "long")
  long$stage <- "selection"
  long$scenario <- "sparse_homogeneous"
  long$replicate <- 1L
  long$method <- "st_bed_bayesc"
  rownames(long) <- NULL
  long
}

test_that("Study 04 spec and profiles preserve the historical matched grid", {
  expect_silent(validate_benchmark_spec(study04_spec))
  incomplete <- study04_spec
  incomplete$diagnostics$thresholds <- NULL
  expect_error(validate_benchmark_spec(incomplete), "diagnostics is missing")
  workshop <- benchmark_convergence_seeds(study04_spec, "workshop")
  benchmark <- benchmark_convergence_seeds(study04_spec, "benchmark")
  expect_equal(nrow(workshop), 4L)
  expect_equal(nrow(benchmark), 24L)
  expect_equal(as.vector(table(benchmark$stage)), c(4L, 20L))
  expect_identical(unique(workshop$simulation_seed), c(6002L, 7002L))
  expect_identical(workshop$chain_seeds[[1L]],
    c(1051100L, 1051200L, 1051300L, 1051400L))
  expect_identical(study04_spec$matched_grid$method,
    c("st_bed_bayesc", "st_csr_sbayesc", "st_bed_bayesr",
      "st_csr_sbayesr"))
})

test_that("true trace extraction preserves chain identity and derived h2", {
  n <- 12L
  values <- array(NA_real_, c(n, 4L, 3L))
  for (chain in 1:4) {
    values[, chain, 1L] <- chain + seq_len(n) / 100
    values[, chain, 2L] <- 2
    values[, chain, 3L] <- 3
  }
  fit <- list(convergence_traces = list(values = values,
    quantities = data.frame(group = c("vbs", "vgs", "ves"))))
  coordinate <- benchmark_convergence_seeds(study04_spec, "workshop")[1L, ]
  draws <- extract_convergence_traces(fit, coordinate,
    study04_spec$diagnostics$registry)
  expect_equal(nrow(draws), n * 4L * 4L)
  expect_equal(unique(draws$value[draws$quantity == "heritability"]), .4)
  expect_error(extract_convergence_traces(list(vgs = 2), coordinate,
    study04_spec$diagnostics$registry), "cannot be substituted")
})

test_that("windows and diagnostics reproduce the historical definitions", {
  draws <- .study04_fixture()
  window <- benchmark_chain_window(draws, 250L, 1000L,
    group_cols = "quantity")
  expect_true(all(table(window$quantity, window$chain) == 1000L))
  diagnostics <- benchmark_convergence_diagnostics(draws, 250L, 1000L,
    study04_spec$diagnostics$registry,
    study04_spec$diagnostics$thresholds)
  expect_equal(nrow(diagnostics), 4L)
  expect_true(all(is.finite(diagnostics$rhat)))
  shifted <- benchmark_convergence_diagnostics(.study04_fixture(shifted = TRUE),
    250L, 1000L, study04_spec$diagnostics$registry,
    study04_spec$diagnostics$thresholds)
  expect_true(any(shifted$rhat > study04_spec$diagnostics$thresholds$rhat))
  constant <- draws
  constant$value <- 1
  status <- benchmark_convergence_diagnostics(constant, 250L, 500L,
    study04_spec$diagnostics$registry,
    study04_spec$diagnostics$thresholds)
  expect_true(all(status$status == "indeterminate"))
})

test_that("candidate stability and recommendation rules remain exact", {
  draws <- .study04_fixture()
  diagnostics <- benchmark_convergence_candidate_grid(draws, study04_spec)
  stability <- benchmark_burnin_stability(draws, study04_spec)
  expect_equal(nrow(diagnostics), 48L)
  expect_equal(nrow(stability), 12L)
  recommendation <- benchmark_convergence_recommendations(diagnostics,
    stability, study04_spec)
  expect_equal(nrow(recommendation), 1L)
  expect_equal(recommendation$recommended_nthin, 1L)
  diagnostics$overall_pass[
    diagnostics$retained_draw_candidate == 500L] <- FALSE
  diagnostics$overall_pass[
    diagnostics$retained_draw_candidate == 1000L] <- TRUE
  later <- benchmark_convergence_recommendations(diagnostics,
    transform(stability, stable = TRUE), study04_spec)
  expect_false(identical(later$recommended_post_burnin_draws, 500))
})

test_that("shared recommendation logic reproduces the frozen selection", {
  capsule <- file.path(study04_root, "results", "reference",
    "04_convergence", "current-selection")
  diagnostics <- read.csv(file.path(capsule,
    "convergence_diagnostics.csv"), stringsAsFactors = FALSE)
  stability <- read.csv(file.path(capsule, "burnin_stability.csv"),
    stringsAsFactors = FALSE)
  names(diagnostics)[names(diagnostics) == "architecture"] <- "scenario"
  names(diagnostics)[names(diagnostics) == "estimand"] <- "quantity"
  names(stability)[names(stability) == "architecture"] <- "scenario"
  names(stability)[names(stability) == "estimand"] <- "quantity"
  observed <- benchmark_convergence_recommendations(diagnostics, stability,
    study04_spec)
  expected <- read.csv(file.path(capsule, "method_recommendations.csv"),
    stringsAsFactors = FALSE)
  columns <- c("method", "recommended_nburn",
    "recommended_post_burnin_draws", "recommended_nit_argument",
    "recommended_nthin", "recommended_nchains", "recommended_ncores")
  observed <- observed[match(expected$method, observed$method), columns]
  expected <- expected[, columns]
  rownames(observed) <- rownames(expected) <- NULL
  expect_equal(observed, expected)
})

test_that("validation-only convergence execution cannot dispatch fits", {
  old <- options(sblrbench.convergence_fit_dispatch = function(...)
    stop("FIT DISPATCH MUST NOT RUN"))
  on.exit(options(old), add = TRUE)
  out <- tempfile("study04-validation-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  spec <- with_current_sblr_description_pin(study04_spec)
  result <- run_benchmark(spec, out, "benchmark",
    validate_only = TRUE)
  expect_equal(nrow(result$status), 24L)
  expect_true(all(result$status$status == "not_run_validate_only"))
})

test_that("convergence plots return named ggplot objects", {
  draws <- .study04_fixture()
  diagnostics <- benchmark_convergence_candidate_grid(draws, study04_spec)
  stability <- benchmark_burnin_stability(draws, study04_spec)
  expect_s3_class(plot_convergence_rhat(diagnostics), "ggplot")
  expect_s3_class(plot_convergence_ess(diagnostics), "ggplot")
  expect_s3_class(plot_convergence_mcse(diagnostics), "ggplot")
  expect_s3_class(plot_convergence_stability(stability), "ggplot")
})

test_that("Study 04 workflow and report contracts are explicit", {
  analysis <- readLines(file.path(study04_root, "studies", "04_convergence",
    "analysis.R"), warn = FALSE)
  report <- readLines(file.path(study04_root, "studies", "04_convergence",
    "report.qmd"), warn = FALSE)
  expect_true(any(grepl("run_benchmark", analysis, fixed = TRUE)))
  for (name in c("rhat_plot", "ess_plot", "mcse_plot",
      "relative_mcse_plot", "stability_plot", "runtime_plot"))
    expect_true(any(grepl(paste0(name, " <-"), analysis, fixed = TRUE)))
  expect_true(any(grepl("results.*reference.*04_convergence", report)))
  expect_false(any(grepl("results/local|readRDS|run_benchmark|tar_make|stblr_bed|stblr_csr",
    report)))
})

test_that("obsolete Study 04 source paths have no active callers", {
  obsolete <- c("config.R", "targets.R", "validation_targets.R",
    "chain_extraction.R", "diagnostic_registry.R", "diagnostics.R",
    "methods.R", "pilot.R", "recommendations.R", "promotion.R",
    "run_convergence_benchmark.R", "worked_convergence_example.R",
    "convergence_contract_smoke_test.R")
  expect_false(any(file.exists(file.path(study04_root, "studies",
    "04_convergence", obsolete))))
})
