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
  expect_error(subset_sblrbench_simulation_samples(s, rep(s$data$sample_ids[[1]], 2)), "unique")
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
  expect_equal(mse$value, unname(colMeans((predicted - s$truth$genetic_values)^2)))
  expect_equal(nmse$value, unname(mse$value / apply(s$truth$genetic_values, 2, var)))
  expect_equal(phen$value, vapply(seq_len(ncol(predicted)), function(j)
    cor(predicted[, j], s$truth$phenotypes[, j]), numeric(1)))
  for (j in seq_along(s$data$trait_names)) {
    expected <- coef(lm(s$truth$genetic_values[, j] ~ predicted[, j]))
    z <- cal[cal$trait == s$data$trait_names[[j]], ]
    expect_equal(z$value, unname(expected))
  }
  flat <- predicted; flat[, 1] <- 1
  failed <- metric_prediction_calibration(s, new_sblrbench_result("method", genetic_value = flat))
  expect_true(all(failed$status[failed$trait == s$data$trait_names[[1]]] == "failed"))
  zero <- s; zero$truth$genetic_values[, 1] <- 1
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
      "prediction_calibration_slope"), 2), value = c(.2, 1.3, .1, .9), status = "ok")
  cmp <- data.frame(comparison_id = "csr_vs_bed", focal_method = "csr",
    comparison_method = "bed")
  z <- paired_method_advantages(metrics, cmp)
  expect_equal(z$advantage[z$paired_metric == "absolute_calibration_intercept_error"], .1)
  expect_equal(z$advantage[z$paired_metric == "absolute_calibration_slope_error"], .2)
})

test_that("single-trait architectures are deterministic and calibrated", {
  pilot_path <- testthat::test_path("..", "..", "studies", "02_prediction", "pilot.R")
  config_path <- testthat::test_path("..", "..", "studies", "02_prediction", "config.R")
  skip_if_not(file.exists(pilot_path) && file.exists(config_path),
    "repository-only Study 02 helpers are excluded from package builds")
  source(pilot_path, local = TRUE)
  config <- source(config_path, local = TRUE)$value
  Z <- matrix(stats::rnorm(3000), 100, 30,
    dimnames = list(paste0("s", 1:100), paste0("m", 1:30)))
  config$simulation$n_causal <- 10L
  for (architecture in names(config$simulation$architectures)) {
    spec <- list(architecture = architecture, replicate = 1L,
      simulation_seed = 5001L + match(architecture, names(config$simulation$architectures)))
    a <- .study02_simulate(spec, Z, config)
    b <- .study02_simulate(spec, Z, config)
    expect_identical(a$truth$effects, b$truth$effects)
    expect_length(a$truth$causal$all, 10L)
    expect_equal(a$truth$parameters$h2_observed, config$simulation$h2,
      tolerance = 1e-12)
    expect_true(check_oracle_genetic_values(a)$ok)
    expect_equal(nrow(a$extras$effect_components), 10L)
    if (architecture == "sparse_mixture")
      expect_true(all(grepl("^variance_", a$extras$effect_components$component)))
  }
})

test_that("active Study 02 method set is exactly the four ST methods", {
  pilot_path <- testthat::test_path("..", "..", "studies", "02_prediction", "pilot.R")
  config_path <- testthat::test_path("..", "..", "studies", "02_prediction", "config.R")
  skip_if_not(file.exists(pilot_path) && file.exists(config_path),
    "repository-only Study 02 helpers are excluded from package builds")
  source(pilot_path, local = TRUE)
  config <- source(config_path, local = TRUE)$value
  specs <- .study02_method_specs(config)
  expect_identical(vapply(specs, `[[`, character(1), "id"), config$methods)
  expect_false(any(grepl("^mt_", config$methods)))
  expect_false("multitrait" %in% names(config))
  expect_identical(config$reference_profiles$one_replicate_development$replicate_count, 1L)
})

test_that("promotion validation rejects active MT rows and partial grids", {
  promotion_path <- testthat::test_path("..", "..", "studies", "02_prediction", "promotion.R")
  config_path <- testthat::test_path("..", "..", "studies", "02_prediction", "config.R")
  skip_if_not(file.exists(promotion_path) && file.exists(config_path),
    "repository-only Study 02 helpers are excluded from package builds")
  source(promotion_path, local = TRUE)
  config <- source(config_path, local = TRUE)$value
  status <- expand.grid(architecture = names(config$simulation$architectures),
    replicate = 1L, method = config$methods, stringsAsFactors = FALSE)
  status$status <- "ok"; status$reason <- ""
  computation <- status; computation$runtime <- 1
  metrics <- merge(status[, c("architecture", "replicate", "method")],
    data.frame(metric = c("prediction_correlation", "prediction_mse",
      "prediction_nmse", "phenotype_prediction_correlation",
      "prediction_calibration_intercept", "prediction_calibration_slope",
      "effect_rmse")), by = NULL)
  metrics$trait <- "trait1"; metrics$value <- 1; metrics$status <- "ok"
  simulations <- data.frame(architecture = names(config$simulation$architectures),
    replicate = 1L, causal_count = 50L, realized_h2 = .3, oracle_ok = TRUE)
  paired <- expand.grid(architecture = names(config$simulation$architectures),
    comparison_id = c("bayesr_vs_bayesc_bed", "sbayesr_vs_sbayesc_csr",
      "csr_vs_bed_bayesc", "csr_vs_bed_bayesr"),
    paired_metric = paste0("metric", 1:7), stringsAsFactors = FALSE)
  paired$replicate <- 1L; paired$focal_method <- "focal"
  paired$comparison_method <- "comparison"; paired$complete_pair <- TRUE
  paired$advantage <- 0
  manifest <- list(task = "single_trait_prediction", replicate_count = 1L,
    benchmark_scope = "one_replicate_development", benchmark_status = "complete",
    active_methods = config$methods, training_sample_count = 3500L,
    test_sample_count = 1500L, canonical_marker_count = 37991L,
    expected_fit_count = 8L, successful_fit_count = 8L, failed_fit_count = 0L,
    qgdata = list(commit = config$example_data$commit, files = config$example_data$files))
  expect_invisible(.study02_validate_promotion_tables(config, metrics, paired,
    computation, status, simulations, manifest))
  bad <- status; bad$method[[1]] <- "mt_bed_bayesr"
  expect_error(.study02_validate_promotion_tables(config, metrics, paired,
    computation, bad, simulations, manifest), "MT rows")
  expect_error(.study02_validate_promotion_tables(config, metrics, paired,
    computation, status[-1, ], simulations, manifest), "eight complete")
  bad_manifest <- manifest; bad_manifest$replicate_count <- 10L
  expect_error(.study02_validate_promotion_tables(config, metrics, paired,
    computation, status, simulations, bad_manifest), "exactly one")
})

test_that("one-replicate summaries retain NA uncertainty", {
  promotion_path <- testthat::test_path("..", "..", "studies", "02_prediction", "promotion.R")
  skip_if_not(file.exists(promotion_path),
    "repository-only Study 02 helpers are excluded from package builds")
  source(promotion_path, local = TRUE)
  metrics <- data.frame(architecture = "sparse_homogeneous", replicate = 1L,
    method = "st_bed_bayesc", metric = "prediction_correlation",
    value = .8, status = "ok")
  summary <- .study02_benchmark_summary(metrics)
  expect_equal(summary$value, .8)
  expect_equal(summary$mean, .8)
  expect_true(is.na(summary$sd))
  expect_equal(summary$replicate_count, 1L)
  paired <- data.frame(architecture = "sparse_homogeneous", replicate = 1L,
    comparison_id = "bayesr_vs_bayesc_bed", focal_method = "st_bed_bayesr",
    comparison_method = "st_bed_bayesc", paired_metric = "prediction_correlation",
    advantage = .1, complete_pair = TRUE)
  paired_summary <- .study02_paired_summary(paired)
  expect_true(is.na(paired_summary$sd_advantage))
  expect_equal(paired_summary$complete_pairs, 1L)
})
