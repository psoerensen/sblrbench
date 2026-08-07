test_that("Study 06 alpha-recovery addendum is registered without reopening closure", {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  metadata_path <- file.path(
    repo_root, "results", "reference", "06_annotation_models", "alpha_recovery_metadata.json"
  )
  final_path <- file.path(
    repo_root, "results", "reference", "06_annotation_models", "final_decision.json"
  )

  metadata <- jsonlite::read_json(metadata_path, simplifyVector = TRUE)
  final <- jsonlite::read_json(final_path, simplifyVector = TRUE)

  expect_identical(metadata$decision, "ALPHA-R1")
  expect_identical(metadata$replicates, 20L)
  expect_identical(metadata$chains, 4L)
  expect_identical(metadata$iterations, 8000L)
  expect_identical(metadata$burn_in, 2000L)
  expect_lt(metadata$fixed_z$maximum_posterior_mean_error, 0.005)
  expect_lt(metadata$fixed_z$maximum_posterior_covariance_error, 0.001)
  expect_lt(metadata$single_fixture_blocked_vs_scalar$maximum_absolute_posterior_mean_difference, 0.01)
  expect_lt(metadata$repeated_recovery$maximum_rhat, 1.0021)

  expect_identical(final$status, "closed")
  expect_identical(final$primary_decision, "EST-R2")
  expect_identical(final$official_qualifier, "EST-R5")
  expect_identical(final$sampler_development_status, "closed_at_PMA-R3")
  supporting_decisions <- if (is.data.frame(final$supporting_evidence)) {
    final$supporting_evidence$decision
  } else {
    vapply(final$supporting_evidence, `[[`, character(1), "decision")
  }
  expect_true("ALPHA-R1" %in% supporting_decisions)
})

test_that("tracked repeated alpha recovery satisfies ALPHA-R1 gates", {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  path <- file.path(
    repo_root, "results", "reference", "06_annotation_models", "alpha_recovery_summary.csv"
  )
  x <- read.csv(path, check.names = FALSE)

  expect_equal(nrow(x), 12L)
  expect_setequal(x$stick, 1:3)
  expect_setequal(
    x$parameter,
    c("Intercept", "enriched_binary", "continuous_signal", "null_annotation")
  )
  expect_lte(max(abs(x$bias)), 0.066)
  expect_gte(min(x$coverage_95), 0.90)
  expect_lte(max(x$max_rhat), 1.0021)

  enriched <- x[x$parameter == "enriched_binary", ]
  enriched <- enriched[order(enriched$stick), ]
  expect_true(all(diff(enriched$mean_post_sd) > 0))
  expect_equal(enriched$mean_eligible, c(1500, 170.7, 84.2), tolerance = 1e-8)
  expect_equal(enriched$mean_enriched_eligible, c(149.05, 77.2, 42.3), tolerance = 1e-8)
  expect_equal(enriched$bias, c(-0.017151477, -0.001762922, -0.012222834), tolerance = 1e-6)
})

test_that("blocked and scalar known-outcome summaries agree", {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
  path <- file.path(
    repo_root, "results", "reference", "06_annotation_models", "alpha_recovery_single_fixture.csv"
  )
  x <- read.csv(path, check.names = FALSE)
  blocked <- x[x$update == "blocked", c("stick", "parameter", "posterior_mean", "posterior_sd")]
  scalar <- x[x$update == "scalar", c("stick", "parameter", "posterior_mean", "posterior_sd")]
  names(blocked)[3:4] <- c("blocked_mean", "blocked_sd")
  names(scalar)[3:4] <- c("scalar_mean", "scalar_sd")
  paired <- merge(blocked, scalar, by = c("stick", "parameter"))

  expect_equal(nrow(paired), 12L)
  expect_lt(max(abs(paired$blocked_mean - paired$scalar_mean)), 0.01)
  expect_lt(max(abs(paired$blocked_sd - paired$scalar_sd)), 0.004)
  expect_lt(max(x$rhat), 1.002)
})
