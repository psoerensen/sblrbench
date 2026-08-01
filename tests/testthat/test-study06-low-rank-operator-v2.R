study06v2_root <- if (file.exists(file.path("studies", "06_ld_operator",
  "v2", "config.R"))) "." else file.path("..", "..")
study06v2_dir <- file.path(study06v2_root, "studies", "06_ld_operator", "v2")
if (dir.exists(study06v2_dir)) {
  study06v2_config <- source(file.path(study06v2_dir, "config.R"),
    local = TRUE)$value
  source(file.path(study06v2_dir, "guard.R"), local = TRUE)
  source(file.path(study06v2_dir, "design_crosswalk.R"), local = TRUE)
  source(file.path(study06v2_dir, "operator_validation.R"), local = TRUE)
  source(file.path(study06v2_dir, "methods.R"), local = TRUE)
  source(file.path(study06v2_dir, "provenance.R"), local = TRUE)
}

test_that("v2 pins the committed BayesR-calibration source revision", {
  skip_if_not(dir.exists(study06v2_dir))
  expect_identical(study06v2_config$required_sblr_sha,
    "96487b3194fc1f8c6789060da5f2e2a0eea89974")
  expect_identical(study06v2_config$required_sblr_version, "0.2.0")
  expect_match(study06v2_config$source_snapshot, "96487b3$")
})

test_that("v2 paths and capsules cannot overwrite v1", {
  skip_if_not(dir.exists(study06v2_dir))
  expect_false(study06v2_config$local_dir ==
    file.path("results", "local", "study06_ld_operator"))
  expect_true(all(grepl("low-rank.*v2|low_rank.*v2",
    c(study06v2_config$local_dir, study06v2_config$convergence_capsule,
      study06v2_config$benchmark_capsule))))
  expect_false(any(study06v2_config$historical_capsules %in%
    c(study06v2_config$convergence_capsule,
      study06v2_config$benchmark_capsule)))
})

test_that("operator-pilot gates are inherited explicitly from v1", {
  expect_true(is.list(study06v2_config$pilot_gate))
  expect_true(all(c("maximum_heritability_difference",
    "maximum_prediction_correlation_difference",
    "minimum_posterior_effect_correlation") %in%
    names(study06v2_config$pilot_gate)))
})

test_that("the exact 60-fit grid has explicit safe low-rank controls", {
  skip_if_not(dir.exists(study06v2_dir))
  grid <- .study06v2_validate_grid(study06v2_config)
  expect_equal(nrow(grid), 60L)
  expect_setequal(grid$configuration, study06v2_config$configurations)
  for (id in grep("^low_rank_", study06v2_config$configurations,
      value = TRUE)) {
    x <- .study06v2_low_rank_configuration(id, study06v2_config)
    expect_identical(x$representation, "low_rank")
    expect_true(is.finite(x$eigen_prop) && x$eigen_prop < 1)
    expect_false(any(c("eigen_filter", "eigen_tau", "eigen_eta") %in%
      names(x)))
  }
})

test_that("production methods explicitly carry retained-low-rank controls", {
  skip_if_not(dir.exists(study06v2_dir))
  grid <- .study06v2_method_grid(study06v2_config)
  low <- grid[startsWith(grid$configuration, "low_rank_"), ]
  expect_true(all(low$interface == "stblr_block_eigen"))
  expect_true(all(low$representation == "low_rank"))
  expect_true(all(is.finite(low$eigen_prop)))
  expect_false(any(c("eigen_filter", "eigen_tau", "eigen_eta") %in%
    names(low)))
})

test_that("legacy and reconstructed-dense controls fail statically", {
  skip_if_not(dir.exists(study06v2_dir))
  expect_error(.study06v2_assert_fit_spec("low_rank_0995",
    list(representation = "dense_reconstructed", eigen_prop = 0.995),
    study06v2_config), "prohibits dense_reconstructed")
  for (field in c("eigen_filter", "eigen_tau", "eigen_eta")) {
    controls <- list(representation = "low_rank", eigen_prop = 0.995)
    controls[[field]] <- 0
    expect_error(.study06v2_assert_fit_spec("low_rank_0995", controls,
      study06v2_config), "prohibit legacy eigen controls")
  }
})

test_that("v1-to-v2 crosswalk preserves the scientific design", {
  skip_if_not(dir.exists(study06v2_dir))
  x <- .study06v2_design_crosswalk(study06v2_config)
  expect_true(all(c("held_identical", "intentionally_changed",
    "no_longer_relevant") %in% x$disposition))
  expect_match(x$v2[x$field == "operator_contract"], "block_low_rank_v1")
})

test_that("Q and w deterministic oracle identities are exact", {
  skip_if_not(dir.exists(study06v2_dir))
  Q1 <- matrix(c(1, 0, .2, .8, -.1, .4), nrow = 2L)
  Q2 <- matrix(c(.7, .1, -.2, .6), nrow = 2L)
  w1 <- c(.4, -.3); w2 <- c(.2, .5)
  beta <- c(.1, -.2, .05, .3, -.1)
  factors <- list(Q1, Q2)
  scores <- list(matrix(w1, 1L), matrix(w2, 1L))
  projected <- c(as.numeric(crossprod(Q1, w1)),
    as.numeric(crossprod(Q2, w2)))
  diagonal <- c(colSums(Q1^2), colSums(Q2^2))
  residuals <- c(w1 - Q1 %*% beta[1:3], w2 - Q2 %*% beta[4:5])
  marker_residual <- c(as.numeric(crossprod(Q1, residuals[1:2])),
    as.numeric(crossprod(Q2, residuals[3:4])))
  inspect <- list(factor = factors, transformed_score = scores,
    projected_score = matrix(projected, 1L), diagonal = matrix(diagonal, 1L),
    residual = matrix(residuals, 1L), marker_residual = matrix(marker_residual, 1L),
    residual_offset = c(0L, 2L),
    quadratic_form = sum((Q1 %*% beta[1:3])^2) +
      sum((Q2 %*% beta[4:5])^2),
    projected_score_dot = sum(beta * projected),
    residual_norm_squared = sum(residuals^2),
    transformed_score_norm_squared = sum(w1^2) + sum(w2^2))
  yy <- 10
  gate <- .study06v2_low_rank_identity_gate(inspect, beta, yy,
    source_blocks = lapply(factors, crossprod), source_score = projected)
  expect_true(gate$summary$pass)
  expect_lte(gate$summary$projected_sse_identity_error, 1e-12)
  expect_true(gate$summary$projected_sse_residual > 0)
})

test_that("full-retention validation requires all positive rank", {
  skip_if_not(dir.exists(study06v2_dir))
  inspect <- list(diagnostics = data.frame(block_size = c(3L, 2L),
    positive_rank = c(3L, 2L), retained_rank = c(3L, 2L)))
  expect_no_error(.study06v2_validate_full_positive_rank(inspect,
    study06v2_config$eigen_prop_full))
  inspect$diagnostics$retained_rank[2] <- 1L
  expect_error(.study06v2_validate_full_positive_rank(inspect,
    study06v2_config$eigen_prop_full), "retain every positive")
})

test_that("runtime block sizes are defined by retained factor columns", {
  factors <- list(matrix(0, 3L, 5L), matrix(0, 2L, 4L))
  expect_identical(vapply(factors, ncol, integer(1)), c(5L, 4L))
})

test_that("named native chain lists retain identifiable indices", {
  chains <- setNames(lapply(1:4, function(i) list(chain_index = i)),
    paste0("task", 1:4))
  index <- unname(sort(vapply(chains, `[[`, integer(1), "chain_index")))
  expect_identical(index, 1:4)
})

test_that("v2 aggregation uses only retained-low-rank comparison labels", {
  source(file.path(study06v2_dir, "phases.R"), local = TRUE)
  registry <- .study06v2_comparison_registry()
  expect_setequal(registry$comparison_id, c(
    "full_csr_minus_bed", "block_csr_minus_full_csr",
    "low_rank_full_minus_block_csr",
    "low_rank_0995_minus_low_rank_full",
    "low_rank_0995_minus_low_rank_0999",
    "low_rank_0999_minus_low_rank_full",
    "low_rank_0995_minus_full_csr"))
  expect_false(any(grepl("block_eigen|hard|ridge", registry$comparison_id)))
})
