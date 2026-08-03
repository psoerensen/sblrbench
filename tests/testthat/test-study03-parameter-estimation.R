study03_root <- testthat::test_path("..", "..")
study03_helpers <- file.path(study03_root, "studies", "03_parameter_estimation", "estimands.R")
if (!file.exists(study03_helpers)) {
  testthat::test_that("repository-only Study 03 helpers are excluded from package builds", {
    testthat::skip("repository-only Study 03 helpers are excluded from package builds")
  })
} else local({
  root <- testthat::test_path("..", "..")
  for (f in c("estimands.R", "simulation.R", "methods.R", "metrics.R", "pilot.R"))
    source(file.path(root, "studies", "03_parameter_estimation", f), local = TRUE)

  testthat::test_that("Study 03 registry is explicit and valid", {
    x <- .study03_validate_registry(.study03_estimand_registry())
    testthat::expect_length(unique(x$estimand_id), 6L)
    testthat::expect_equal(sum(x$primary), 4L)
    testthat::expect_true(all(grepl("st_bed|st_csr", x$available_methods)))
  })

  testthat::test_that("burn-in and nonlinear transforms are draw-wise", {
    fit <- list(input = list(nburn = 2L, nit = 4L, nthin = 1L),
      pi_trace = matrix(c(0, 0, .1, .9)), vbs = matrix(c(0, 0, 9, 1)),
      vgs = matrix(c(0, 0, 1, 3)), ves = matrix(c(0, 0, 3, 1)))
    d <- .study03_extract_draws(fit, "st_bed_bayesc", .study03_estimand_registry(), 10)
    tm <- d$value[d$estimand_id == "total_marker_effect_variance"]
    testthat::expect_equal(mean(tm), mean(c(9 * .1 * 10, 1 * .9 * 10)))
    testthat::expect_false(isTRUE(all.equal(mean(tm), mean(c(9, 1)) * mean(c(.1, .9)) * 10)))
    testthat::expect_equal(d$value[d$estimand_id == "heritability"], c(.25, .75))
    bad <- fit; bad$ves <- matrix(1, nrow = 3)
    testthat::expect_error(.study03_extract_draws(bad, "x", .study03_estimand_registry(), 10), "aligned")
  })

  testthat::test_that("truth construction preserves distinct realized quantities", {
    Z <- matrix(c(-1, 0, 1, -1, 1, 0), 3, 2,
      dimnames = list(paste0("i", 1:3), paste0("m", 1:2)))
    config <- source(file.path(root, "studies", "03_parameter_estimation", "config.R"), local = TRUE)$value
    config$simulation$n_causal <- 1L
    sim <- .study03_simulate(list(architecture = "sparse_homogeneous", replicate = 1L,
      simulation_seed = 9L), Z, config)
    truth <- .study03_truth(sim, config)
    testthat::expect_equal(truth$truth[truth$estimand_id == "causal_proportion"], .5)
    testthat::expect_equal(truth$truth[truth$estimand_id == "heritability"], .3, tolerance = 1e-12)
  })

  testthat::test_that("seeds ignore method order", {
    config <- source(file.path(root, "studies", "03_parameter_estimation", "config.R"), local = TRUE)$value
    specs <- .study03_replicate_specs(config); methods <- .study03_method_specs(config)
    a <- .study03_seed_registry(specs, methods, config)
    b <- .study03_seed_registry(specs, rev(methods), config)
    testthat::expect_equal(a[order(a$method), c("method", "fit_seed")],
      b[order(b$method), c("method", "fit_seed")], ignore_attr = TRUE)
  })

  testthat::test_that("paired differences require exact complete keys", {
    methods <- c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr")
    x <- data.frame(architecture = "sparse_homogeneous", replicate = 1L,
      estimand_id = "heritability", method = methods,
      posterior_mean = 1:4, absolute_error = 4:1, interval_width_95 = 1:4,
      covered_95 = c(TRUE, TRUE, FALSE, FALSE))
    testthat::expect_equal(nrow(.study03_paired(x)), 4L)
    testthat::expect_error(.study03_paired(x[-1, ]), "Incomplete")
  })

  testthat::test_that("tracked Study 03 capsule and report contract validate", {
    source(file.path(root, "studies", "03_parameter_estimation", "promotion.R"), local = TRUE)
    capsule <- file.path(root, "results", "reference", "03_parameter_estimation",
      "current")
    testthat::expect_silent(.study03_validate_capsule(capsule))
    report <- readLines(file.path(root, "studies", "03_parameter_estimation",
      "parameter-estimation.qmd"), warn = FALSE)
    testthat::expect_true(any(grepl("results.*reference", report)))
    testthat::expect_false(any(grepl("stblr_bed|stblr_csr|tar_make", report)))
  })
})
