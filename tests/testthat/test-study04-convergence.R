study04_root <- testthat::test_path("..", "..")
study04_dir <- file.path(study04_root, "studies", "04_convergence")
if (!file.exists(file.path(study04_dir, "diagnostics.R"))) {
  testthat::test_that("repository-only Study 04 helpers are excluded from package builds", {
    testthat::skip("repository-only Study 04 helpers are excluded from package builds")
  })
} else local({
for (file in c("config.R", "diagnostic_registry.R", "chain_extraction.R",
               "diagnostics.R", "recommendations.R", "methods.R", "pilot.R")) {
  source(file.path(study04_dir, file), local = TRUE)
}

.study04_fixture <- function(n = 3000L, shifted = FALSE) {
  set.seed(22)
  values <- replicate(4L, stats::arima.sim(list(ar = .25), n = n))
  if (shifted) values[, 4L] <- values[, 4L] + 2
  data.frame(
    architecture = "sparse_homogeneous", method = "st_bed_bayesc",
    chain = rep(1:4, each = n), raw_iteration = rep(seq_len(n), 4L),
    effect_variance = as.vector(exp(values) / 10),
    genetic_variance = as.vector(exp(values) / 3),
    residual_variance = as.vector(exp(-values)), stringsAsFactors = FALSE
  ) |>
    transform(heritability = genetic_variance / (genetic_variance + residual_variance))
}

test_that("chain identity and window contracts are strict", {
  x <- .study04_fixture(50L)
  expect_equal(as.vector(table(.study04_window(x, 10L, 20L)$chain)), rep(20L, 4L))
  expect_error(.study04_window(x[x$chain != 4L, ], 10L, 20L), "chain identity")
  expect_error(.study04_window(rbind(x, x[1, ]), 10L, 20L), "chain identity")
  expect_error(.study04_window(x, 40L, 20L), "Insufficient")
})

test_that("derived quantities and native extraction preserve chain boundaries", {
  n <- 12L
  a <- array(NA_real_, c(n, 4L, 3L))
  for (ch in 1:4) { a[, ch, 1] <- ch + seq_len(n) / 100; a[, ch, 2] <- 2; a[, ch, 3] <- 3 }
  fit <- list(convergence_traces = list(values = a,
    quantities = data.frame(group = c("vbs", "vgs", "ves"))))
  x <- .study04_extract_chain_draws(fit, "sparse_homogeneous", "st_bed_bayesc")
  expect_equal(nrow(x), n * 4L)
  expect_equal(unique(x$heritability), .4)
  expect_equal(x$effect_variance[x$chain == 3][1], 3.01)
  expect_error(.study04_extract_chain_draws(fit, "a", "m", expected_chains = 2L), "Unexpected")
})

test_that("rank diagnostics respond to shifted and autocorrelated chains", {
  cfg <- source(file.path(study04_dir, "config.R"), local = TRUE)$value
  good <- .study04_diagnostics(.study04_fixture(), 250L, 1000L, .study04_registry(), cfg$thresholds)
  bad <- .study04_diagnostics(.study04_fixture(shifted = TRUE), 250L, 1000L, .study04_registry(), cfg$thresholds)
  expect_true(all(is.finite(good$rhat)))
  expect_true(any(bad$rhat > cfg$thresholds$rhat))
  reordered <- .study04_fixture(); reordered$chain <- c(4L, 2L, 1L, 3L)[reordered$chain]
  reordered <- reordered[order(reordered$chain, reordered$raw_iteration), ]
  again <- .study04_diagnostics(reordered, 250L, 1000L, .study04_registry(), cfg$thresholds)
  expect_equal(good$rhat, again$rhat, tolerance = 1e-12)
  constant <- .study04_fixture(); constant[5:8] <- 1
  ind <- .study04_diagnostics(constant, 250L, 500L, .study04_registry(), cfg$thresholds)
  expect_true(all(ind$status == "indeterminate"))
})

test_that("seeds and matched grid do not depend on order", {
  cfg <- source(file.path(study04_dir, "config.R"), local = TRUE)$value
  a <- .study04_chain_seeds("sparse_homogeneous", "st_bed_bayesc", cfg)
  b <- .study04_chain_seeds("sparse_homogeneous", rev(cfg$methods)[4], cfg)
  expect_identical(a, b)
  expect_equal(length(unique(a)), 4L)
  expect_equal(nrow(cfg$matched_grid), 4L)
  expect_equal(as.vector(table(cfg$matched_grid$architecture)), c(2L, 2L))
})

test_that("burn-in stability and stable checkpoint selection are explicit", {
  cfg <- source(file.path(study04_dir, "config.R"), local = TRUE)$value
  x <- .study04_fixture()
  d <- .study04_diagnostic_grid(x, cfg)
  s <- .study04_burnin_stability(x, d, cfg)
  expect_equal(nrow(s), 12L)
  expect_true(all(is.finite(s$standardized_mean_shift)))
  rec <- .study04_recommend(d, s, cfg)
  expect_equal(nrow(rec), 1L)
  expect_equal(rec$recommended_nthin, 1L)
  d$overall_pass[d$retained_draw_candidate == 500L] <- FALSE
  d$overall_pass[d$retained_draw_candidate == 1000L] <- TRUE
  rec2 <- .study04_recommend(d, transform(s, stable = TRUE), cfg)
  expect_false(identical(rec2$recommended_post_burnin_draws, 500))
})

test_that("the frozen capsule and report contract are strict", {
  old <- setwd(study04_root); on.exit(setwd(old), add = TRUE)
  source(file.path("studies", "04_convergence", "promotion.R"), local = TRUE)
  capsule <- file.path("results", "reference", "04_convergence",
    "st-multichain-convergence-development-v1")
  expect_true(.study04_validate_capsule(capsule))
  status <- read.csv(file.path(capsule, "chain_status.csv"))
  expect_equal(nrow(status), 16L)
  expect_equal(length(unique(status$chain)), 4L)
  report <- readLines(file.path("studies", "04_convergence", "convergence-development-pilot.qmd"), warn = FALSE)
  expect_true(any(grepl("results.*reference.*04_convergence", report)))
  expect_false(any(grepl("tar_make|stblr_bed|stblr_csr", report)))
})
})
