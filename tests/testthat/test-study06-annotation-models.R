study06_root <- if (file.exists(file.path("studies", "06_annotation_models",
  "config.R"))) "." else testthat::test_path("..", "..")
study06_dir <- file.path(study06_root, "studies", "06_annotation_models")
if (!file.exists(file.path(study06_dir, "config.R")))
  testthat::skip("repository-only Study 06 files are absent from package build")
cfg05 <- source(file.path(study06_dir, "config.R"), local = TRUE)$value
for (f in c("annotation_design.R", "simulation.R", "methods.R",
            "chain_extraction.R", "diagnostics.R", "metrics.R", "pilot.R",
            "promotion.R"))
  source(file.path(study06_dir, f), local = TRUE)

test_that("annotation construction is deterministic, aligned and standardized", {
  ids <- paste0("m", seq_len(1000))
  A1 <- .study06_annotation_design(ids, cfg05)
  A2 <- .study06_annotation_design(ids, cfg05)
  expect_identical(A1, A2)
  expect_identical(rownames(A1), ids)
  expect_identical(colnames(A1),
    c("Intercept", "enriched_binary", "continuous_signal", "null_annotation"))
  expect_equal(mean(A1[, "enriched_binary"]), .1)
  expect_equal(unname(colMeans(A1[, 3:4])), c(0, 0), tolerance = 1e-12)
  expect_equal(unname(apply(A1[, 3:4], 2, sd)), c(1, 1), tolerance = 1e-12)
  expect_equal(qr(A1)$rank, 4L)
  expect_error(.study06_validate_annotation(A1[-1, ], ids, cfg05),
    "contract failed")
})

test_that("exact marker probability transformation is bounded and calibrated", {
  ids <- paste0("m", seq_len(5000))
  A <- .study06_annotation_design(ids, cfg05)
  alpha <- .study06_true_alpha(A, cfg05)
  p1 <- .study06_marker_probabilities(A,
    alpha$informative_annotations, cfg05$mixture_var)
  p0 <- .study06_marker_probabilities(A,
    alpha$uninformative_annotations, cfg05$mixture_var)
  expect_equal(p1, sblr::sbayesrc_marker_pi(A,
    alpha$informative_annotations, cfg05$mixture_var))
  expect_true(all(p1 >= 0 & p1 <= 1))
  expect_equal(as.numeric(rowSums(p1)), rep(1, nrow(A)), tolerance = 1e-12)
  expect_equal(colMeans(p1), colMeans(p0), tolerance = 1e-10)
  expect_equal(sum(1 - p1[, 1]), 50, tolerance = .5)
  expect_true(alpha$enriched_expected_nonnull_share >= .5)
  expect_true(alpha$enriched_expected_nonnull_share <= .7)
  expect_equal(alpha$uninformative_annotations[-1, ], matrix(0, 3, 3,
    dimnames = dimnames(alpha$uninformative_annotations[-1, ])))
})

test_that("seed hierarchy is deterministic and method-order invariant", {
  a <- .study06_fit_seeds("informative_annotations", 2,
    "st_bed_bayesrc", cfg05)
  b <- .study06_fit_seeds("informative_annotations", 2,
    "st_bed_bayesrc", cfg05)
  expect_identical(a, b)
  expect_length(unique(a[paste0("chain_", 1:4)]), 4L)
  expect_false(identical(a, .study06_fit_seeds("informative_annotations", 2,
    "st_csr_sbayesrc", cfg05)))
  expect_identical(.study06_simulation_seeds("informative_annotations", 2, cfg05),
    .study06_simulation_seeds("informative_annotations", 2, cfg05))
})

test_that("simulation truth is internally consistent and targets training h2", {
  set.seed(1)
  Z <- matrix(rnorm(120 * 5000), 120, 5000,
    dimnames = list(paste0("i", 1:120), paste0("m", 1:5000)))
  A <- .study06_annotation_design(colnames(Z), cfg05)
  alpha <- .study06_true_alpha(A, cfg05)
  sim <- .study06_simulate("informative_annotations", 1, Z, 1:84,
    A, alpha, cfg05)
  expect_equal(sim$truth$heritability, .3, tolerance = 1e-10)
  expect_equal(as.numeric(Z %*% sim$effect), unname(sim$genetic_value),
    tolerance = 1e-10)
  expect_true(sim$truth$nonnull_count >= 20)
  expect_true(sim$truth$nonnull_count <= 100)
  expect_equal(as.numeric(rowSums(sim$prior_probability)), rep(1, ncol(Z)),
    tolerance = 1e-12)
})

test_that("metric orientation and zero-variance handling are explicit", {
  expect_equal(.study06_rmse(c(1, 2), c(2, 4)), sqrt(2.5))
  expect_equal(.study06_cor(1:4, 1:4), 1)
  expect_true(is.na(.study06_cor(rep(1, 4), 1:4)))
  expect_equal(.study06_auroc(c(.1, .9), c(FALSE, TRUE)), 1)
  expect_equal(.study06_auprc(c(.1, .9), c(FALSE, TRUE)), 1)
})

test_that("alpha and chain dimensions are preserved", {
  q <- data.frame(parameter_name = c("alpha", "sigmaSqAlpha"),
    annotation_name = c("enriched_binary", NA),
    stick_name = c("step_1", "step_1"),
    component_name = NA, pattern_name = NA)
  expect_identical(.study06_quantity_id(q),
    c("alpha:enriched_binary:step_1", "sigmaSqAlpha:step_1"))
  x <- expand.grid(iteration = 1:5, chain = 1:4,
    quantity = c("a", "b"))
  x$value <- seq_len(nrow(x))
  expect_equal(nrow(.study06_chain_window(x, 1, 4)), 32L)
  expect_error(.study06_chain_window(x[x$chain != 4, ], 1, 4),
    "four chains")
})

test_that("method and scenario grids are exact", {
  expect_identical(vapply(.study06_method_specs(cfg05), `[[`, "", "id"),
    cfg05$methods)
  specs <- .study06_specs(cfg05)
  expect_length(specs, 10L)
  expect_equal(as.numeric(table(vapply(specs, `[[`, "", "scenario"))), c(5, 5))
  expect_equal(nrow(.study06_seed_registry(cfg05)), 160L)
})

test_that("optional manifest booleans are report-safe", {
  manifest <- list()
  expect_false(isTRUE(manifest$some_optional_field))
  manifest$some_optional_field <- TRUE
  expect_true(isTRUE(manifest$some_optional_field))
})

test_that("capsule validators fail closed on missing fields", {
  path <- withr::local_tempdir()
  expect_error(.study06_validate_convergence_capsule(path), "incomplete")
  expect_error(.study06_validate_benchmark_capsule(path), "incomplete")
})

test_that("sampler-free installed annotation interface contract passes", {
  required <- c("stblr_bed", "stblr_csr", "stblr_csr_annot",
    "sbayesrc_marker_pi", "make_sbayesrc_alpha_init")
  expect_true(all(required %in% getNamespaceExports("sblr")))
  expect_true("bayesrc" %in% eval(formals(sblr::stblr_bed)$method))
  expect_true("annotation_probit_stick" %in%
    eval(formals(sblr::stblr_csr_annot)$annotation_model))
  A <- cbind(Intercept = 1, signal = c(-1, 0, 1))
  rownames(A) <- paste0("m", 1:3)
  alpha <- matrix(c(-3, -.25, -.67, .5, 0, 0), 2L, 3L,
    byrow = TRUE, dimnames = list(colnames(A), paste0("step_", 1:3)))
  p <- sblr::sbayesrc_marker_pi(A, alpha, c(0, .01, .1, 1))
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(unname(rowSums(p)), rep(1, nrow(A)), tolerance = 1e-12)
})

test_that("Study 06 remains explicitly in development", {
  readme <- paste(readLines(file.path(study06_dir, "README.md"), warn = FALSE),
    collapse = "\n")
  expect_match(readme, "Status: In development", fixed = TRUE)
  expect_match(readme, "not a completed benchmark", ignore.case = TRUE)
})
