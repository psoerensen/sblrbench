study05_root <- if (file.exists(file.path("studies", "05_annotation_models",
  "config.R"))) "." else testthat::test_path("..", "..")
study05_dir <- file.path(study05_root, "studies", "05_annotation_models")
if (!file.exists(file.path(study05_dir, "config.R")))
  testthat::skip("repository-only Study 05 files are absent from package build")
cfg05 <- source(file.path(study05_dir, "config.R"), local = TRUE)$value
for (f in c("annotation_design.R", "simulation.R", "methods.R",
            "chain_extraction.R", "diagnostics.R", "metrics.R", "pilot.R",
            "promotion.R"))
  source(file.path(study05_dir, f), local = TRUE)

test_that("annotation construction is deterministic, aligned and standardized", {
  ids <- paste0("m", seq_len(1000))
  A1 <- .study05_annotation_design(ids, cfg05)
  A2 <- .study05_annotation_design(ids, cfg05)
  expect_identical(A1, A2)
  expect_identical(rownames(A1), ids)
  expect_identical(colnames(A1),
    c("Intercept", "enriched_binary", "continuous_signal", "null_annotation"))
  expect_equal(mean(A1[, "enriched_binary"]), .1)
  expect_equal(unname(colMeans(A1[, 3:4])), c(0, 0), tolerance = 1e-12)
  expect_equal(unname(apply(A1[, 3:4], 2, sd)), c(1, 1), tolerance = 1e-12)
  expect_equal(qr(A1)$rank, 4L)
  expect_error(.study05_validate_annotation(A1[-1, ], ids, cfg05),
    "contract failed")
})

test_that("exact marker probability transformation is bounded and calibrated", {
  ids <- paste0("m", seq_len(5000))
  A <- .study05_annotation_design(ids, cfg05)
  alpha <- .study05_true_alpha(A, cfg05)
  p1 <- .study05_marker_probabilities(A,
    alpha$informative_annotations, cfg05$mixture_var)
  p0 <- .study05_marker_probabilities(A,
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
  a <- .study05_fit_seeds("informative_annotations", 2,
    "st_bed_bayesrc", cfg05)
  b <- .study05_fit_seeds("informative_annotations", 2,
    "st_bed_bayesrc", cfg05)
  expect_identical(a, b)
  expect_length(unique(a[paste0("chain_", 1:4)]), 4L)
  expect_false(identical(a, .study05_fit_seeds("informative_annotations", 2,
    "st_csr_sbayesrc", cfg05)))
  expect_identical(.study05_simulation_seeds("informative_annotations", 2, cfg05),
    .study05_simulation_seeds("informative_annotations", 2, cfg05))
})

test_that("simulation truth is internally consistent and targets training h2", {
  set.seed(1)
  Z <- matrix(rnorm(120 * 5000), 120, 5000,
    dimnames = list(paste0("i", 1:120), paste0("m", 1:5000)))
  A <- .study05_annotation_design(colnames(Z), cfg05)
  alpha <- .study05_true_alpha(A, cfg05)
  sim <- .study05_simulate("informative_annotations", 1, Z, 1:84,
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
  expect_equal(.study05_rmse(c(1, 2), c(2, 4)), sqrt(2.5))
  expect_equal(.study05_cor(1:4, 1:4), 1)
  expect_true(is.na(.study05_cor(rep(1, 4), 1:4)))
  expect_equal(.study05_auroc(c(.1, .9), c(FALSE, TRUE)), 1)
  expect_equal(.study05_auprc(c(.1, .9), c(FALSE, TRUE)), 1)
})

test_that("alpha and chain dimensions are preserved", {
  q <- data.frame(parameter_name = c("alpha", "sigmaSqAlpha"),
    annotation_name = c("enriched_binary", NA),
    stick_name = c("step_1", "step_1"),
    component_name = NA, pattern_name = NA)
  expect_identical(.study05_quantity_id(q),
    c("alpha:enriched_binary:step_1", "sigmaSqAlpha:step_1"))
  x <- expand.grid(iteration = 1:5, chain = 1:4,
    quantity = c("a", "b"))
  x$value <- seq_len(nrow(x))
  expect_equal(nrow(.study05_chain_window(x, 1, 4)), 32L)
  expect_error(.study05_chain_window(x[x$chain != 4, ], 1, 4),
    "four chains")
})

test_that("method and scenario grids are exact", {
  expect_identical(vapply(.study05_method_specs(cfg05), `[[`, "", "id"),
    cfg05$methods)
  specs <- .study05_specs(cfg05)
  expect_length(specs, 10L)
  expect_equal(as.numeric(table(vapply(specs, `[[`, "", "scenario"))), c(5, 5))
  expect_equal(nrow(.study05_seed_registry(cfg05)), 160L)
})

test_that("optional manifest booleans are report-safe", {
  manifest <- list()
  expect_false(isTRUE(manifest$some_optional_field))
  manifest$some_optional_field <- TRUE
  expect_true(isTRUE(manifest$some_optional_field))
})

test_that("capsule validators fail closed on missing fields", {
  path <- withr::local_tempdir()
  expect_error(.study05_validate_convergence_capsule(path), "incomplete")
  expect_error(.study05_validate_benchmark_capsule(path), "incomplete")
})

test_that("sampler-free installed interface contract passes", {
  source(file.path(study05_dir, "contract_smoke_test.R"), local = TRUE)
  expect_true(isTRUE(run_study05_contract_smoke_test()))
})

test_that("overnight runner validation mode cannot sample or mutate Git", {
  runner <- readLines(file.path(study05_root, "scripts",
    "run_study05_annotation_models.R"), warn = FALSE)
  launcher <- readLines(file.path(study05_root, "scripts",
    "run_study05_annotation_models.ps1"), warn = FALSE)
  validation_exit <- grep("if \\(validate_only\\)", runner)
  target_call <- grep("targets::tar_make", runner)
  expect_length(validation_exit, 1L)
  expect_length(target_call, 1L)
  expect_lt(validation_exit, target_call)
  expect_true(any(grepl("atomic_csv", runner, fixed = TRUE)))
  expect_true(any(grepl("study05_status.csv", runner, fixed = TRUE)))
  expect_true(any(grepl("Validated existing frozen capsule; target execution skipped",
    runner, fixed = TRUE)))
  prohibited <- "(git (add|commit|push|tag)|targets::tar_make.*preflight)"
  expect_false(any(grepl(prohibited, c(runner, launcher),
    ignore.case = TRUE)))
})

test_that("runner records failed phases in the atomic status contract", {
  runner <- paste(readLines(file.path(study05_root, "scripts",
    "run_study05_annotation_models.R"), warn = FALSE), collapse = "\n")
  expect_match(runner, 'if \\(ans\\$ok\\) "completed" else "failed"')
  expect_match(runner, "error_message = error")
  expect_match(runner, "file.rename\\(tmp, path\\)")
})
