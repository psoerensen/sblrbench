study06_source_dir <- testthat::test_path("..", "..", "studies",
  "06_ld_operator")
if (!dir.exists(study06_source_dir)) {
  testthat::test_that("Study 06 sources are excluded from package builds", {
    testthat::skip("repository-only Study 06 helpers are excluded from package builds")
  })
} else local({
study06_files <- c("blocks.R", "operators.R", "operator_validation.R",
  "simulation.R", "methods.R", "chain_extraction.R",
  "diagnostics.R", "metrics.R", "pilot.R")
for (file in study06_files)
  source(testthat::test_path("..", "..", "studies",
    "06_ld_operator", file), local = TRUE)

study06_config <- source(testthat::test_path("..", "..", "studies",
  "06_ld_operator", "config.R"), local = TRUE)$value

test_that("Study 06 block construction is deterministic and complete", {
  ids <- paste0("rs", seq_len(1L + 2L * 500L))
  a <- .study06_blocks(ids, 500L)
  b <- .study06_blocks(ids, 500L)
  expect_identical(a, b)
  expect_silent(.study06_validate_blocks(a, ids))
  expect_identical(unlist(Map(seq.int, a$start, a$end),
    use.names = FALSE), seq_along(ids))
  expect_identical(a$size, c(500L, 500L, 1L))
  broken <- a
  broken$start[2L] <- broken$start[2L] + 1L
  expect_error(.study06_validate_blocks(broken, ids),
    "gap, overlap|coverage")
})

test_that("Study 06 method and seed grids are fixed and order invariant", {
  grid <- .study06_method_grid(study06_config)
  expect_equal(nrow(grid), 12L)
  expect_setequal(grid$configuration,
    study06_config$configurations)
  expected <- expand.grid(architecture = study06_config$architectures,
    replicate = 1:5,
    configuration = study06_config$configurations,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  expect_equal(nrow(expected), 60L)
  seed_a <- .study06_seed_registry(study06_config)
  reversed <- study06_config
  reversed$configurations <- rev(reversed$configurations)
  seed_b <- .study06_seed_registry(reversed)
  key <- c("architecture", "replicate", "configuration", "chain")
  seed_a <- seed_a[do.call(order, seed_a[key]), ]
  seed_b <- seed_b[do.call(order, seed_b[key]), ]
  expect_identical(seed_a$simulation_seed, seed_b$simulation_seed)
  expect_identical(anyDuplicated(seed_a[c(key, "chain_seed")]), 0L)
  expect_true(all(vapply(split(seed_a$chain_seed,
    interaction(seed_a$architecture, seed_a$replicate,
      seed_a$configuration, drop = TRUE)),
    function(x) length(unique(x)) == 4L, logical(1))))
})

test_that("packed reconstruction, products, and quadratic forms agree", {
  A <- matrix(c(4, 1, .5, 1, 3, -.2, .5, -.2, 2),
    3, 3)
  packed <- .study06_pack_triangle(A)
  expect_equal(.study06_unpack_triangle(packed, 3L), A)
  blocks <- list(A, diag(c(2, 5)))
  x <- c(1, -2, .5, 3, -1)
  expected <- c(A %*% x[1:3], diag(c(2, 5)) %*% x[4:5])
  expect_equal(.study06_apply_blocks(blocks, x),
    as.numeric(expected))
  expect_equal(drop(crossprod(x,
    .study06_apply_blocks(blocks, x))),
    drop(crossprod(x[1:3], A %*% x[1:3])) +
      drop(crossprod(x[4:5], diag(c(2, 5)) %*% x[4:5])))
})

test_that("cross-block LD accounting uses the fixed orientation", {
  csr <- list(row_ptr = c(0, 2, 3, 4, 4),
    col_idx = c(1L, 2L, 3L, 3L),
    values = c(.5, .2, .3, .4))
  x <- .study06_cross_block_summary(csr, 4L, 2L)
  expect_equal(x$cross_block_edge_count_removed, 2L)
  expect_equal(x$cross_block_squared_ld_mass_removed,
    .2^2 + .3^2)
  expect_true(x$fraction_full_csr_squared_ld_mass_retained < 1)
})

test_that("the frozen public hard route accepts an effective no-op", {
  x <- data.frame(eigen_tau = c(.02, .10),
    retained_rank = c(100, 100), original_rank = 100,
    retained_rank_proportion = c(1, 1),
    retained_positive_eigenvalue_mass = c(100, 100),
    positive_eigenvalue_mass = 100,
    retained_mass_proportion = c(1, 1),
    minimum_runtime_diagonal = 1,
    operator_frobenius_maximum_error = 0,
    matrix_vector_maximum_error = 0,
    quadratic_form_maximum_error = 0,
    eigenvalues_removed = 0,
    filter_activity_status = "effective_no_op",
    prediction_check_status = "fixture",
    pass = TRUE)
  selected <- .study06_select_hard_filter(x, .10)
  expect_equal(selected$eigen_tau, .10)
  expect_equal(selected$filter_activity_status, "effective_no_op")
})

test_that("synthetic spectra exercise active hard truncation exactly", {
  x <- .study06_synthetic_filter_validation()
  expect_true(all(x$pass))
  controlled <- x[x$matrix_id == "controlled_small", ]
  expect_equal(controlled$retained_rank[
    controlled$effective_threshold == .01][1L], 4L)
  expect_equal(controlled$retained_rank[
    controlled$effective_threshold == .10], 3L)
  expect_true(all(controlled$zero_ridge_maximum_error < 1e-10))
  expect_true(all(controlled$filtered_reconstruction_maximum_error < 1e-10))
  expect_true(all(controlled$minimum_diagonal > 0))
})

test_that("well-conditioned, near-singular, and multiblock contracts are stable", {
  A <- .study06_known_spectrum_matrix(c(5, 3, 2, 1, .5))
  no_op <- .study06_filter_symmetric_spectrum(A, .1)
  expect_equal(no_op$retained_rank, 5L)
  expect_equal(no_op$matrix, A, tolerance = 1e-10)
  C <- .study06_known_spectrum_matrix(c(5, 2, 1, .05, 1e-8))
  fixed <- .study06_filter_dense_matrix(C, "ridge_fixed", eta = .01 / .99)
  lw <- .study06_filter_dense_matrix(C, "ridge_lw", lw_shrinkage = .75)
  expect_true(all(is.finite(fixed$matrix)) && all(diag(fixed$matrix) > 0))
  expect_true(all(is.finite(lw$matrix)) && all(diag(lw$matrix) > 0))
  blocks <- list(A, C)
  v <- seq_len(10) / 10
  expect_equal(.study06_apply_blocks(blocks, v),
    c(A %*% v[1:5], C %*% v[6:10]))
  dense <- matrix(0, 10, 10)
  dense[1:5, 1:5] <- A; dense[6:10, 6:10] <- C
  expect_equal(.study06_apply_blocks(blocks, v), as.numeric(dense %*% v))
  expect_equal(dense[1:5, 6:10], matrix(0, 5, 5))
})

test_that("fixed ridge and posterior equivalence rules are frozen", {
  candidates <- data.frame(shrinkage_weight = c(.001, .01, .05),
    eigen_eta = c(.001 / .999, .01 / .99, .05 / .95), pass = TRUE)
  expect_equal(.study06_select_fixed_ridge(candidates, .01)$eigen_eta,
    .01 / .99)
  exact <- .study06_scalar_equivalence(1, 1, .01, .01, .9, 1.1, .9, 1.1)
  expect_equal(exact$classification, "numerically_equivalent")
  mc <- .study06_scalar_equivalence(1.01, 1, .01, .01, .9, 1.1, .9, 1.1)
  expect_equal(mc$classification, "consistent_with_monte_carlo_error")
  material <- .study06_scalar_equivalence(1.5, 1, .01, .01,
    1.4, 1.6, .9, 1.1)
  expect_equal(material$classification, "material_difference")
})

test_that("prediction and paired metric orientation is explicit", {
  run <- list(status = "ok", architecture = "sparse_homogeneous",
    replicate = 1L,
    method = list(configuration = "bed", native_method = "bayesc"),
    fit = list(bm = matrix(c(1, 0), 2, 1)))
  simulation <- list(phenotype = matrix(c(1, 0), 2, 1,
      dimnames = list(c("a", "b"), "trait1")),
    genetic_values = matrix(c(1, 1), 2, 1,
      dimnames = list(c("a", "b"), "trait1")))
  Z <- matrix(c(1, 0, 2, 0), 2, 2,
    dimnames = list(c("a", "b"), c("m1", "m2")))
  split <- list(test_ids = c("a", "b"))
  m <- .study06_prediction_metrics(run, simulation, Z, split)
  expect_equal(m$value[m$metric == "phenotype_prediction_rmse"], 0)
  table <- expand.grid(architecture = "sparse_homogeneous",
    replicate = 1:5, configuration = c("block_csr",
      "block_eigen_unfiltered"), metric = "x",
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  table$method <- "sbayesc"
  table$value <- ifelse(table$configuration ==
    "block_eigen_unfiltered", 2, 1)
  paired <- .study06_paired_differences(table)
  z <- paired[paired$comparison_id ==
    "unfiltered_block_eigen_minus_block_csr", ]
  expect_equal(z$difference, rep(1, 5))
})

test_that("optional fields and manifest flags are safe", {
  safe_flag <- function(x, name)
    isTRUE(length(x[[name]]) == 1L && !is.na(x[[name]]) &&
      x[[name]])
  expect_false(safe_flag(list(), "optional"))
  expect_false(safe_flag(list(optional = NULL), "optional"))
  expect_true(safe_flag(list(optional = TRUE), "optional"))
})

test_that("Study 06 fallback preserves posterior draw vectors", {
  draws <- c(0.1, 0.2, 0.3)
  expect_identical(draws %||% 0, draws)
  expect_identical(NULL %||% 0, 0)
  expect_identical(numeric() %||% 0, 0)
  expect_identical(NA_real_ %||% 0, 0)
})

test_that("Study 06 sampler-free contract smoke test passes", {
  source(testthat::test_path("..", "..", "studies",
    "06_ld_operator", "contract_smoke_test.R"), local = TRUE)
  expect_true(run_study06_contract_smoke_test())
})

test_that("Study 06 runner contains no Git mutation or network command", {
  path <- testthat::test_path("..", "..", "scripts",
    "run_study06_ld_operator.R")
  if (!file.exists(path)) skip("runner not created yet")
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_false(grepl("git (add|commit|push|tag)", text,
    ignore.case = TRUE))
  expect_false(grepl("install\\.packages|download\\.file|curl|wget",
    text, ignore.case = TRUE))
  expect_false(grepl("requireNamespace\\(['\\\"]tarchetypes|library\\(['\\\"]tarchetypes",
    text, ignore.case = TRUE))
})
})
