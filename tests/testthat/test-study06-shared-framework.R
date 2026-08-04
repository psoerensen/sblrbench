study06_spec <- function() read_benchmark_spec(test_path("..", "..", "studies",
  "06_ld_operator", "spec.R"))

study06_logic <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(test_path("..", "..", "studies", "06_ld_operator",
    "operator-design.R"), environment)
  environment
}

test_that("Study 06 spec and coordinate contract are preserved", {
  spec <- study06_spec()
  expect_invisible(validate_benchmark_spec(spec))
  benchmark <- benchmark_seeds(spec, "benchmark")
  workshop <- benchmark_seeds(spec, "workshop")
  expect_equal(nrow(benchmark), 60L)
  expect_equal(nrow(workshop), 12L)
  expect_identical(unique(benchmark$scenario),
    c("sparse_homogeneous", "sparse_mixture"))
  expect_identical(unique(benchmark$configuration),
    c("bed", "full_csr", "block_csr", "low_rank_full",
      "low_rank_0999", "low_rank_0995"))
  expect_equal(benchmark$simulation_seed[1:6], rep(61100L, 6L))
  expect_equal(benchmark$fit_seed[1:6],
    61100L + 600000L + (1:6) * 1000L)
  expect_equal(benchmark$chain_seeds[[1]],
    benchmark$fit_seed[1] + (1:4) * 101L)
})

test_that("Study 06 operator, prior, block, and retention policies are exact", {
  spec <- study06_spec()
  expect_equal(spec$controls$simulation$h2, 0.30)
  expect_equal(spec$controls$simulation$n_causal, 50L)
  expect_equal(spec$scenarios$sparse_mixture$mixture_prob, c(.60, .30, .10))
  expect_equal(spec$scenarios$sparse_mixture$mixture_var, c(.01, .1, 1))
  expect_equal(spec$controls$bayesr$pi,
    c(.99, rep(.01 / 3, 3L)))
  expect_equal(spec$controls$bayesr$mixture_var, c(0, .01, .1, 1))
  expect_equal(spec$operators$block$size, 1000L)
  expect_equal(spec$operators$eigen$proportions[c("low_rank_0999",
    "low_rank_0995")], c(low_rank_0999 = .999, low_rank_0995 = .995))
  expect_equal(unlist(spec$operators$equivalence_tolerances),
    c(absolute = 2e-3, relative = 2e-6, product_absolute = 2e-2,
      quadratic_absolute = 5e-2, probability = 1e-10))
  controls <- spec$controls$benchmark$recommendations
  homogeneous <- controls[controls$scenario == "sparse_homogeneous", ]
  mixture <- controls[controls$scenario == "sparse_mixture", ]
  expect_equal(homogeneous$nit, rep(250L, 6L))
  expect_equal(homogeneous$nburn, c(250L, 250L, 250L, 250L, 500L, 250L))
  expect_equal(mixture$nit, c(2000L, 2000L, 1000L, 1000L, 1000L, 1000L))
  expect_equal(mixture$nburn, c(250L, 250L, 250L, 250L, 250L, 500L))
})

test_that("supplemental Study 06 window and SBayesR controls are preserved", {
  supplemental <- study06_spec()$supplemental
  expect_equal(supplemental$marker_window$source_rows, 13392:14891)
  expect_equal(supplemental$block_starts,
    c(1L, 251L, 501L, 624L, 751L, 1001L, 1251L))
  expect_equal(supplemental$retained_block_ranks,
    c(248L, 248L, 123L, 127L, 248L, 248L, 248L))
  expect_equal(supplemental$retained_total_rank, 1490L)
  expect_equal(supplemental$simulation_seed, 17002L)
  expect_equal(supplemental$chain_seeds,
    c(150104L, 250104L, 350104L, 450104L))
  expect_equal(unname(unlist(supplemental$controls[c("nburn", "nit", "nthin",
    "nchains", "ncores", "seed")])), c(250, 1000, 1, 4, 4, 50104))
})

test_that("Study 06 block and low-rank fixtures are deterministic", {
  logic <- study06_logic()
  ids <- paste0("m", seq_len(37991L))
  blocks <- logic$study06_blocks(ids, 1000L)
  expect_equal(nrow(blocks), 38L)
  expect_equal(tail(blocks$size, 1L), 991L)
  expect_invisible(logic$study06_validate_blocks(blocks, ids))
  operator <- diag(c(5, 3, 2, 1))
  retained <- logic$study06_retain_eigensystem(operator, .80)
  expect_equal(retained$retained_rank, 3L)
  expect_equal(retained$retained_mass, 10 / 11)
  expect_equal(crossprod(retained$factor), retained$reconstructed)
})

test_that("operator metrics reproduce deterministic known values", {
  reference <- diag(3)
  candidate <- diag(c(1, .5, 1))
  metrics <- operator_matrix_metrics(reference, candidate)
  expect_equal(metrics$maximum_absolute_error, .5)
  expect_equal(metrics$diagonal_error, .5)
  expect_equal(metrics$symmetry_error, 0)
  action <- operator_action_metrics(reference, candidate,
    matrix(c(1, 2, 3), ncol = 1L))
  expect_equal(action$product_maximum_absolute_error, 1)
  expect_equal(action$quadratic_form_maximum_absolute_error, 2)
})

test_that("Study 06 semantic checkpoint identities include operator settings", {
  spec <- study06_spec()
  inputs <- list(study = spec$study, component = "main", scenario = "sparse_mixture",
    replicate = 1L, method = "sparse_mixture__low_rank_0995",
    operator = list(contract = spec$operators$contract,
      block = spec$operators$block, eigen = spec$operators$eigen,
      retained_mass = .995),
    controls = spec$controls$bayesr, seeds = benchmark_seeds(spec)[6, ],
    package_sha = spec$packages$sblr$sha)
  first <- benchmark_semantic_checkpoint_hash(
    benchmark_semantic_checkpoint_identity("study06-main", inputs))
  inputs$operator$retained_mass <- .999
  second <- benchmark_semantic_checkpoint_hash(
    benchmark_semantic_checkpoint_identity("study06-main", inputs))
  expect_false(identical(first, second))
})

test_that("Study 06 reporting helpers return named ggplot objects", {
  comparison <- data.frame(configuration = c("a", "b"),
    maximum_corrected_score_error = c(0, .1))
  eigen <- data.frame(configuration = c("a", "a"), block_start = c(0, 10),
    retained_rank = c(10, 9), retained_mass_fraction = c(1, .995))
  recovery <- data.frame(architecture = "sparse_homogeneous",
    configuration = c("a", "b"), estimand = "heritability", rmse = c(.1, .2))
  expect_s3_class(plot_operator_errors(comparison), "ggplot")
  expect_s3_class(plot_operator_retained_rank(eigen), "ggplot")
  expect_s3_class(plot_operator_spectrum(eigen), "ggplot")
  expect_s3_class(plot_operator_recovery(recovery), "ggplot")
})

test_that("Study 06 validation-only execution cannot dispatch fits", {
  spec <- study06_spec()
  calls <- 0L
  old <- options(sblrbench.ld_operator_runner = function(...) {
    calls <<- calls + 1L
    stop("fit runner must not be called")
  })
  on.exit(options(old), add = TRUE)
  output <- tempfile("study06-validate-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  result <- run_benchmark(spec, output, "workshop", validate_only = TRUE)
  expect_equal(calls, 0L)
  expect_equal(nrow(result$status), 12L)
  expect_true(all(result$status$status == "not_run_validate_only"))
})

test_that("Study 06 workflows preserve capsule-only reports and named plots", {
  root <- test_path("..", "..")
  analysis <- readLines(file.path(root, "studies/06_ld_operator/analysis.R"),
    warn = FALSE)
  expect_true(any(grepl("results <- run_benchmark", analysis, fixed = TRUE)))
  for (name in c("operator_error_plot", "retained_rank_plot",
      "eigenvalue_plot", "recovery_plot", "runtime_plot"))
    expect_true(any(grepl(paste0(name, " <-"), analysis, fixed = TRUE)))
  for (report in c("studies/06_ld_operator/report.qmd",
      "studies/06_ld_operator/sbayesr_ld_robustness/report.qmd")) {
    text <- readLines(file.path(root, report), warn = FALSE)
    expect_false(any(grepl("results/local|readRDS|run_benchmark|stblr_", text)))
  }
})
