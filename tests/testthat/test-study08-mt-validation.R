study08_root <- if (file.exists(file.path("studies", "08_mt_validation",
  "config.R"))) "." else file.path("..", "..")
study08_dir <- file.path(study08_root, "studies", "08_mt_validation")
if (dir.exists(study08_dir)) {
  for (f in c("interface_audit.R", "state_contract.R", "simulation.R", "alignment.R",
    "operators.R", "methods.R", "chain_extraction.R", "diagnostics.R",
    "runtime_scaling.R", "metrics.R", "pilot.R", "promotion.R"))
    source(file.path(study08_dir, f), local = TRUE)
  study08_config <- source(file.path(study08_dir, "config.R"),
    local = TRUE)$value
}

test_that("Study 08 source is available in repository tests", {
  skip_if_not(dir.exists(study08_dir))
  expect_true(file.exists(file.path(study08_dir, "targets.R")))
  readme <- paste(readLines(file.path(study08_dir, "README.md"), warn = FALSE),
    collapse = "\n")
  expect_match(readme, "Status: In development", fixed = TRUE)
  expect_match(readme, "no completed benchmark", ignore.case = TRUE)
})

test_that("two-trait joint states have exact stable semantics", {
  skip_if_not(dir.exists(study08_dir))
  expect_no_error(.study08_validate_state_contract(study08_config))
  states <- .study08_state_models(study08_config$trait_names)
  expect_identical(.study08_state_id(states), 0:3)
  expect_identical(unname(.study08_state_inclusion(0:3)), unname(states))
  p <- rbind(c(.1, .2, .3, .4), c(.25, .25, .25, .25))
  pip <- .study08_trait_pip_from_state_probabilities(p, states)
  expect_equal(pip[, 1L], p[, 2L] + p[, 4L])
  expect_equal(pip[, 2L], p[, 3L] + p[, 4L])
  swapped <- .study08_permute_states(p)
  expect_equal(swapped[, 2L], p[, 3L])
  expect_equal(swapped[, 3L], p[, 2L])
  expect_equal(swapped[, 4L], p[, 4L])
  expect_equal(rowSums(p), rep(1, 2L))
})

test_that("simulation truth and covariance identities are exact", {
  skip_if_not(dir.exists(study08_dir))
  set.seed(71)
  Z <- scale(matrix(rnorm(600 * 120), 600, 120))
  rownames(Z) <- paste0("i", seq_len(nrow(Z)))
  colnames(Z) <- paste0("m", seq_len(ncol(Z)))
  config <- study08_config
  for (a in config$contract_architectures)
    config$simulation$state_counts[[a]] <- pmin(
      config$simulation$state_counts[[a]], c(8L, 8L, 12L))
  for (a in config$contract_architectures) {
    sim <- .study08_simulate(Z, a, 1L, config)
    expect_no_error(.study08_validate_simulation(sim, Z, config))
    expect_equal(stats::cov(sim$genetic_values), sim$truth$cov_g,
      tolerance = 1e-12)
    expect_equal(sim$truth$genetic_correlation,
      sim$truth$cov_g[1L, 2L] / sqrt(prod(diag(sim$truth$cov_g))),
      tolerance = 1e-12)
    expect_equal(sim$truth$cov_e[1L, 2L], 0, tolerance = 1e-10)
    expect_true(all(eigen(sim$truth$cov_g, symmetric = TRUE,
      only.values = TRUE)$values >= -1e-10))
  }
})

test_that("trait marker and sample permutations preserve deterministic values", {
  skip_if_not(dir.exists(study08_dir))
  set.seed(72)
  Z <- scale(matrix(rnorm(80 * 30), 80, 30))
  B <- matrix(rnorm(60), 30, 2); Y <- Z %*% B
  evidence <- .study08_permutation_contract(Z, Y, B, 73L)
  expect_true(all(evidence$passed))
  expect_lte(max(evidence$maximum_absolute_error), 1e-10)
  ids <- paste0("m", seq_len(ncol(Z)))
  expect_identical(.study08_marker_subset(ids, 20L), ids[1:20])
})

test_that("chain extraction derives heritability and genetic correlation drawwise", {
  skip_if_not(dir.exists(study08_dir))
  nit <- 8L; chains <- 4L
  values <- array(NA_real_, c(nit, chains, 6L))
  values[, , 1L] <- 2; values[, , 2L] <- 8
  values[, , 3L] <- 3; values[, , 4L] <- 2
  values[, , 5L] <- 1; values[, , 6L] <- .8
  q <- data.frame(group = c("vgs", "vgs", "ves", "ves", "cov_g",
    "pattern_pi"), trait_index = c(1L, 2L, 1L, 2L, 1L, -1L),
    trait2_index = c(-1L, -1L, -1L, -1L, 2L, -1L),
    pattern_name = c(rep(NA, 5L), "both"), stringsAsFactors = FALSE)
  fit <- list(bm = matrix(0, 5, 2,
    dimnames = list(paste0("m", 1:5), c("trait1", "trait2"))),
    convergence_traces = list(values = values, quantities = q))
  run <- list(status = "ok", fit = fit,
    implementation = list(id = "mt_csr_sbayesc"))
  draws <- .study08_extract_draws(run, "partially_shared", 1L)
  expect_equal(unique(draws$value[draws$estimand == "heritability__trait1"]),
    2 / 5)
  expect_equal(unique(draws$value[draws$estimand == "heritability__trait2"]),
    8 / 10)
  expect_equal(unique(draws$value[draws$estimand == "genetic_correlation"]),
    1 / sqrt(16))
  expect_true(all(table(draws$estimand, draws$chain) == nit))
})

test_that("runtime marker selection enforces all frozen limits", {
  skip_if_not(dir.exists(study08_dir))
  impl <- study08_config$runtime_implementations
  projection <- expand.grid(implementation = impl,
    marker_count = c(2000L, 4000L), stringsAsFactors = FALSE)
  projection$feasible <- TRUE
  projection$projected_10_fit_seconds <- 100
  selected <- .study08_select_marker_count(projection, study08_config)
  expect_identical(selected$marker_count, 4000L)
  projection$feasible[projection$marker_count == 4000L &
    projection$implementation == impl[[1L]]] <- FALSE
  expect_identical(.study08_select_marker_count(projection,
    study08_config)$marker_count, 2000L)
})

test_that("optional manifest flags are report safe", {
  skip_if_not(dir.exists(study08_dir))
  expect_false(.study08_optional_flag(NULL))
  expect_false(.study08_optional_flag(logical()))
  expect_false(.study08_optional_flag(NA))
  expect_true(.study08_optional_flag(TRUE))
  expect_error(.study08_validate_capsule(tempfile(), "contract"),
    "incomplete")
})

test_that("interfaces and policies are frozen without unsupported scope", {
  skip_if_not(dir.exists(study08_dir))
  audit <- .study08_interface_audit(study08_config)
  expect_setequal(audit$interface,
    c("mtblr_bed", "mtblr_csr", "mtblr_block_eigen"))
  expect_true(all(!audit$sampled_maf_s))
  expect_identical(study08_config$trait_count, 2L)
  expect_identical(study08_config$implementations,
    c("mt_bed_bayesc", "mt_csr_sbayesc", "mt_block_eigen_sbayesc"))
})

test_that("reconstructed-dense MT block-eigen execution is explicitly paused", {
  skip_if_not(dir.exists(study08_dir))
  expect_no_error(.study08_assert_execution_allowed("mt_bed_bayesc"))
  expect_no_error(.study08_assert_execution_allowed("mt_csr_sbayesc"))
  expect_no_error(.study08_assert_execution_allowed("mt_block_csr_sbayesc"))
  expect_no_error(.study08_assert_phase_allowed("contract"))
  expect_error(.study08_assert_phase_allowed("runtime"),
    "Study 08 MT block-eigen execution is paused")
  expect_error(.study08_assert_phase_allowed("all"),
    "Study 08 MT block-eigen execution is paused")
  expect_error(
    .study08_assert_execution_allowed("mt_block_eigen_sbayesc"),
    paste(
      "Study 08 MT block-eigen execution is paused",
      "current mtblr_block_eigen\\(\\) backend is reconstructed dense",
      sep = ".*"
    )
  )
  expect_identical(study08_config$execution_status,
    "paused_pending_retained_low_rank_mt_operator")
  expect_false("mt_block_eigen_sbayesc" %in%
    study08_config$contract_implementations)
})
