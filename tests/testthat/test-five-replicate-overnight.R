root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
study_sources_available <- file.exists(file.path(root, "studies", "five_replicate_helpers.R"))
if (study_sources_available)
  source(file.path(root, "studies", "five_replicate_helpers.R"), local = TRUE)

test_that("active pipelines do not require tarchetypes", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  files <- c("_targets.R", "studies/02_prediction/targets.R",
    "studies/03_parameter_estimation/targets.R",
    "studies/04_convergence/targets.R",
    "studies/04_convergence/validation_targets.R",
    "scripts/run_five_replicate_overnight.R")
  text <- unlist(lapply(file.path(root, files), readLines, warn = FALSE))
  expect_false(any(grepl("tarchetypes::|library\\(tarchetypes\\)|require\\(tarchetypes\\)", text)))
})

test_that("frozen recommendations drive exact four-chain settings", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  rec <- .five_replicate_recommendations(file.path(root,
    "results/reference/04_convergence/st-multichain-convergence-development-v1/method_recommendations.csv"))
  expect_equal(rec$recommended_nburn, rep(250L, 4L))
  expect_equal(rec$recommended_post_burnin_draws, c(250L, 1000L, 250L, 1000L))
  expect_equal(rec$recommended_nchains, rep(4L, 4L))
  expect_equal(rec$recommended_nthin, rep(1L, 4L))
})

test_that("Study 02 five-replicate grid and seeds are order invariant", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  source(file.path(root, "studies/02_prediction/pilot.R"), local = TRUE)
  cfg <- source(file.path(root, "studies/02_prediction/config.R"), local = TRUE)$value
  cfg$profile <- "five_replicate_development"
  specs <- .study02_replicate_specs(cfg, "", "")
  methods <- .study02_method_specs(cfg)
  expect_length(specs, 10L)
  expect_equal(sort(unique(vapply(specs, `[[`, integer(1), "replicate"))), 1:5)
  grid <- expand.grid(spec = seq_along(specs), method = seq_along(methods))
  expect_equal(nrow(grid), 40L)
  seed <- function(s, m) cfg$mcmc$seed_offset +
    match(s$architecture, names(cfg$simulation$architectures)) * 10000L +
    s$replicate * 100L + m$method_index
  forward <- sort(vapply(seq_len(nrow(grid)), function(i)
    seed(specs[[grid$spec[i]]], methods[[grid$method[i]]]), integer(1)))
  reverse <- sort(vapply(rev(seq_len(nrow(grid))), function(i)
    seed(specs[[grid$spec[i]]], methods[[grid$method[i]]]), integer(1)))
  expect_identical(forward, reverse)
  expect_equal(anyDuplicated(.five_replicate_chain_seeds(forward[1L])), 0L)
})

test_that("Study 02 summaries retain failures and use successful-replicate SD", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  source(file.path(root, "studies/02_prediction/promotion.R"), local = TRUE)
  x <- data.frame(architecture = "sparse_homogeneous", method = "st_bed_bayesc",
    metric = "prediction_correlation", replicate = 1:5,
    value = c(1, 2, 3, 4, NA), status = c(rep("ok", 4), "failed"))
  z <- .study02_benchmark_summary(x)
  expect_equal(z$replicate_count, 5L)
  expect_equal(z$successful_replicates, 4L)
  expect_equal(z$sd, stats::sd(1:4))
})

test_that("Study 03 validates and pools four chains after draw-wise transforms", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  for (f in c("estimands.R", "methods.R", "metrics.R"))
    source(file.path(root, "studies/03_parameter_estimation", f), local = TRUE)
  registry <- .study03_estimand_registry()
  nit <- 3L; nburn <- 2L
  a <- array(NA_real_, dim = c(nit, 4L, 3L))
  for (ch in 1:4) {
    a[, ch, 1] <- c(1, 2, 4) + ch
    a[, ch, 2] <- c(2, 3, 5) + ch
    a[, ch, 3] <- c(3, 4, 6) + ch
  }
  chains <- lapply(1:4, function(ch) list(pi_trace = matrix(c(99, 99, .1, .2, .4), ncol = 1)))
  fit <- list(input = list(nit = nit, nburn = nburn), chains = chains,
    convergence_traces = list(values = a,
      quantities = data.frame(group = c("vbs", "vgs", "ves"))))
  draws <- .study03_extract_multichain_draws(fit, "st_bed_bayesc", registry, 10L)
  expect_equal(unique(table(draws$estimand_id, draws$chain)), nit)
  total <- draws[draws$estimand_id == "total_marker_effect_variance", ]$value
  vbs <- draws[draws$estimand_id == "effect_variance", ]$value
  expect_equal(total[1:3], vbs[1:3] * c(.1, .2, .4) * 10)
  summary <- .study03_summarise_draws(draws)
  expect_true(all(summary$chain_count == 4L))
  expect_true(all(summary$draws_per_chain == nit))
})

test_that("Study 04 validation grid contains 20 fits and 80 unique chains", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  for (f in c("pilot.R", "methods.R"))
    source(file.path(root, "studies/04_convergence", f), local = TRUE)
  cfg <- source(file.path(root, "studies/04_convergence/config.R"), local = TRUE)$value
  cfg$profile <- "five_replicate_validation"
  specs <- .study04_specs(cfg)
  expect_length(specs, 20L)
  seeds <- do.call(rbind, lapply(specs, function(x) data.frame(
    architecture = x$architecture, replicate = x$replicate, method = x$method,
    chain = 1:4, seed = .study04_chain_seeds(x$architecture, x$method, cfg, x$replicate))))
  expect_equal(nrow(seeds), 80L)
  expect_equal(anyDuplicated(seeds$seed), 0L)
  counts <- table(paste(seeds$architecture, seeds$replicate, seeds$method, sep = "|"))
  expect_true(all(counts == 4L))
})

test_that("orchestrator validation mode exits before sampling and contains no Git mutation", {
  skip_if_not(study_sources_available, "repository-only study files are excluded from package builds")
  text <- readLines(file.path(root, "scripts/run_five_replicate_overnight.R"), warn = FALSE)
  validate_line <- grep("if \\(validate_only\\) quit", text)
  study_line <- grep("phase_ok\\[\"study02\"\\]", text)
  expect_length(validate_line, 1L)
  expect_true(validate_line < study_line)
  expect_true(any(grepl("record\\(phase", text)))
  expect_false(any(grepl("git (add|commit|push|tag)|c\\(\"(add|commit|push|tag)\"", text)))
})
