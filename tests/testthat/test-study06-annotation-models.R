study06_root <- test_path("..", "..")
study06_spec <- read_benchmark_spec(file.path(study06_root, "studies",
  "06_annotation_models", "spec.R"))
study06_logic <- new.env(parent = globalenv())
sys.source(file.path(study06_root, "studies", "06_annotation_models",
  "annotation-design.R"), envir = study06_logic)

study06_output <- function(label) file.path(tempdir(),
  study06_spec$study_version, label)

test_that("Study 06 v1 current-stop evidence is byte-for-byte preserved", {
  expected <- c(
    annotation_design_summary.csv = "8efdbb05f149527896afb858736c325c4e296b77508cafa943826a74701a56ea",
    benchmark_manifest.json = "cb0a014d52eb8798b6374196899abd1dd4dd5f85432d34d43fdec8b3e093fc66",
    candidate_settings.csv = "cacb65d1ef146b20290981417464ae87cd8a89bed5059e1292c793dd7a8aac01",
    checksums.csv = "d9bbb5274ef260581fee5b0a7e4f3da0f7220e78de8ae8f22a7b8b2d05288e64",
    computational_summary.csv = "d4b40171cdcde2bb5a81f42ec18ea430ef5f838cea574facb7c17d7fb8b1efdf",
    config.R = "11550e5b35c4ce01f4499689730ccfde3136920d56f9e5f05f006bfc3a9f0f00",
    contract_smoke_test.R = "09cabb5e5304852f33a528221d025c91384dab96669cc2b7b706863661c17b7e",
    convergence_diagnostics.csv = "f872c5bb184a072daf7d9a83961cfcdf0da746ebe5963c0c9bb06e7318829c00",
    example_data_manifest.csv = "06b24f8659bb86c84ebd2c288cbe33e82df1711de9d6fd358e9b96de0510a897",
    fit_status.csv = "d4b40171cdcde2bb5a81f42ec18ea430ef5f838cea574facb7c17d7fb8b1efdf",
    interface_audit_sources.csv = "4399351f9e458992f825766402065dc94a794d3a5579308124487d663613f05f",
    method_recommendations.csv = "cc9ef74cfbe919f48f96a7bb563f347a74815d1bb35a1769fee675b7b5b4d5e9",
    README.md = "f22f1712a9edc041e95d516c3944a2f55202334886d0a94e04a9adbdfba94bbe",
    reproduce.R = "8b498e44ad1a3663aac3bab4b485c42e5b66da66bb27a6db922593fbe515afef",
    scalar_chain_draws.csv = "d96e65518cf692e694571d42b25506f60c0f92c45fd0a4a9757677a79b7b6bbf",
    seed_registry.csv = "9548ca6f1d6fce06dda39ea746746194332bb4fbc55bef16a44922ccc21a196f",
    session_info.txt = "5c04881a2f59d736448437c86d0b5ee14e4f29e67f02dabadc0ef835d60f57da",
    source_files.csv = "015b3df1005b47e3becce0548cfdaf616e53132d5196def9ad1e3e1fe8b69d1a",
    true_alpha.csv = "b6768e971f8fa75ff7f74df2925b591301cb37cd84ad055e4969fe74e919f7bd")
  capsule <- file.path(study06_root, study06_spec$frozen_capsule$current_stop)
  files <- list.files(capsule, recursive = TRUE, all.files = FALSE)
  expect_setequal(files, names(expected))
  observed <- vapply(file.path(capsule, names(expected)), digest::digest,
    character(1), algo = "sha256", file = TRUE, serialize = FALSE)
  names(observed) <- names(expected)
  expect_identical(observed, expected)
})

test_that("Study 06 v2 spec and four-route registry are prespecified", {
  expect_invisible(validate_benchmark_spec(study06_spec))
  expect_identical(study06_spec$study_version,
    "v2_identifiable_qualification")
  expect_identical(names(study06_spec$versioning$profiles), c(
    "v1_sparse_qualification_failed", "v1_sparse_stress",
    "v2_identifiable_qualification"))
  expect_identical(names(study06_spec$methods), c("st_bed_bayesr",
    "st_bed_bayesrc", "st_block_eigen_sbayesr",
    "st_block_eigen_sbayesrc"))
  block <- study06_spec$methods[3:4]
  expect_true(all(vapply(block, function(x)
    identical(x$interface, "stblr_block_eigen") &&
      identical(x$operator_representation, "low_rank") &&
      identical(x$eigen_prop, .995), logical(1))))
  expect_null(study06_spec$controls$priors$intercept_flat)
  text <- paste(capture.output(str(study06_spec$methods)), collapse = " ")
  expect_false(grepl("CSR|pair|block.allocation", text, ignore.case = TRUE))
})

test_that("v2 marker selection is deterministic, separated, and exhaustive", {
  ids <- sprintf("m%05d", seq_len(37991L))
  positions <- cumsum(rep(1000, length(ids)))
  one <- study06_logic$select_study06_v2_blocks(ids, positions, study06_spec)
  two <- study06_logic$select_study06_v2_blocks(ids, positions, study06_spec)
  expect_identical(one, two)
  expect_equal(nrow(one), 1500L)
  expect_identical(as.integer(table(one$block_id)), rep(100L, 15L))
  expect_equal(anyDuplicated(one$marker_id), 0L)
  expect_true(all(diff(one$qc_index) > 0L))
  expect_identical(one$within_block_index, rep(seq_len(100L), 15L))
  expect_gt(min(diff(tapply(one$position_bp, one$block_id, min))), 100L * 1000L)
})

test_that("block audit checks rank, conditioning, and method marker equality", {
  ids <- sprintf("m%05d", seq_len(37991L))
  panel <- study06_logic$select_study06_v2_blocks(ids, seq_along(ids),
    study06_spec)
  set.seed(41)
  x <- matrix(rnorm(220L * 1500L), 220L, 1500L,
    dimnames = list(paste0("s", 1:220), panel$marker_id))
  audit <- study06_logic$audit_study06_v2_blocks(panel, x,
    list(BED = panel$marker_id, block_eigen = panel$marker_id))
  expect_equal(nrow(audit$blocks), 15L)
  expect_true(all(audit$blocks$marker_count == 100L))
  expect_true(all(audit$blocks$within_block_rank == 100L))
  expect_true(all(audit$blocks$minimum_positive_eigenvalue > 0))
  expect_true(all(audit$blocks$retained_positive_mass >= .995))
  expect_true(audit$summary$method_marker_identity_and_order_equal)
  expect_true(is.finite(audit$summary$maximum_cross_block_absolute_correlation))
})

test_that("split and annotation designs are deterministic and non-oracle", {
  sample_ids <- sprintf("s%04d", 1:2000)
  split1 <- make_prediction_split(sample_ids, .70, 3101L)
  split2 <- make_prediction_split(sample_ids, .70, 3101L)
  expect_identical(split1, split2)
  expect_length(split1$train_ids, 1400L)
  expect_length(split1$test_ids, 600L)
  expect_length(intersect(split1$train_ids, split1$test_ids), 0L)
  ids <- sprintf("m%04d", 1:1500)
  A1 <- study06_logic$construct_annotation_design(ids, study06_spec)
  A2 <- study06_logic$construct_annotation_design(ids, study06_spec)
  expect_identical(A1, A2)
  expect_true(all(A1[, "Intercept"] == 1))
  expect_equal(sum(A1[, "enriched_binary"]), 225L)
  expect_equal(qr(A1)$rank, 4L)
  expect_equal(A1[, "continuous_signal"], A2[, "continuous_signal"])
})

test_that("informative and uninformative marginals match at the v2 targets", {
  ids <- sprintf("m%04d", 1:1500)
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  truth <- study06_logic$construct_annotation_truth(A, study06_spec)
  pi_i <- study06_logic$annotation_marker_probabilities(A,
    truth$informative_annotations, c(0, .01, .1, 1))
  pi_u <- study06_logic$annotation_marker_probabilities(A,
    truth$uninformative_annotations, c(0, .01, .1, 1))
  expect_equal(unname(colMeans(pi_i)), c(.88, .06, .036, .024), tolerance = 1e-6)
  expect_equal(colMeans(pi_i), colMeans(pi_u), tolerance = 1e-10)
  expect_equal(sum(1 - pi_i[, 1L]), 180, tolerance = .01)
  expect_equal(unname(colSums(pi_i)[-1L]), c(90, 54, 36), tolerance = .01)
  expect_true(all(truth$uninformative_annotations[-1L, ] == 0))
  expect_true(all(truth$informative_annotations["null_annotation", ] == 0))
  enriched <- A[, "enriched_binary"] == 1
  share <- sum((1 - pi_i[, 1L])[enriched]) / sum(1 - pi_i[, 1L])
  expect_gte(share, .50)
  expect_lte(share, .60)
  expect_gt(mean(1 - pi_i[enriched, 1L]), mean(1 - pi_i[!enriched, 1L]))
})

study06_simulation_fixture <- function(scenario = "informative_annotations",
                                       component_seed = NULL) {
  ids <- sprintf("m%04d", 1:1500)
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  truth <- study06_logic$construct_annotation_truth(A, study06_spec)
  coordinate <- benchmark_annotation_seeds(study06_spec, "benchmark",
    "qualification")
  coordinate <- coordinate[coordinate$scenario == scenario, , drop = FALSE][1L, ]
  if (!is.null(component_seed)) coordinate$component_seed <- component_seed
  set.seed(91)
  W <- matrix(rnorm(200L * 1500L), 200L, 1500L,
    dimnames = list(sprintf("s%03d", 1:200), ids))
  simulation <- study06_logic$simulate_annotation_architecture(
    as.list(coordinate), W, 1:140, A, truth, study06_spec)
  list(simulation = simulation, A = A, truth = truth, coordinate = coordinate)
}

test_that("truth identifiability, component minima, and heritability are enforced", {
  for (scenario in names(study06_spec$scenarios)) {
    fixture <- study06_simulation_fixture(scenario)
    counts <- fixture$simulation$extras$component_counts
    expect_gte(sum(counts[-1L]), 160L)
    expect_lte(sum(counts[-1L]), 200L)
    expect_true(all(counts[-1L] >= c(60L, 30L, 20L)))
    sticks <- fixture$simulation$extras$stick_counts
    expect_true(all(sticks$eligible_count >= c(500L, 120L, 50L)))
    expect_true(all(pmin(sticks$outcome_zero_count,
      sticks$outcome_one_count) >= c(100L, 25L, 20L)))
    expect_true(all(sticks$enriched_outcome_cell_minimum > 0L))
    expect_true(is.finite(fixture$simulation$extras$realized_heritability))
    expect_equal(fixture$simulation$extras$realized_heritability, .5,
      tolerance = .02)
    marker <- fixture$simulation$extras$marker_truth
    expect_true(any(marker$true_nonnull & marker$enriched_binary == 0L))
    expect_true(any(!marker$true_nonnull & marker$enriched_binary == 1L))
  }
})

test_that("replicate seeds reproduce and change only replicate truth", {
  one <- study06_simulation_fixture()
  again <- study06_simulation_fixture()
  changed <- study06_simulation_fixture(component_seed =
    one$coordinate$component_seed + 1L)
  expect_identical(one$simulation$extras$marker_truth,
    again$simulation$extras$marker_truth)
  expect_identical(one$simulation$truth$phenotypes,
    again$simulation$truth$phenotypes)
  expect_false(identical(one$simulation$extras$marker_truth$component_index,
    changed$simulation$extras$marker_truth$component_index))
  expect_identical(one$simulation$extras$study_version,
    study06_spec$study_version)
})

test_that("diagnostic isolation modes reflect public API availability", {
  coordinate <- benchmark_annotation_seeds(study06_spec, "benchmark",
    "qualification")[1L, , drop = FALSE]
  standard <- annotation_method_controls(study06_spec, coordinate,
    mode = "qualification")
  fixed <- annotation_method_controls(study06_spec, coordinate,
    mode = "qualification",
    diagnostic_mode = "fixed_true_annotation_coefficients")
  expect_false("intercept_flat" %in% names(standard))
  expect_true(standard$updateAlpha)
  expect_false(fixed$updateAlpha)
  expect_error(annotation_method_controls(study06_spec, coordinate,
    mode = "qualification",
    diagnostic_mode = "fixed_true_marker_component_probabilities"),
    "No current public")
})

test_that("block-eigen SBayesRC uses its supported public gamma control", {
  controls <- list(mixture_var = c(0, .01, .1, 1), h2 = .5)
  block <- sblrbench:::.annotation_public_api_controls(
    study06_spec$methods$st_block_eigen_sbayesrc, controls)
  bed <- sblrbench:::.annotation_public_api_controls(
    study06_spec$methods$st_bed_bayesrc, controls)
  expect_equal(block$gamma, controls$mixture_var)
  expect_false("mixture_var" %in% names(block))
  expect_equal(bed$mixture_var, controls$mixture_var)
  expect_false("gamma" %in% names(bed))
})

test_that("v2 outputs cannot collide with v1 and qualification is opt-in", {
  expect_error(run_benchmark(study06_spec,
    file.path(study06_root, study06_spec$frozen_capsule$current_stop),
    validate_only = TRUE), "cannot overlap v1")
  expect_error(run_benchmark(study06_spec, tempfile("unversioned-study06-"),
    validate_only = TRUE), "versioned namespace")
  old <- options(sblrbench.annotation_fit_dispatch = function(...)
    stop("qualification sampler called"))
  on.exit(options(old), add = TRUE)
  out <- run_benchmark(study06_spec, study06_output("validate-only"),
    profile = "benchmark", validate_only = TRUE)
  expect_equal(nrow(out$qualification_coordinates), 4L)
  expect_equal(nrow(out$coordinate_grid), 40L)
  expect_identical(out$status$status, rep("not_run_validate_only", 40L))
  expect_true(all(out$status$study_version == study06_spec$study_version))
})

test_that("qualification decisions and checkpoints carry v2 identity", {
  schema <- annotation_qualification_artifact_schema()
  expect_identical(schema$schema, "sblrbench-annotation-qualification-v2")
  expect_true("study_version" %in% schema$required_fields)
  coordinate <- benchmark_annotation_seeds(study06_spec, "benchmark",
    "qualification")[1L, , drop = FALSE]
  controls <- annotation_method_controls(study06_spec, coordinate,
    mode = "qualification")
  ids <- paste0("m", 1:12)
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  data <- list(split = list(train_ids = paste0("s", 1:7),
    test_ids = paste0("s", 8:10)), markers = list(marker_ids = ids))
  bundle <- list(annotations = A, simulation = list(truth = list(
    phenotypes = matrix(seq_len(10), ncol = 1),
    effects = matrix(seq_len(12), ncol = 1))))
  identity <- sblrbench:::.annotation_checkpoint_identities(study06_spec,
    coordinate, controls, data, bundle, "qualification")
  expect_identical(identity$checkpoint$scientific_inputs$study_version,
    study06_spec$study_version)
  expect_identical(identity$checkpoint$scientific_inputs$operator_settings$eigen_prop,
    .995)
})

annotation_trace_fixture <- function(A) {
  annotations <- colnames(A)
  sticks <- paste0("component_", 0:2, "_stick")
  alpha_q <- expand.grid(annotation_name = annotations, stick_name = sticks,
    stringsAsFactors = FALSE)
  q <- rbind(data.frame(parameter_name = "alpha", alpha_q),
    data.frame(parameter_name = "sigmaSqAlpha", annotation_name = annotations,
      stick_name = "", stringsAsFactors = FALSE))
  values <- array(0, dim = c(3L, 4L, nrow(q)))
  values[, , seq_len(nrow(alpha_q))] <- rep(seq(-.4, .4,
    length.out = nrow(alpha_q)), each = 12L)
  values[, , (nrow(alpha_q) + 1L):nrow(q)] <- 1
  list(convergence_traces = list(values = values, quantities = q),
    chains = as.list(seq_len(4L)), alpha = matrix(99, 4, 3),
    alpha_final = matrix(-99, 4, 3), annotation_prior = matrix(1, nrow(A), 4))
}

test_that("true retained traces remain required for draw-wise priors", {
  ids <- paste0("m", seq_len(20L))
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  fixture <- annotation_trace_fixture(A)
  traces <- extract_annotation_coefficient_traces(fixture, expected_chains = 4L)
  prior <- summarise_drawwise_annotation_prior(traces, A, c(0, .01, .1, 1))
  expect_identical(prior$status, "ok")
  expect_equal(nrow(prior$draws), 12L)
  missing <- fixture
  missing$convergence_traces <- NULL
  unavailable <- extract_annotation_coefficient_traces(missing)
  expect_identical(unavailable$status, "unavailable")
  expect_match(unavailable$reason, "not substitutes")
})

test_that("qualification parameter truth uses the simulation bundle", {
  simulation <- bench_fixture()
  simulation$data$train_ids <- simulation$data$sample_ids[1:3]
  quantities <- data.frame(group = c("vgs", "ves"))
  values <- array(NA_real_, dim = c(10L, 4L, 2L))
  values[, , 1L] <- 2
  values[, , 2L] <- 1
  result <- new_sblrbench_result("fixture", native_fit = list(
    convergence_traces = list(values = values, quantities = quantities)))
  estimates <- sblrbench:::.annotation_parameter_estimates(result,
    list(scenario = "informative_annotations", replicate = 1L,
      method = "st_bed_bayesrc"),
    list(simulation = simulation))
  expect_equal(estimates$parameter,
    c("genetic_variance", "residual_variance", "heritability"))
  expected_vg <- var(simulation$truth$genetic_values[1:3, 1L])
  expected_ve <- var(simulation$truth$residuals[1:3, 1L])
  expect_equal(estimates$truth,
    c(expected_vg, expected_ve, expected_vg / (expected_vg + expected_ve)))
})

test_that("qualification scientific and route gates are deterministic", {
  ids <- paste0("m", 1:20)
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  fixture <- annotation_trace_fixture(A)
  traces <- extract_annotation_coefficient_traces(fixture, expected_chains = 4L)
  traces$value[traces$parameter == "alpha" &
    traces$annotation %in% c("enriched_binary", "continuous_signal")] <- 1
  null <- traces$parameter == "alpha" & traces$annotation == "null_annotation"
  traces$value[null] <- rep(c(-.2, .2), length.out = sum(null))
  prior <- summarise_drawwise_annotation_prior(traces, A,
    c(0, .01, .1, 1), retain_marker_summary = FALSE)
  checks <- sblrbench:::.annotation_scientific_checks(traces, prior,
    "informative_annotations", study06_spec)
  expect_true(all(checks$pass))
  marker <- data.frame(marker_id = ids,
    posterior_mean_effect = seq(-1, 1, length.out = length(ids)))
  record <- function(method, shift = 0) list(
    scenario = "informative_annotations", method = method,
    marker = transform(marker, posterior_mean_effect =
      posterior_mean_effect + shift),
    parameters = data.frame(parameter = "heritability",
      posterior_mean = .5 + shift),
    prediction = data.frame(metric = "genetic_value_correlation",
      value = .8 + shift))
  records <- list(record("st_bed_bayesrc"),
    record("st_block_eigen_sbayesrc", .01),
    within(record("st_bed_bayesrc"), scenario <- "uninformative_annotations"),
    within(record("st_block_eigen_sbayesrc", .01),
      scenario <- "uninformative_annotations"))
  route <- sblrbench:::.annotation_route_checks(records, study06_spec)
  expect_equal(nrow(route), 6L)
  expect_true(all(route$pass))
})

test_that("paired power isolation registry and shuffle are deterministic", {
  environment <- new.env(parent = globalenv())
  sys.source(file.path(study06_root, "studies", "06_annotation_models",
    "power-isolation.R"), envir = environment)
  registry <- environment$study06_power_registry(study06_spec)
  expect_equal(nrow(registry), 8L)
  expect_identical(as.integer(table(registry$route)), c(4L, 4L))
  expect_identical(as.integer(table(registry$condition)), rep(2L, 4L))
  ids <- sprintf("marker_%04d", 1:1500)
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  first <- environment$study06_shuffle_annotations(A, 6201L)
  second <- environment$study06_shuffle_annotations(A, 6201L)
  expect_identical(first, second)
  expect_true(all(first$annotations[, "Intercept"] == 1))
  expect_lte(max(abs(stats::cor(first$annotations[, -1L]) -
    stats::cor(A[, -1L]))), 1e-14)
})

test_that("paired power controls preserve qualification history and requested traces", {
  environment <- new.env(parent = globalenv())
  sys.source(file.path(study06_root, "studies", "06_annotation_models",
    "power-isolation.R"), envir = environment)
  ids <- sprintf("marker_%04d", 1:1500)
  controls <- environment$study06_power_controls(study06_spec, ids,
    701020L, c(701121L, 701222L, 701323L, 701424L), TRUE, FALSE)
  expect_identical(controls$nit, 9000L)
  expect_identical(controls$nburn, 0L)
  expect_identical(controls$nchains, 4L)
  expect_identical(controls$chain_seeds,
    c(701121L, 701222L, 701323L, 701424L))
  expect_identical(controls$convergence_control$selected_markers, ids)
  expect_identical(controls$convergence_control$selected_marker_quantities,
    "component")
  expect_false(controls$updateAlpha)
  expect_identical(controls$diagnostic_mode,
    "fixed_true_annotation_coefficients")
  unavailable <- environment$study06_power_controls(study06_spec,
    character(), 701020L, c(701121L, 701222L, 701323L, 701424L),
    FALSE, FALSE)
  expect_null(unavailable$convergence_control$selected_markers)
  expect_null(unavailable$convergence_control$selected_marker_quantities)
})

test_that("paired power trace fallback is deterministic and retains every causal marker", {
  environment <- new.env(parent = globalenv())
  sys.source(file.path(study06_root, "studies", "06_annotation_models",
    "power-isolation.R"), envir = environment)
  truth <- data.frame(marker_id = sprintf("m%04d", 1:1500),
    true_nonnull = seq_len(1500) <= 171)
  first <- environment$study06_trace_marker_set(truth)
  second <- environment$study06_trace_marker_set(truth)
  expect_identical(first, second)
  expect_equal(first$marker_count, 496L)
  expect_equal(first$causal_count, 171L)
  expect_true(all(truth$marker_id[truth$true_nonnull] %in% first$marker_ids))
  expect_false(first$complete_genomewide_occupancy)
  expect_lte(first$estimated_extended_gib, 1)
})

test_that("paired power convergence traces include expected active count", {
  environment <- new.env(parent = globalenv())
  sys.source(file.path(study06_root, "studies", "06_annotation_models",
    "power-isolation.R"), envir = environment)
  ids <- paste0("m", seq_len(20L))
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  groups <- c("vbs", "vgs", "ves", rep("component_pi", 4L))
  quantities <- data.frame(group = groups,
    component_name = c(rep(NA_character_, 3L), paste0("component_", 0:3)))
  values <- array(1, dim = c(5L, 4L, length(groups)))
  values[, , 4L] <- .88
  values[, , 5L] <- .06
  values[, , 6L] <- .04
  values[, , 7L] <- .02
  result <- new_sblrbench_result("fixture", native_fit = list(
    convergence_traces = list(values = values, quantities = quantities)))
  draws <- environment$study06_power_required_traces(result, NULL, NULL, A,
    data.frame(marker_id = ids, true_nonnull = FALSE),
    data.frame(fit_id = "baseline--bed", update_alpha = FALSE), study06_spec)
  expected <- draws[draws$quantity == "expected_active_count", ]
  expect_equal(nrow(expected), 20L)
  expect_equal(expected$value, rep(2.4, 20L))
})
