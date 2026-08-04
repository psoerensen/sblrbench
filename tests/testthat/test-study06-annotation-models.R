study06_root <- test_path("..", "..")
study06_spec <- read_benchmark_spec(file.path(study06_root, "studies",
  "06_annotation_models", "spec.R"))
study06_logic <- new.env(parent = globalenv())
sys.source(file.path(study06_root, "studies", "06_annotation_models",
  "annotation-design.R"), envir = study06_logic)

test_that("Study 06 spec preserves the audited grids and qualification gate", {
  expect_invisible(validate_benchmark_spec(study06_spec))
  final <- benchmark_annotation_seeds(study06_spec, "benchmark", "final")
  qualification <- benchmark_annotation_seeds(study06_spec, "benchmark",
    "qualification")
  expect_equal(nrow(final), 40L)
  expect_equal(nrow(qualification), 4L)
  expect_equal(sum(lengths(final$chain_seeds)), 160L)
  expect_identical(qualification[c("scenario", "replicate", "method")],
    study06_spec$qualification$entries)
  expect_equal(anyDuplicated(final$fit_seed), 0L)
  expect_equal(anyDuplicated(unlist(final$chain_seeds)), 0L)
})

test_that("invalid annotation specs fail clearly", {
  broken <- study06_spec
  broken$annotation_design$columns <- rev(broken$annotation_design$columns)
  expect_error(validate_benchmark_spec(broken), "exact order")
  broken <- study06_spec
  broken$qualification$entries <- broken$qualification$entries[-1L, ]
  expect_error(validate_benchmark_spec(broken), "four audited")
  broken <- study06_spec
  broken$qualification$thresholds$rhat <- 1.05
  expect_error(validate_benchmark_spec(broken), "thresholds changed")
  broken <- study06_spec
  broken$annotation_design <- NULL
  expect_error(validate_benchmark_spec(broken), "missing required")
})

test_that("annotation construction is deterministic and aligned", {
  ids <- sprintf("m%05d", seq_len(study06_spec$data$expected_marker_count))
  A1 <- study06_logic$construct_annotation_design(ids, study06_spec)
  A2 <- study06_logic$construct_annotation_design(ids, study06_spec)
  expect_identical(A1, A2)
  expect_identical(rownames(A1), ids)
  expect_identical(colnames(A1), study06_spec$annotation_design$columns)
  expect_true(all(A1[, "Intercept"] == 1))
  expect_equal(sum(A1[, "enriched_binary"]), 3799L)
  expect_equal(colMeans(A1[, c("continuous_signal", "null_annotation")]),
    c(continuous_signal = 0, null_annotation = 0), tolerance = 1e-12)
  expect_equal(apply(A1[, c("continuous_signal", "null_annotation")], 2, sd),
    c(continuous_signal = 1, null_annotation = 1), tolerance = 1e-12)
  truth <- study06_logic$construct_annotation_truth(A1, study06_spec)
  expect_true(all(truth$uninformative_annotations[-1L, ] == 0))
  informative <- study06_logic$annotation_marker_probabilities(A1,
    truth$informative_annotations,
    study06_spec$controls$simulation$mixture_var)
  expect_gt(mean(1 - informative[A1[, "enriched_binary"] == 1, 1]),
    mean(1 - informative[A1[, "enriched_binary"] == 0, 1]))
})

test_that("Study 06 controls and semantic identities preserve annotation inputs", {
  coordinate <- benchmark_annotation_seeds(study06_spec, "benchmark",
    "qualification")[1L, , drop = FALSE]
  controls <- annotation_method_controls(study06_spec, coordinate,
    "benchmark", "qualification")
  expect_identical(controls$nit, 9000L)
  expect_identical(controls$nburn, 0L)
  expect_identical(controls$nchains, 4L)
  expect_equal(controls$mixture_var, c(0, .01, .1, 1))
  ids <- paste0("m", 1:12)
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  data <- list(split = list(train_ids = paste0("s", 1:7),
    test_ids = paste0("s", 8:10)), markers = list(marker_ids = ids))
  bundle <- list(annotations = A, simulation = list(truth = list(
    phenotypes = matrix(seq_len(10), ncol = 1),
    effects = matrix(seq_len(12), ncol = 1))))
  first <- sblrbench:::.annotation_checkpoint_identities(study06_spec,
    coordinate, controls, data, bundle, "qualification")
  bundle$annotations <- A[, c(1, 3, 2, 4), drop = FALSE]
  reordered <- sblrbench:::.annotation_checkpoint_identities(study06_spec,
    coordinate, controls, data, bundle, "qualification")
  expect_false(identical(first$checkpoint_hash, reordered$checkpoint_hash))
  bundle$annotations <- A
  bundle$annotations[1, 2] <- 1 - bundle$annotations[1, 2]
  changed <- sblrbench:::.annotation_checkpoint_identities(study06_spec,
    coordinate, controls, data, bundle, "qualification")
  expect_false(identical(first$checkpoint_hash, changed$checkpoint_hash))
})

annotation_trace_fixture <- function(A) {
  annotations <- colnames(A)
  sticks <- paste0("component_", 0:2, "_stick")
  alpha_q <- expand.grid(annotation_name = annotations, stick_name = sticks,
    stringsAsFactors = FALSE)
  q <- rbind(data.frame(parameter_name = "alpha", alpha_q),
    data.frame(parameter_name = "sigmaSqAlpha",
      annotation_name = annotations, stick_name = "",
      stringsAsFactors = FALSE))
  values <- array(0, dim = c(3L, 4L, nrow(q)))
  values[, , seq_len(nrow(alpha_q))] <- rep(seq(-.4, .4,
    length.out = nrow(alpha_q)), each = 12L)
  values[, , (nrow(alpha_q) + 1L):nrow(q)] <- 1
  list(convergence_traces = list(values = values, quantities = q),
    chains = as.list(seq_len(4L)), alpha = matrix(99, 4, 3),
    alpha_final = matrix(-99, 4, 3), annotation_prior = matrix(1, nrow(A), 4))
}

test_that("BED and CSR use true comparable draw-wise annotation priors", {
  ids <- paste0("m", seq_len(20L))
  A <- study06_logic$construct_annotation_design(ids, study06_spec)
  fixture <- annotation_trace_fixture(A)
  for (interface in c("BED", "CSR")) {
    fit <- structure(fixture, interface = interface)
    traces <- extract_annotation_coefficient_traces(fit, expected_chains = 4L)
    expect_true(all(traces$status == "ok"))
    expect_setequal(unique(traces$parameter), c("alpha", "sigmaSqAlpha"))
    prior <- summarise_drawwise_annotation_prior(traces, A,
      study06_spec$controls$simulation$mixture_var)
    expect_identical(prior$status, "ok")
    expect_equal(nrow(prior$marker), 20L)
    expect_equal(nrow(prior$draws), 12L)
    first <- traces[traces$parameter == "alpha" & traces$chain == 1L &
      traces$iteration == 1L, , drop = FALSE]
    coefficient <- matrix(first$value,
      nrow = length(unique(first$annotation)), byrow = FALSE,
      dimnames = list(unique(first$annotation), unique(first$stick)))
    coefficient <- coefficient[colnames(A), paste0("component_", 0:2,
      "_stick"), drop = FALSE]
    reference <- sblr::sbayesrc_marker_pi(A, coefficient,
      study06_spec$controls$simulation$mixture_var)
    observed <- as.matrix(prior$marker[paste0(
      "posterior_mean_prior_component_", 0:3)])
    expect_equal(observed, reference, tolerance = 1e-12,
      ignore_attr = TRUE)
    gate_prior <- summarise_drawwise_annotation_prior(traces, A,
      study06_spec$controls$simulation$mixture_var,
      retain_marker_summary = FALSE)
    expect_null(gate_prior$marker)
    expect_equal(gate_prior$draws, prior$draws)
  }
  missing <- fixture
  missing$convergence_traces <- NULL
  unavailable <- extract_annotation_coefficient_traces(missing)
  expect_identical(unavailable$status, "unavailable")
  expect_match(unavailable$reason, "not substitutes")
})

test_that("annotation recovery metrics reproduce deterministic fixtures", {
  metadata <- list(study = "06_annotation_models",
    scenario = "informative_annotations", replicate = 1L,
    method = "st_bed_bayesrc")
  truth <- data.frame(marker_id = paste0("m", 1:4),
    true_nonnull = c(TRUE, FALSE, TRUE, FALSE), effect = c(1, 0, -1, 0))
  marker <- data.frame(marker_id = truth$marker_id,
    posterior_inclusion_probability = c(.9, .1, .8, .2),
    posterior_mean_effect = c(.8, .1, -.9, 0))
  metrics <- annotation_marker_recovery(marker, truth, metadata, top_k = 2L)
  expect_equal(metrics$value[metrics$metric == "causal_marker_rank"], 1.5)
  expect_equal(metrics$value[metrics$metric == "top_2_recovery"], 1)
  expect_equal(metrics$value[metrics$metric == "pip_auprc"], 1)
  expect_s3_class(plot_annotation_marker_recovery(metrics), "ggplot")
  prior_metrics <- data.frame(scenario = "informative_annotations",
    method = "st_bed_bayesrc", metric = "enriched_prior_contrast",
    value = .2, status = "ok")
  expect_s3_class(plot_annotation_prior_recovery(prior_metrics), "ggplot")
  parameter <- data.frame(scenario = "informative_annotations",
    method = "st_bed_bayesrc", parameter = "heritability", truth = .3,
    posterior_mean = .31, lower_95 = .25, upper_95 = .36, status = "ok")
  expect_s3_class(plot_annotation_parameter_recovery(parameter), "ggplot")
})

test_that("validate-only resolves Study 06 without any fit dispatch", {
  old <- options(sblrbench.annotation_fit_dispatch = function(...)
    stop("sampler called"))
  on.exit(options(old), add = TRUE)
  out <- run_benchmark(study06_spec, tempfile("study06-validation-"),
    profile = "benchmark", validate_only = TRUE)
  expect_equal(nrow(out$qualification_coordinates), 4L)
  expect_equal(nrow(out$coordinate_grid), 40L)
  expect_equal(nrow(out$status), 40L)
  expect_identical(out$status$status, rep("not_run_validate_only", 40L))
})

test_that("final execution requires an identity-matched qualification", {
  output <- tempfile("study06-final-")
  expect_error(run_benchmark(study06_spec, output, mode = "final"),
    "qualification decision")
  decision <- annotation_qualification_artifact_schema()
  expect_true(all(c("semantic_checkpoint_hash", "reusable_history_hash") %in%
    decision$entry_fields))
  entry <- list(scenario = "informative_annotations", replicate = 1L,
    method = "st_bed_bayesrc", available_history = 9000L,
    selected_burnin = 1000L, selected_retained = 2000L, rhat = 1,
    ess_bulk = 500, ess_tail = 500, relative_mcse = .01,
    all_quantities_pass = TRUE,
    quantity_decisions = list(list(quantity = "alpha", rhat = 1,
      ess_bulk = 500, ess_tail = 500, relative_mcse = .01, pass = TRUE)),
    semantic_checkpoint_hash = "checkpoint", reusable_history_hash = "history")
  entries <- lapply(seq_len(nrow(study06_spec$qualification$entries)),
    function(i) utils::modifyList(entry, as.list(
      study06_spec$qualification$entries[i, ])))
  stale <- list(schema = decision$schema, study = study06_spec$study,
    spec_hash = "stale", sblr_sha = study06_spec$packages$sblr$sha,
    qgdata_sha = study06_spec$packages$qgdata$sha, entries = entries,
    overall_decision = "passed", created_at = "metadata")
  expect_error(validate_annotation_qualification_decision(stale,
    study06_spec), "stale or incompatible")
})

test_that("Study 06 workflow and report retain safe boundaries", {
  analysis <- readLines(file.path(study06_root, "studies",
    "06_annotation_models", "analysis.R"), warn = FALSE)
  template <- readLines(file.path(study06_root, "inst", "templates",
    "annotation-analysis.R"), warn = FALSE)
  report <- readLines(file.path(study06_root, "studies",
    "06_annotation_models", "report.qmd"), warn = FALSE)
  expect_true(any(grepl('SBLR_BENCH_MODE", "validate_only', analysis,
    fixed = TRUE)))
  expect_true(any(grepl('SBLR_BENCH_MODE", "validate_only', template,
    fixed = TRUE)))
  expect_true(any(grepl("annotation_prior_plot <-", analysis, fixed = TRUE)))
  expect_true(any(grepl("qualification failed", report,
    ignore.case = TRUE)))
  expect_true(any(grepl("Status: In development", report,
    fixed = TRUE)))
  forbidden <- c("^\\s*run_benchmark\\(", "^\\s*stblr_bed\\(",
    "^\\s*stblr_csr")
  expect_false(any(vapply(forbidden, function(x) any(grepl(x, report)),
    logical(1))))
  expect_true(any(grepl("current-stop", report, fixed = TRUE)))
})
