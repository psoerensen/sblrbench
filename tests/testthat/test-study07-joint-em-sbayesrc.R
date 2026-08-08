study07_root <- if (file.exists(file.path("studies",
    "07_joint_em_sbayesrc", "spec.R"))) "." else file.path("..", "..")
study07_dir <- file.path(study07_root, "studies", "07_joint_em_sbayesrc")

test_that("Study 07 freezes one Study 06 replicate and five method routes", {
  expect_true(file.exists(file.path(study07_dir, "spec.R")))
  spec <- source(file.path(study07_dir, "spec.R"), local = TRUE)$value
  expect_identical(spec$study, "07_joint_em_sbayesrc")
  expect_identical(spec$scope$simulated_data_replicates, 1L)
  expect_true(spec$scope$inference_replicates_are_not_data_replicates)
  expect_identical(spec$source$study, "06_annotation_models")
  expect_identical(spec$source$scenario, "informative_annotations")
  expect_identical(spec$source$replicate, 1L)
  expect_identical(spec$methods, c("SBayesR", "SBayesRC", "SBayesRC-EM",
    "SBayesRC-S", "SBayesRC-S-EM"))
  expect_identical(spec$result$decision, "STUDY07-R5")
  expect_true(spec$result$one_replicate_result_is_clear)
  expect_false(spec$result$additional_simulation_replicates_required)
  expect_false(spec$result$capsule_promoted)
  expect_identical(spec$source$specification_hash,
    "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56")
  expect_identical(spec$source$truth_hash,
    "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb")
})

test_that("Study 07 preserves distinct joint and MCEM contracts", {
  spec <- source(file.path(study07_dir, "spec.R"), local = TRUE)$value
  expect_identical(spec$joint$nchains, 4L)
  expect_identical(spec$em$starts,
    c("baseline", "positive", "negative", "mixed"))
  expect_equal(unname(spec$model$sigmaSqAlpha_fixed_em), rep(1, 3))
  expect_equal(spec$model$selection$pi_A_fixed_em, 0.30)
  expect_identical(unname(spec$model$selection$joint_delta_init),
    rep(1L, 3L))
  expect_false(identical(spec$model$selection$joint_delta_init,
    spec$model$selection_truth))
  expect_equal(unname(spec$model$selection$tau2_fixed_em), rep(1, 3))
  expect_false(spec$scope$population_calibration_claim)
  expect_false(spec$scope$scaling_claim)
})

test_that("Study 07 documents the current internal EM routes explicitly", {
  text <- paste(readLines(file.path(study07_dir, "implementation-map.md"),
    warn = FALSE), collapse = "\n")
  expect_match(text, "0.2.0", fixed = TRUE)
  expect_match(text, ".stblr_mcem_sbayesrc_block_eigen", fixed = TRUE)
  expect_match(text, ".stblr_mcem_sbayesrc_s_block_eigen", fixed = TRUE)
  expect_match(text, "annotation_pip_eb", fixed = TRUE)
})

test_that("Study 07 entry points do not invoke a simulation", {
  files <- c("data.R", "methods.R", "extraction.R", "spec.R", "README.md")
  text <- paste(unlist(lapply(file.path(study07_dir, files), readLines,
    warn = FALSE)), collapse = "\n")
  expect_false(grepl("simulate_benchmark_data\\s*\\(", text))
  expect_false(grepl("simulate_annotation", text, fixed = TRUE))
  expect_match(text, "one simulated-data replicate", ignore.case = TRUE)
})

test_that("Study 07 navigation reports the review-pending R5 result", {
  readme <- paste(readLines(file.path(study07_dir, "README.md"), warn = FALSE),
    collapse = "\n")
  report <- paste(readLines(file.path(study07_dir, "report.qmd"), warn = FALSE),
    collapse = "\n")
  expect_match(readme, "STUDY07-R5", fixed = TRUE)
  expect_match(report, "STUDY07-R5 — EM demonstration failed", fixed = TRUE)
  expect_match(report, "one frozen simulation replicate", ignore.case = TRUE)
  expect_match(report, "No capsule was promoted", fixed = TRUE)
})
