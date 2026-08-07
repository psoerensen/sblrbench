test_that("Study 06 GCTB block contract is pinned and route-specific", {
  env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "studies",
    "06_annotation_models", "gctb-block-contract-validation.R"), envir = env)
  cfg <- env$study06_gctb_block_constants(normalizePath(
    testthat::test_path("..", ".."), winslash = "/"))
  expect_identical(cfg$sblr_sha,
    "0c89234273389e14112ba0e08ef9d11d3e1819dc")
  expect_identical(cfg$official_sha,
    "b95d3fcbad8ff358290922a58fff893439296138")
  expect_identical(cfg$gctb_sha,
    "cc7fa7d765c83a89c6375946cf77fe50ba1a317e")
  registry <- env$study06_gctb_block_registry(cfg)
  expect_identical(registry$fit_id, c("S0_new", "S1_new"))
  expect_true(all(registry$residual_policy == "gctb_block"))
  expect_true(all(registry$block_ve_mode == "allMixVe"))
  expect_true(all(registry$resam_thresh == 1.1))
  expect_true(all(registry$minimum_ve_ratio == .7))
  expect_true(all(registry$niter - registry$burn == 6000L))
})

test_that("Study 06 residual semantic crosswalk does not relabel history", {
  env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "studies",
    "06_annotation_models", "gctb-block-contract-validation.R"), envir = env)
  x <- env$study06_residual_semantic_crosswalk()
  expect_identical(x$policy, c("global_projected_legacy", "gctb_block",
    "global individual-level residual", "global operator residual"))
  expect_false(x$historical_relabeling_allowed[[1L]])
  expect_match(x$ves[[2L]], "mean block residual")
  expect_match(x$heritability[[2L]], "sum\\(block Vg\\)")
})

test_that("frozen large identity retains explicit new block execution contract", {
  env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "studies",
    "06_annotation_models", "large-feasibility.R"), envir = env)
  spec <- env$study06_large_spec()
  expect_identical(env$study06_large_hash(spec),
    "b001bc36a5531e5e6b342286a253fc1fd34dad4265359d89d2feaa026d4533df")
  contract <- env$study06_large_block_contract()
  expect_identical(contract$residual_policy, "gctb_block")
  expect_identical(contract$block_ve_mode, "allMixVe")
  expect_equal(contract$resam_thresh, 1.1)
  expect_equal(contract$minimum_ve_ratio, .7)
  registry <- env$study06_large_registry(spec)
  learned <- env$study06_large_controls(registry[registry$fit_id == "B1", ],
    matrix(0, 4, 3), spec, smoke = FALSE)
  expect_identical(learned$nit, 9000L)
  expect_identical(learned$nburn, 3000L)
  expect_identical(learned$chain_seeds,
    c(760121L, 760222L, 760323L, 760424L))
})

test_that("installed block API exposes GCTB policy and compact outputs", {
  f <- formals(sblr::stblr_block_eigen)
  expect_true(all(c("residual_policy", "block_ve_mode", "resam_thresh",
    "minimum_ve_ratio", "block_ve_keep_history") %in% names(f)))
  expect_equal(eval(f$resam_thresh), 1.1)
  expect_equal(eval(f$minimum_ve_ratio), .7)
  expect_identical(eval(f$block_ve_mode)[[1L]], "allMixVe")
})

test_that("Study 06 contract and large decisions preserve historical meaning", {
  root <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
  small <- jsonlite::read_json(file.path(root, "docs", "dev",
    "study06_gctb_block_contract_validation_decision.json"),
    simplifyVector = TRUE)
  large <- jsonlite::read_json(file.path(root, "docs", "dev",
    "study06_large_feasibility_decision.json"), simplifyVector = TRUE)
  expect_identical(small$decision, "SMALL-G1")
  expect_true(small$s0_all_gates_passed)
  expect_true(small$s1_all_gates_passed)
  expect_false(small$formal_qualification_changed)
  expect_identical(large$decision, "LARGE-G2")
  expect_identical(large$historical_phase$decision, "LARGE-F6")
  expect_identical(large$small_contract_gate, "SMALL-G1")
  expect_identical(large$smoke_gate, "SMOKE-G1")
  expect_equal(large$additional_replicates, 0)
  expect_equal(large$changed_seed_retries, 0)
  expect_false(large$final_benchmark_authorized)
})

test_that("large completed-analysis helpers are read-only and deterministic", {
  env <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "studies",
    "06_annotation_models", "large-feasibility-analysis.R"), envir = env)
  x <- matrix(c(1, 2, 3, 4), 2L)
  expect_equal(unname(env$study06_large_draw_summary(x)["mean"]), 2.5)
  pip <- c(.9, .8, .2, .1)
  expect_identical(env$study06_large_bfdr_set(pip, .20), 1:2)
})
