.study03_spec <- function() read_benchmark_spec(testthat::test_path("..","..",
  "studies","03_parameter_estimation","spec.R"))

test_that("Study 03 spec, profiles, coordinates, and seeds are exact", {
  spec <- .study03_spec(); expect_invisible(validate_benchmark_spec(spec))
  expect_equal(resolve_benchmark_profile(spec,"workshop")$replicate_count,1L)
  expect_equal(resolve_benchmark_profile(spec,"benchmark")$replicate_count,5L)
  x <- benchmark_seeds(spec,"benchmark")
  expect_equal(nrow(x),40L)
  expect_identical(unique(x$scenario),c("sparse_homogeneous","sparse_mixture"))
  expect_identical(unique(x$method),names(spec$methods))
  expect_identical(x$architecture_seed[1],6001L)
  expect_identical(x$simulation_seed[1],6002L)
  expect_identical(x$fit_seed[1],30101L)
  expect_identical(x$chain_seeds[[1]],c(130101L,230101L,330101L,430101L))
  bad <- spec; bad$estimands <- NULL
  expect_error(validate_benchmark_spec(bad),"estimands")
})

test_that("Study 03 simulation and realized truth are deterministic", {
  spec <- .study03_spec(); spec$controls$simulation$n_causal <- 10L
  set.seed(9182L); z <- matrix(rnorm(3000),100,30,
    dimnames=list(paste0("s",1:100),paste0("m",1:30)))
  expected <- c(sparse_homogeneous=
    "b88f71bf594303ad79d7b7f8924d5116ba3ceeb304bd836a300c1b3a1f873401",
    sparse_mixture=
    "0cfef6247c447a95ba0ee94c1728d1197965db74d278296a393c861e07150671")
  for(scenario in names(spec$scenarios)) {
    sim <- simulate_prediction_architecture(list(scenario=scenario,replicate=1L,
      simulation_seed=5001L+match(scenario,names(spec$scenarios))),z,spec)
    value <- digest::digest(list(effects=sim$truth$effects,
      phenotypes=sim$truth$phenotypes,causal=sim$truth$causal,
      extras=sim$extras),algo="sha256")
    expect_identical(value,unname(expected[[scenario]]))
    truth <- sblrbench:::parameter_estimand_truth(sim,spec)
    expect_equal(truth$truth[truth$estimand_id=="causal_proportion"],1/3)
    expect_equal(truth$truth[truth$estimand_id=="heritability"],.3,
      tolerance=1e-12)
  }
})

test_that("parameter draw transformations and metrics preserve definitions", {
  fit <- list(input=list(nburn=2L,nit=4L,nthin=1L),
    pi_trace=matrix(c(0,0,.1,.9)),vbs=matrix(c(0,0,9,1)),
    vgs=matrix(c(0,0,1,3)),ves=matrix(c(0,0,3,1)))
  spec <- .study03_spec()
  d <- sblrbench:::extract_parameter_draws(fit,"st_bed_bayesc",
    spec$estimands,10,1L)
  expect_equal(d$value[d$estimand_id=="heritability"],c(.25,.75))
  expect_equal(d$value[d$estimand_id=="total_marker_effect_variance"],c(9,9))
  summary <- sblrbench:::summarise_parameter_draws(d)
  truth <- data.frame(estimand_id=summary$estimand_id,
    architecture="sparse_homogeneous",replicate=1L,truth=summary$posterior_mean,
    truth_type="realized_quantity",status="ok",reason="")
  method <- sblrbench:::resolve_benchmark_methods(spec)[[1L]]
  recovery <- sblrbench:::parameter_recovery_metrics(summary,truth,method,
    spec$estimands)
  expect_true(all(recovery$bias==0)); expect_true(all(recovery$covered_95))
})

test_that("Study 03 controls and validate-only dispatch are preserved", {
  spec <- with_current_sblr_description_pin(.study03_spec())
  coords <- benchmark_seeds(spec,"benchmark")
  for(id in names(spec$methods)) {
    row <- coords[coords$scenario=="sparse_homogeneous" &
      coords$replicate==1L & coords$method==id,]
    ctl <- sblrbench:::benchmark_method_controls(spec,id,"benchmark",
      row$fit_seed,row$chain_seeds[[1]])
    expect_identical(ctl$nburn,250L); expect_identical(ctl$nchains,4L)
    expect_identical(ctl$nit,if(grepl("bayesr",id))2000L else 250L)
  }
  withr::local_options(list(sblrbench.parameter_fit_dispatch=function(...)
    stop("FIT MUST NOT RUN")))
  out <- tempfile("study03-validation-")
  result <- run_benchmark(spec,out,"benchmark",validate_only=TRUE)
  expect_equal(nrow(result$status),40L)
  expect_true(all(result$status$status=="not_run_validate_only"))
})

test_that("Study 03 frozen capsule and report contract validate", {
  root <- testthat::test_path("..","..")
  spec <- .study03_spec(); capsule <- file.path(root,"results","reference",
    "03_parameter_estimation","current")
  expect_invisible(sblrbench:::validate_parameter_capsule(capsule,spec))
  report <- readLines(file.path(root,"studies","03_parameter_estimation",
    "report.qmd"),warn=FALSE)
  fences <- cumsum(grepl("^```",report)); executable <- paste(
    report[fences%%2L==1L & !grepl("^```",report)],collapse="\n")
  expect_false(grepl("results/local|readRDS|run_benchmark|targets::|tar_read|stblr_(bed|csr)\\(",executable))
})
