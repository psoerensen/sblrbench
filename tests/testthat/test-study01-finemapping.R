study01_spec <- function() read_benchmark_spec(test_path("..","..","studies",
  "01_finemapping","spec.R"))

test_that("Study 01 spec, coordinates, and historical seeds are preserved", {
  spec <- study01_spec()
  expect_invisible(validate_benchmark_spec(spec))
  coordinates <- benchmark_seeds(spec,"benchmark")
  expect_equal(nrow(coordinates),40L)
  expect_equal(unique(coordinates$scenario),"separated")
  expect_equal(unique(coordinates$replicate),1:10)
  expect_equal(coordinates$causal_seed[1:4],rep(2002L,4))
  expect_equal(coordinates$simulation_seed[1:4],rep(3002L,4))
  expect_equal(coordinates$fit_seed[1:4],10000L+100L+1:4)
  expect_equal(nrow(benchmark_coordinates(spec,"workshop")),4L)
})

test_that("Study 01 controls and credible-set policy remain exact", {
  spec <- study01_spec()
  row <- benchmark_seeds(spec,"benchmark")[1,,drop=FALSE]
  controls <- benchmark_method_controls(spec,row$method,"benchmark",
    row$fit_seed,row$chain_seeds[[1L]])
  expect_equal(unname(unlist(controls[c("nit","nburn","nthin","nchains",
    "ncores")])),c(500,250,1,1,1))
  expect_equal(spec$locus_design$credible_set_target,.95)
  expect_equal(spec$locus_design$min_r2,.5)
  expect_identical(spec$locus_design$algorithm,"sblr::make_credible_sets")
})

test_that("separated causal selection preserves marker order and distance", {
  environment <- new.env(parent=baseenv())
  sys.source(test_path("..","..","studies","01_finemapping",
    "locus-design.R"),environment)
  glist <- list(rsids=list(paste0("m",1:6)),pos=list(seq(0,by=1e6,length.out=6)),
    maf=list(rep(.2,6)))
  selection <- environment$select_separated_causal_markers(glist,1L,
    glist$rsids[[1]],3L,1e6,2002L,.05,.5)
  expect_equal(selection$marker_ids,
    glist$rsids[[1]][sort(match(selection$marker_ids,glist$rsids[[1]]))])
  expect_gte(selection$pairwise_min_distance_bp,1e6)
})

test_that("fine-mapping metric and plotting fixtures are deterministic", {
  marker <- data.frame(study="01_finemapping",scenario="separated",replicate=1L,
    method="st_bed_bayesc",marker=paste0("m",1:4),trait="trait1",
    true_effect=c(1,0,0,0),posterior_mean_effect=c(.8,.1,0,0),
    posterior_inclusion_probability=c(.8,.1,.05,.05),causal=c(TRUE,FALSE,FALSE,FALSE),
    causal_rank=1:4)
  metrics <- evaluate_metrics_from_tables(marker)
  expect_equal(metrics$value[metrics$metric=="pip_brier"],.01375)
  expect_s3_class(plot_causal_marker_pip(marker),"ggplot")
  expect_s3_class(plot_causal_marker_rank(marker),"ggplot")
})

test_that("Study 01 validation-only execution cannot dispatch a fit", {
  spec <- study01_spec()
  old <- options(sblrbench.finemapping_fit_dispatch=function(...)
    stop("fit dispatch must not run"))
  on.exit(options(old),add=TRUE)
  out <- tempfile("study01-validation-")
  result <- run_benchmark(spec,out,"workshop",validate_only=TRUE)
  expect_equal(nrow(result$status),4L)
  expect_true(all(result$status$status=="not_run_validate_only"))
})

test_that("Study 01 report remains capsule-only", {
  text <- paste(readLines(test_path("..","..","studies","01_finemapping",
    "report.qmd"),warn=FALSE),collapse="\n")
  expect_match(text,"results/reference/01_finemapping/current",fixed=TRUE)
  expect_false(grepl("run_benchmark\\(|readRDS\\(|targets::tar_",text))
})
