reporting_path <- testthat::test_path("..", "..", "R", "benchmark-reporting.R")
if (!file.exists(reporting_path)) {
  testthat::test_that("reporting helpers are excluded from package builds", {
    testthat::skip("website-only helpers are excluded from package builds")
  })
} else local({
  source(reporting_path, local = TRUE)
  testthat::test_that("method and architecture mappings are complete and stable", {
    ids <- c("st_bed_bayesc","st_bed_bayesr","st_csr_sbayesc","st_csr_sbayesr")
    testthat::expect_identical(sblrbench_method_levels(), ids)
    for (x in list(sblrbench_method_labels(),sblrbench_method_colours(),
                   sblrbench_method_shapes(),sblrbench_method_linetypes()))
      testthat::expect_setequal(names(x), ids)
    testthat::expect_identical(sblrbench_architecture_levels(),
      c("sparse_homogeneous","sparse_mixture"))
    testthat::expect_identical(as.character(sblrbench_method_factor(rev(ids))),
      unname(sblrbench_method_labels()[rev(ids)]))
  })
  testthat::test_that("formatters preserve missingness and readable units", {
    testthat::expect_equal(format_sblrbench_proportion(.12345,3),"0.123")
    testthat::expect_equal(format_sblrbench_runtime(c(30,120)),c("30.0 s","2.0 min"))
    testthat::expect_match(format_sblrbench_interval(1,.5,1.5,2),"1.00")
    testthat::expect_equal(format_sblrbench_number(NA),"—")
  })
  testthat::test_that("replicate preparation never fabricates one-replicate SD", {
    one <- data.frame(method="a",replicate=1,value=2)
    p1 <- prepare_sblrbench_replicates(one,"method")
    testthat::expect_true(is.na(p1$summary$sd))
    many <- data.frame(method="a",replicate=1:3,value=c(1,2,3))
    p2 <- prepare_sblrbench_replicates(many,"method")
    testthat::expect_equal(p2$summary$n_replicates,3)
    testthat::expect_equal(p2$summary$mean,2)
    testthat::expect_equal(p2$summary$minimum,1)
  })
  testthat::test_that("code display preserves text and rejects missing files", {
    f <- tempfile(fileext=".R"); writeBin(charToRaw("x <- 1\r\ny <- 2"),f)
    out <- testthat::capture_output(display_capsule_script(f))
    testthat::expect_match(out,"x <- 1")
    testthat::expect_match(out,"y <- 2")
    collapsed <- testthat::capture_output(display_capsule_script_collapsed(f))
    testthat::expect_match(collapsed,"<details>",fixed=TRUE)
    testthat::expect_error(display_capsule_script(paste0(f,"-missing")),"does not exist")
  })
})
