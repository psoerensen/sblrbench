testthat::test_that("Study 06 stick transform preserves the probability contract", {
  source(file.path("..", "..", "studies", "06_annotation_models",
    "estimability-and-contrasts.R"), local = TRUE)
  eta <- matrix(c(-2, 0, 2, 1, -1, .5), nrow = 2, byrow = TRUE)
  pi <- study06_pi(eta)
  testthat::expect_true(all(is.finite(pi)))
  testthat::expect_true(all(pi >= 0 & pi <= 1))
  testthat::expect_equal(rowSums(pi), c(1, 1), tolerance = 1e-14)
  testthat::expect_equal(pi[, 1], 1 - stats::pnorm(eta[, 1]))
})

testthat::test_that("Study 06 ranking metrics are exact on simple inputs", {
  source(file.path("..", "..", "studies", "06_annotation_models",
    "estimability-and-contrasts.R"), local = TRUE)
  testthat::expect_equal(study06_auc(4:1, c(TRUE, TRUE, FALSE, FALSE)), 1)
  testthat::expect_equal(study06_auprc(4:1, c(TRUE, TRUE, FALSE, FALSE)), 1)
  metrics <- study06_pair_metrics(1:100, 1:100, top = c(25L, 50L, 100L))
  testthat::expect_equal(unname(metrics[c("pearson", "spearman", "rmse", "mae")]), c(1, 1, 0, 0))
  testthat::expect_equal(unname(metrics[c("top25_overlap", "top50_overlap", "top100_overlap")]), c(1, 1, 1))
})

testthat::test_that("raw-alpha summaries consume complete matrix columns", {
  source(file.path("..", "..", "studies", "06_annotation_models",
    "estimability-and-contrasts.R"), local = TRUE)
  descriptor <- data.frame(annotation_index = 1:2, stick_index = c(1L, 1L),
    annotation_name = c("a", "b"), stick_name = c("s", "s"))
  alpha <- list(descriptor = descriptor,
    traces = list(cbind(1:20, 101:120), cbind(21:40, 121:140),
      cbind(41:60, 141:160), cbind(61:80, 161:180)))
  truth <- matrix(0, 2, 1)
  result <- study06_raw_alpha("test", alpha, truth)$summary
  testthat::expect_equal(result$posterior_mean, c(40.5, 140.5))
})
