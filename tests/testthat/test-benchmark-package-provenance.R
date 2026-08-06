current_experiment_provenance_fixture <- function(...) {
  modifyList(list(
    source_sha = "a165fb0635afcb8a712e8658175dfbb19896b3c3",
    source_clean = TRUE, source_status = character(),
    remote_sha = NA_character_, remote_sha_available = FALSE,
    installed_tree_sha256 = paste(rep("a", 64L), collapse = "")), list(...))
}

test_that("historical package pins validate independently of the loaded package", {
  historical <- "f2e3647920ed7e8b1ea9d47a6571b3753285682a"
  expect_invisible(benchmark_assert_recorded_package_sha(historical,
    historical))
  expect_error(benchmark_assert_recorded_package_sha(
    paste0("0", substring(historical, 2L)), historical),
    "Recorded package SHA mismatch")
})

test_that("current local installs use source identity and installed-tree hash", {
  current <- "a165fb0635afcb8a712e8658175dfbb19896b3c3"
  provenance <- current_experiment_provenance_fixture()
  expect_invisible(benchmark_assert_experiment_package(provenance, current))
  expect_true(is.na(provenance$remote_sha))
  expect_false(provenance$remote_sha_available)
})

test_that("a conflicting installed RemoteSha remains a hard mismatch", {
  current <- "a165fb0635afcb8a712e8658175dfbb19896b3c3"
  provenance <- current_experiment_provenance_fixture(
    remote_sha = "f2e3647920ed7e8b1ea9d47a6571b3753285682a",
    remote_sha_available = TRUE)
  expect_error(benchmark_assert_experiment_package(provenance, current),
    "RemoteSha mismatch")
})

test_that("current experiment rejects source mismatch, dirt, and missing identity", {
  current <- "a165fb0635afcb8a712e8658175dfbb19896b3c3"
  expect_error(benchmark_assert_experiment_package(
    current_experiment_provenance_fixture(source_sha = paste0("0",
      substring(current, 2L))), current), "source SHA mismatch")
  expect_error(benchmark_assert_experiment_package(
    current_experiment_provenance_fixture(source_clean = FALSE,
      source_status = " M src/file.cpp"), current), "repository is dirty")
  expect_error(benchmark_assert_experiment_package(
    current_experiment_provenance_fixture(source_sha = NULL), current),
    "source identity is missing")
  expect_error(benchmark_assert_experiment_package(
    current_experiment_provenance_fixture(installed_tree_sha256 = NA_character_),
    current), "tree identity is missing")
})
