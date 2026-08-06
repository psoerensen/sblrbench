test_that("validated package provenance is strict", {
  provenance <- benchmark_package_provenance("sblr")
  expect_identical(provenance$version, "0.2.0")
  expect_match(provenance$installed_tree_sha256, "^[0-9a-f]{64}$")
  if (isTRUE(provenance$remote_sha_available)) {
    expect_invisible(benchmark_assert_package_sha("sblr", provenance$sha))
    expect_error(benchmark_assert_package_sha(
      "sblr", paste0("x", provenance$sha)), "SHA mismatch")
  } else {
    expect_true(is.na(provenance$sha))
  }
})

test_that("checkpoint save and strict reuse preserve identity", {
  path <- tempfile(fileext = ".rds")
  value <- list(input_hash = benchmark_hash_object(list(seed = 17L)),
    payload = 1:3)
  expect_identical(benchmark_atomic_save_rds(value, path), path)
  loaded <- benchmark_load_checkpoint(path, value$input_hash)
  expect_true(loaded$reused)
  expect_identical(loaded$value, value)
  expect_error(benchmark_load_checkpoint(path, "different"),
    "refusing reuse")
})

test_that("capsule checksums detect post-inventory changes", {
  path <- tempfile("capsule-")
  dir.create(path)
  writeLines("stable", file.path(path, "table.csv"))
  writeLines("{}", file.path(path, "benchmark_manifest.json"))
  required <- c("table.csv", "benchmark_manifest.json", "checksums.csv")
  utils::write.csv(benchmark_capsule_checksums(path),
    file.path(path, "checksums.csv"), row.names = FALSE)
  expect_invisible(benchmark_validate_capsule_checksums(path, required))
  writeLines("changed", file.path(path, "table.csv"))
  expect_error(benchmark_validate_capsule_checksums(path, required),
    "checksum validation failed")
})

test_that("shared scalar diagnostics match the direct calculation", {
  draws <- expand.grid(iteration = seq_len(40L), chain = seq_len(4L))
  draws$value <- sin(draws$iteration / 5) + draws$chain / 100
  thresholds <- list(rhat = 2, ess_bulk = 1, ess_tail = 1,
    relative_mcse = 2, chain_count = 4L)
  observed <- benchmark_scalar_diagnostics(draws, thresholds)
  wide <- reshape(draws, idvar = "iteration", timevar = "chain",
    direction = "wide")
  matrix <- as.matrix(wide[-1L])
  expect_equal(observed$rhat, posterior::rhat(matrix))
  expect_equal(observed$ess_bulk, posterior::ess_bulk(matrix))
  expect_equal(observed$ess_tail, posterior::ess_tail(matrix))
  expect_equal(observed$mcse_mean, posterior::mcse_mean(matrix))
  expect_identical(observed$chain_count, 4L)
  expect_identical(observed$draws_per_chain, 40L)
})
