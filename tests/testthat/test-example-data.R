test_that("example-data URLs are pinned and names are stable", {
  urls <- sapply(sblrbench:::.sblrbench_example_data$filename,
                 sblrbench:::.sblrbench_example_data_url)
  expect_identical(names(urls), sblrbench:::.sblrbench_example_data$filename)
  expect_true(all(grepl(sblrbench:::.sblrbench_qgdata_commit, urls, fixed = TRUE)))
  expect_false(any(grepl("/main/", urls, fixed = TRUE)))
  expect_identical(basename(urls), sblrbench:::.sblrbench_example_data$filename)
})

test_that("checksum validation detects changed fixtures", {
  path <- tempfile(); writeBin(charToRaw("known fixture"), path)
  size <- unname(file.info(path)$size); md5 <- unname(tools::md5sum(path))
  expect_true(sblrbench:::.sblrbench_validate_example_file(path, size, md5))
  writeBin(charToRaw("changed fixture"), path)
  expect_false(sblrbench:::.sblrbench_validate_example_file(path, size, md5))
})

test_that("invalid existing files are protected from overwrite", {
  destination <- tempfile(); dir.create(destination)
  writeLines("corrupt", file.path(destination, "human.bed"))
  expect_error(download_sblrbench_example_data(destination, overwrite = FALSE, quiet = TRUE),
               "failed size/checksum")
})
