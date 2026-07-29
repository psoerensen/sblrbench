study01_helper <- testthat::test_path(
  "..", "..", "studies", "01_finemapping", "setup_example_data.R"
)
.study01_helpers_available <- file.exists(study01_helper)
if (.study01_helpers_available) source(study01_helper, local = TRUE)

.mock_glist <- function() {
  list(
    ids = c("s1", "s2", "s3"), n = 3L,
    rsids = list(c("m3", "m1", "m2", "m4")),
    rsidsLD = NULL,
    af = list(c(0.20, 0.30, 0.40, 0.10)),
    maf = list(c(0.20, 0.30, 0.40, 0.10))
  )
}

test_that("Study 01 preserves chromosome marker order", {
  skip_if_not(.study01_helpers_available, "repository-only Study 01 helpers are excluded from package builds")
  Glist <- .mock_glist()
  info <- .study01_filter_markers(
    Glist, 1L, c("m2", "m3", "m1"), 0.05, 1000L
  )
  expect_identical(info$marker_ids, c("m3", "m1", "m2"))
  expect_identical(info$marker_count_before, 4L)
  expect_identical(info$marker_count_after, 3L)
})

test_that("Study 01 marker validation rejects invalid sets", {
  skip_if_not(.study01_helpers_available, "repository-only Study 01 helpers are excluded from package builds")
  Glist <- .mock_glist()
  expect_error(
    .study01_filter_markers(Glist, 1L, c("m1", "m1"), 0.05, 2L),
    "unique"
  )
  expect_error(
    .study01_filter_markers(Glist, 1L, c("absent", "m1"), 0.05, 2L),
    "absent"
  )
  expect_error(
    .study01_filter_markers(Glist, 1L, character(), 0.05, 2L),
    "No markers"
  )
})

test_that("working Glist receives rsidsLD without altering input", {
  skip_if_not(.study01_helpers_available, "repository-only Study 01 helpers are excluded from package builds")
  Glist <- .mock_glist()
  original <- Glist
  working <- .study01_set_rsids_ld(Glist, 1L, c("m3", "m1"))
  expect_identical(working$rsidsLD[[1L]], c("m3", "m1"))
  expect_identical(Glist, original)
})

test_that("qgg extraction preserves requested identities", {
  skip_if_not(.study01_helpers_available, "repository-only Study 01 helpers are excluded from package builds")
  skip_if_not_installed("qgg")
  bed <- system.file("extdata", "sample_chr1.bed", package = "qgg")
  bim <- system.file("extdata", "sample_chr1.bim", package = "qgg")
  fam <- system.file("extdata", "sample_chr1.fam", package = "qgg")
  skip_if(!nzchar(bed) || !file.exists(bed),
          "Installed qgg has no bundled PLINK example")
  Glist <- qgg::gprep(study = "sblrbench-test", bedfiles = bed,
                      bimfiles = bim, famfiles = fam)
  ids <- rev(head(Glist$ids, 4L))
  rsids <- rev(head(Glist$rsids[[1L]], 5L))
  Z <- .study01_extract_genotypes(Glist, 1L, ids, rsids)
  expect_identical(rownames(Z), ids)
  expect_identical(colnames(Z), rsids)
  expect_true(all(is.finite(Z)))
})

test_that("Study 01 simulation keeps qgg scale and oracle detects changes", {
  skip_if_not(.study01_helpers_available, "repository-only Study 01 helpers are excluded from package builds")
  skip_if_not_installed("sblr")
  set.seed(42)
  Z <- scale(matrix(sample(0:2, 120, replace = TRUE), 20, 6))
  dimnames(Z) <- list(paste0("sample", seq_len(nrow(Z))),
                      paste0("marker", seq_len(ncol(Z))))
  raw <- sblr::mtsim(W = Z, standardize_W = FALSE, nt = 1L,
                     n_shared = 2L, n_specific = 0L, h2 = 0.2,
                     seed = 1001L)
  simulation <- as_sblrbench_simulation(
    raw, study = "finemapping", architecture = "scale_validation",
    replicate = 0L, seed = 1001L
  )
  expect_identical(simulation$data$sample_ids, rownames(Z))
  expect_identical(simulation$data$marker_ids, colnames(Z))
  expect_true(check_oracle_genetic_values(simulation, tolerance = 1e-10)$ok)
  changed <- simulation
  changed$data$genotypes <- changed$data$genotypes * 1.01
  expect_error(
    check_oracle_genetic_values(changed, tolerance = 1e-10,
                                stop_on_failure = TRUE),
    "failed"
  )
})
