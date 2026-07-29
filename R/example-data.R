.sblrbench_example_data <- data.frame(
  filename = c("human.bed", "human.bim", "human.fam", "human.pheno", "human.covar"),
  size_bytes = c(62500003, 1882359, 117786, 92786, 641513),
  md5 = c("e89bea9a6cedd9eeef3fd0a5c807db81", "0105119b04c67b7ac7f66cc5e6680963",
          "3c5db3d9eb7f3fc893c75f6f2b89836d", "6a9e7cb1162e43999c170a363863176d",
          "d06002aa2b1b79bdc4c0e92c21f27ae5"),
  stringsAsFactors = FALSE
)

.sblrbench_qgdata_repository <- "psoerensen/qgdata"
.sblrbench_qgdata_commit <- "6cca5819e711d326cfb2614d7e9d9f34942612cd"
.sblrbench_qgdata_subdirectory <- "simulated_human_data"

.sblrbench_example_data_url <- function(filename) {
  paste0("https://raw.githubusercontent.com/", .sblrbench_qgdata_repository,
         "/", .sblrbench_qgdata_commit, "/",
         .sblrbench_qgdata_subdirectory, "/", filename)
}

.sblrbench_validate_example_file <- function(path, expected_size, expected_md5) {
  if (!file.exists(path)) return(FALSE)
  size_ok <- identical(unname(file.info(path)$size), as.numeric(expected_size))
  md5_ok <- identical(unname(tools::md5sum(path)), expected_md5)
  isTRUE(size_ok && md5_ok)
}

#' Download the public simulated genotype example
#'
#' Downloads the five Study 01 PLINK, phenotype, and covariate files from a
#' pinned revision of `psoerensen/qgdata`. Existing files that match the
#' recorded sizes and MD5 checksums are reused.
#'
#' @param destination Directory for downloaded files.
#' @param overwrite Whether a corrupt or incomplete existing file may be
#'   replaced. Valid existing files are always reused.
#' @param quiet Passed to [utils::download.file()].
#' @return A named character vector of validated local paths. Attributes record
#'   the source repository, commit, and subdirectory.
#' @export
download_sblrbench_example_data <- function(
    destination = file.path("results", "local", "example_data"),
    overwrite = FALSE,
    quiet = FALSE) {
  if (!is.character(destination) || length(destination) != 1L || is.na(destination) || !nzchar(destination))
    stop("destination must be one non-empty path.", call. = FALSE)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite))
    stop("overwrite must be TRUE or FALSE.", call. = FALSE)
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet))
    stop("quiet must be TRUE or FALSE.", call. = FALSE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(destination, .sblrbench_example_data$filename)
  for (i in seq_len(nrow(.sblrbench_example_data))) {
    valid <- .sblrbench_validate_example_file(paths[[i]],
      .sblrbench_example_data$size_bytes[[i]], .sblrbench_example_data$md5[[i]])
    if (valid) next
    if (file.exists(paths[[i]]) && !overwrite)
      stop("Existing example-data file failed size/checksum validation: ", paths[[i]],
           ". Set overwrite = TRUE to replace it.", call. = FALSE)
    temporary <- tempfile(pattern = paste0(.sblrbench_example_data$filename[[i]], "-"),
                          tmpdir = destination)
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    utils::download.file(.sblrbench_example_data_url(.sblrbench_example_data$filename[[i]]),
                         temporary, mode = "wb", quiet = quiet)
    if (!.sblrbench_validate_example_file(temporary,
          .sblrbench_example_data$size_bytes[[i]], .sblrbench_example_data$md5[[i]]))
      stop("Downloaded example-data file failed size/checksum validation: ",
           .sblrbench_example_data$filename[[i]], call. = FALSE)
    if (file.exists(paths[[i]])) unlink(paths[[i]])
    if (!file.rename(temporary, paths[[i]]))
      stop("Could not move validated download into place: ", paths[[i]], call. = FALSE)
  }
  names(paths) <- .sblrbench_example_data$filename
  attr(paths, "repository") <- .sblrbench_qgdata_repository
  attr(paths, "commit") <- .sblrbench_qgdata_commit
  attr(paths, "subdirectory") <- .sblrbench_qgdata_subdirectory
  paths
}
