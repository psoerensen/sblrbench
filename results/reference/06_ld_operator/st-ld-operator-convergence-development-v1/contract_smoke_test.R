run_study06_contract_smoke_test <- function() {
  write_bed <- function(path, dosage) {
    code <- function(x) ifelse(is.na(x), 1L,
      c(`0` = 3L, `1` = 2L, `2` = 0L)[as.character(x)])
    packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
      values <- as.integer(code(dosage[marker, ]))
      values <- c(values, rep(0L, (-length(values)) %% 4L))
      vapply(seq(1L, length(values), by = 4L), function(i)
        sum(values[i:(i + 3L)] * c(1L, 4L, 16L, 64L)),
        integer(1))
    }))
    con <- file(path, "wb")
    on.exit(close(con), add = TRUE)
    writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), con)
  }
  root <- tempfile("study06-contract-")
  dir.create(root)
  bed <- file.path(root, "tiny.bed")
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 0, 1),
    c(0, 1, 2, 0, 1, 2, 0, 2),
    c(2, 1, 0, 2, 1, 0, 2, 1),
    c(0, 0, 1, 1, 2, 2, 1, 0))
  write_bed(bed, dosage)
  marker_ids <- paste0("rs", 1:4)
  af <- rowMeans(dosage) / 2
  Glist <- list(n = ncol(dosage), bedfiles = bed,
    rsids = list(marker_ids), rsidsLD = list(marker_ids),
    af = list(af))
  stats <- list(
    wy = list(T1 = stats::setNames(c(1, -.5, .25, .75),
      marker_ids)),
    ww = list(T1 = stats::setNames(rep(ncol(dosage), 4),
      marker_ids)),
    yy = c(T1 = 8), n = ncol(dosage), n_bed = ncol(dosage),
    m = 4L, chr = 1L, bed_files = bed, cls = list(1:4),
    af = list(af), rows = 1:ncol(dosage),
    marker_names = marker_ids, trait_names = "T1",
    marker_metadata = data.frame(marker_id = marker_ids,
      chromosome_or_file = 1L, bed_column = 1:4,
      allele_frequency = af))
  blocks <- .study06_blocks(marker_ids, 2L)
  .study06_validate_blocks(blocks, marker_ids)
  exact <- .study06_inspect_operator(Glist, stats, blocks,
    filter = "ridge_fixed", eta = 0)
  lw <- .study06_inspect_operator(Glist, stats, blocks,
    filter = "ridge_lw")
  hard <- .study06_inspect_operator(Glist, stats, blocks,
    filter = "hard_truncate", tau = 0)
  stopifnot(
    identical(exact$filter$mode, "ridge_fixed"),
    identical(exact$filter$eta, 0),
    identical(hard$filter$mu_floor, 0.01),
    identical(lw$filter$mode, "ridge_lw"),
    all(exact$diagonal > 0), all(lw$diagonal > 0),
    all(hard$diagonal > 0))
  prefix <- file.path(root, "runtime")
  csr <- .study06_write_runtime_csr(exact, prefix, marker_ids)
  from_csr <- .study06_inspect_from_csr(prefix, exact$diagonal,
    blocks, exact)
  config_path <- file.path("studies", "06_ld_operator", "config.R")
  if (!file.exists(config_path))
    config_path <- file.path("..", "..", "studies",
      "06_ld_operator", "config.R")
  config <- source(config_path, local = TRUE)$value
  config$operator_tolerance$absolute <- 1e-3
  config$operator_tolerance$relative <- 1e-5
  config$operator_tolerance$product_absolute <- 1e-3
  config$operator_tolerance$quadratic_absolute <- 1e-2
  gate <- .study06_equivalence_gate(from_csr, exact, blocks, config)
  synthetic <- .study06_synthetic_filter_validation()
  stopifnot(isTRUE(gate$summary$pass),
    csr$nnz == 2L, nrow(gate$block_metrics) == 2L,
    all(synthetic$pass))
  unlink(root, recursive = TRUE)
  TRUE
}
