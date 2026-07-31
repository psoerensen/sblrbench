.study06_blocks <- function(marker_ids, target_size) {
  stopifnot(is.character(marker_ids), length(marker_ids) > 0L,
    !anyNA(marker_ids), !anyDuplicated(marker_ids),
    length(target_size) == 1L, is.finite(target_size), target_size > 0)
  m <- length(marker_ids)
  starts <- seq.int(1L, m, by = as.integer(target_size))
  ends <- c(starts[-1L] - 1L, m)
  data.frame(
    block_id = sprintf("block_%04d", seq_along(starts)),
    start = starts, end = ends, size = ends - starts + 1L,
    first_marker = marker_ids[starts], last_marker = marker_ids[ends],
    stringsAsFactors = FALSE)
}

.study06_validate_blocks <- function(blocks, marker_ids) {
  required <- c("block_id", "start", "end", "size",
    "first_marker", "last_marker")
  if (!is.data.frame(blocks) || !all(required %in% names(blocks)) ||
      !nrow(blocks) || anyNA(blocks[required]) ||
      blocks$start[1L] != 1L ||
      !identical(blocks$start[-1L], blocks$end[-nrow(blocks)] + 1L) ||
      blocks$end[nrow(blocks)] != length(marker_ids) ||
      any(blocks$start > blocks$end) ||
      any(blocks$size != blocks$end - blocks$start + 1L) ||
      anyDuplicated(blocks$block_id) ||
      !identical(blocks$first_marker, marker_ids[blocks$start]) ||
      !identical(blocks$last_marker, marker_ids[blocks$end]))
    stop("Study 06 block coverage or marker order is invalid.", call. = FALSE)
  covered <- unlist(Map(seq.int, blocks$start, blocks$end),
    use.names = FALSE)
  if (!identical(covered, seq_along(marker_ids)))
    stop("Study 06 blocks contain a gap, overlap, or duplicate.", call. = FALSE)
  invisible(TRUE)
}

.study06_block_candidates <- function(marker_ids, config) {
  do.call(rbind, lapply(config$block_candidates, function(size) {
    b <- .study06_blocks(marker_ids, size)
    .study06_validate_blocks(b, marker_ids)
    data.frame(candidate_id = paste0("equal_", size),
      target_size = size, block_count = nrow(b),
      minimum_block_size = min(b$size),
      median_block_size = stats::median(b$size),
      maximum_block_size = max(b$size),
      dense_upper_entries = sum(b$size * (b$size + 1) / 2),
      estimated_dense_storage_bytes = 4 *
        sum(b$size * (b$size + 1) / 2), stringsAsFactors = FALSE)
  }))
}

.study06_cross_block_summary <- function(csr, marker_count, target_size) {
  if (length(csr$row_ptr) != marker_count + 1L ||
      length(csr$col_idx) != length(csr$values))
    stop("Full CSR arrays have inconsistent dimensions.", call. = FALSE)
  counts <- diff(csr$row_ptr)
  row <- rep.int(seq_len(marker_count), counts)
  col <- as.integer(csr$col_idx) + 1L
  if (length(row) != length(col) || any(col <= row) ||
      any(col > marker_count))
    stop("Full CSR upper-triangle indexing is invalid.", call. = FALSE)
  block_row <- (row - 1L) %/% as.integer(target_size)
  block_col <- (col - 1L) %/% as.integer(target_size)
  removed <- block_row != block_col
  total_sq <- sum(csr$values^2)
  total_abs <- sum(abs(csr$values))
  data.frame(
    target_size = as.integer(target_size),
    cross_block_edge_count_removed = sum(removed),
    cross_block_squared_ld_mass_removed = sum(csr$values[removed]^2),
    cross_block_absolute_ld_mass_removed = sum(abs(csr$values[removed])),
    fraction_full_csr_squared_ld_mass_retained =
      if (total_sq > 0) 1 - sum(csr$values[removed]^2) / total_sq else 1,
    fraction_full_csr_absolute_ld_mass_retained =
      if (total_abs > 0) 1 - sum(abs(csr$values[removed])) / total_abs else 1,
    maximum_removed_cross_block_edge =
      if (any(removed)) max(abs(csr$values[removed])) else 0,
    stringsAsFactors = FALSE)
}

.study06_seed <- function(config, architecture, replicate = 1L,
                          configuration = NULL, chain = NULL) {
  a <- match(architecture, config$architectures)
  if (is.na(a) || replicate < 1L || replicate > config$replicate_count)
    stop("Invalid Study 06 seed coordinates.", call. = FALSE)
  out <- config$seeds$simulation_base +
    (a - 1L) * config$seeds$architecture_stride +
    as.integer(replicate) * config$seeds$replicate_stride
  if (!is.null(configuration)) {
    canonical_configurations <- c("bed", "full_csr", "block_csr",
      "block_eigen_unfiltered", "block_eigen_hard",
      "block_eigen_ridge_fixed")
    k <- match(configuration, canonical_configurations)
    if (is.na(k)) stop("Unknown Study 06 configuration.", call. = FALSE)
    out <- out + k * config$seeds$configuration_stride +
      config$seeds$fit_base
  }
  if (!is.null(chain)) {
    if (!chain %in% 1:4) stop("Study 06 chain must be 1--4.", call. = FALSE)
    out <- out + as.integer(chain) * config$seeds$chain_stride
  }
  as.integer(out %% .Machine$integer.max)
}
