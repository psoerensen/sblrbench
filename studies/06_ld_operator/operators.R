.study06_unpack_triangle <- function(packed, size) {
  expected <- size * (size + 1L) / 2L
  if (length(packed) != expected || any(!is.finite(packed)))
    stop("Invalid packed block-eigen upper triangle.", call. = FALSE)
  out <- matrix(0, size, size)
  k <- 1L
  for (i in seq_len(size)) for (j in i:size) {
    out[i, j] <- out[j, i] <- packed[k]
    k <- k + 1L
  }
  out
}

.study06_inspect_operator <- function(Glist, stats, blocks,
                                      filter = "ridge_fixed",
                                      tau = 0, eta = 0,
                                      effects = NULL) {
  .study06_validate_blocks(blocks, stats$marker_names)
  m <- length(stats$marker_names)
  if (is.null(effects)) effects <- numeric(m)
  if (length(effects) != m) stop("effects length differs from marker count.")
  cls <- stats$cls
  af <- unlist(stats$af, use.names = FALSE)
  if (!identical(stats$marker_names,
      unlist(Map(function(cc, ii) Glist$rsids[[cc]][ii],
        stats$chr, cls), use.names = FALSE)))
    stop("BED, statistic, and operator marker order differ.", call. = FALSE)
  sblr:::stblr_block_eigen_contract_internal(
    bed_files = stats$bed_files, n_bed = as.integer(stats$n_bed),
    cls = cls, rows = stats$rows, af = af,
    block_start = as.integer(blocks$start - 1L),
    wy = do.call(rbind, stats$wy), effects = effects,
    eigen_filter = filter, eigen_tau = tau, eigen_eta = eta,
    validation_mutation = "")
}

.study06_dense_blocks <- function(inspect) {
  Map(.study06_unpack_triangle, inspect$packed_upper_triangle,
    inspect$block_size)
}

.study06_pack_triangle <- function(A) {
  unlist(lapply(seq_len(nrow(A)), function(i) A[i, i:ncol(A)]),
    use.names = FALSE)
}

.study06_apply_blocks <- function(blocks, x) {
  if (sum(vapply(blocks, nrow, integer(1))) != length(x))
    stop("Operator and vector dimensions differ.", call. = FALSE)
  out <- numeric(length(x))
  offset <- 0L
  for (A in blocks) {
    idx <- offset + seq_len(nrow(A))
    out[idx] <- as.numeric(A %*% x[idx])
    offset <- offset + nrow(A)
  }
  out
}

.study06_apply_csr_crossproduct <- function(csr, diagonal, x) {
  m <- length(diagonal)
  if (length(x) != m || length(csr$row_ptr) != m + 1L)
    stop("CSR operator action dimensions differ.", call. = FALSE)
  out <- diagonal * x
  counts <- diff(csr$row_ptr)
  row <- rep.int(seq_len(m), counts)
  col <- as.integer(csr$col_idx) + 1L
  value <- as.numeric(csr$values) * sqrt(diagonal[row] * diagonal[col])
  for (k in seq_along(value)) {
    out[row[k]] <- out[row[k]] + value[k] * x[col[k]]
    out[col[k]] <- out[col[k]] + value[k] * x[row[k]]
  }
  out
}

.study06_action_comparison <- function(reference_action, candidate_action,
                                       vector_id, reference_id,
                                       candidate_id, input_vector = NULL) {
  a <- as.numeric(reference_action); b <- as.numeric(candidate_action)
  if (length(a) != length(b) || any(!is.finite(c(a, b))))
    stop("Operator action comparison is invalid.", call. = FALSE)
  data.frame(vector_id = vector_id, reference_operator_id = reference_id,
    operator_id = candidate_id,
    maximum_absolute_error = max(abs(b - a)),
    relative_norm_difference = sqrt(sum((b - a)^2)) /
      max(sqrt(sum(a^2)), .Machine$double.eps),
    cosine_similarity = sum(a * b) /
      max(sqrt(sum(a^2) * sum(b^2)), .Machine$double.eps),
    quadratic_form_difference = if (is.null(input_vector)) NA_real_ else
      abs(sum(as.numeric(input_vector) * b) -
        sum(as.numeric(input_vector) * a)), stringsAsFactors = FALSE)
}

.study06_write_u32 <- function(path, x) {
  x <- as.double(x)
  if (any(!is.finite(x)) || any(x < 0) || any(x > 2^32 - 1) ||
      any(x != floor(x))) stop("Invalid uint32 values.", call. = FALSE)
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  lo <- as.integer(x %% 256)
  b2 <- as.integer(floor(x / 256) %% 256)
  b3 <- as.integer(floor(x / 65536) %% 256)
  b4 <- as.integer(floor(x / 16777216) %% 256)
  writeBin(as.raw(rbind(lo, b2, b3, b4)), con)
  invisible(path)
}

.study06_write_float32 <- function(path, x) {
  if (any(!is.finite(x))) stop("Non-finite CSR values.", call. = FALSE)
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.numeric(x), con, size = 4L, endian = "little")
  invisible(path)
}

.study06_write_runtime_csr <- function(inspect, prefix, marker_ids,
                                       zero_tolerance = 0) {
  dense <- .study06_dense_blocks(inspect)
  m <- length(marker_ids)
  if (sum(vapply(dense, nrow, integer(1))) != m)
    stop("Runtime CSR marker count mismatch.", call. = FALSE)
  row_cols <- vector("list", m)
  row_vals <- vector("list", m)
  offset <- 0L
  for (A in dense) {
    for (i in seq_len(nrow(A))) {
      if (i >= nrow(A)) next
      j <- seq.int(i + 1L, nrow(A))
      r <- A[i, j] / sqrt(A[i, i] * diag(A)[j])
      keep <- is.finite(r) & abs(r) > zero_tolerance
      row <- offset + i
      row_cols[[row]] <- as.double(offset + j[keep] - 1L)
      row_vals[[row]] <- as.numeric(r[keep])
    }
    offset <- offset + nrow(A)
  }
  counts <- lengths(row_cols)
  row_ptr <- c(0, cumsum(counts))
  col_idx <- unlist(row_cols, use.names = FALSE)
  values <- unlist(row_vals, use.names = FALSE)
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  row_file <- paste0(prefix, ".row_ptr.u64.bin")
  col_file <- paste0(prefix, ".col_idx.u32.0based.bin")
  val_file <- paste0(prefix, ".values.f32.bin")
  sblr:::.stblr_write_uint64_file(row_file, row_ptr)
  .study06_write_u32(col_file, col_idx)
  .study06_write_float32(val_file, values)
  writeLines(c("format=sparse_ld_csr",
    "storage=streamed_upper_triangle", "n_bed=0", "n_used=0",
    "n_samples_used=0", paste0("n_variants=", m),
    paste0("nnz=", length(values)), "triangle=upper",
    "diagonal=implicit_1", paste0("row_ptr_file=", row_file),
    paste0("col_idx_file=", col_file), paste0("values_file=", val_file),
    "row_ptr_type=uint64", "col_idx_type=uint32",
    "values_type=float32", "index_base=0", "value=r"),
    paste0(prefix, ".meta.txt"), useBytes = TRUE)
  writeLines(marker_ids, paste0(prefix, ".markers.txt"), useBytes = TRUE)
  list(prefix = prefix, row_ptr = row_ptr, col_idx = col_idx,
    values = values, nnz = length(values))
}

.study06_runtime_glist <- function(Glist, stats, csr) {
  out <- Glist
  out$sparseLD <- list(prefix = csr$prefix, chr = stats$chr,
    bed_files = stats$bed_files, cls = stats$cls, af = stats$af,
    marker_names = stats$marker_names,
    scale = "standardized_genotype", source = "study06_runtime_matched",
    reference_n = stats$n, rows = stats$rows,
    orientation_status = "bed_coding_by_construction")
  out$rsidsLD <- lapply(stats$chr, function(cc)
    stats$marker_names[stats$marker_metadata$chromosome_or_file == cc])
  out
}

.study06_inspect_from_csr <- function(prefix, runtime_diagonal, blocks,
                                      template) {
  csr <- sblr::sparseLD_read_CSR(prefix, one_based = FALSE)
  m <- length(runtime_diagonal)
  if (csr$nrow != m || length(csr$row_ptr) != m + 1L)
    stop("Runtime CSR dimensions differ from the operator.", call. = FALSE)
  counts <- diff(csr$row_ptr)
  row <- rep.int(seq_len(m), counts)
  col <- as.integer(csr$col_idx) + 1L
  dense <- lapply(seq_len(nrow(blocks)), function(k) {
    idx <- blocks$start[k]:blocks$end[k]
    A <- diag(runtime_diagonal[idx], nrow = length(idx))
    keep <- row %in% idx & col %in% idx
    ii <- row[keep] - blocks$start[k] + 1L
    jj <- col[keep] - blocks$start[k] + 1L
    v <- csr$values[keep] *
      sqrt(runtime_diagonal[row[keep]] * runtime_diagonal[col[keep]])
    A[cbind(ii, jj)] <- v
    A[cbind(jj, ii)] <- v
    A
  })
  out <- template
  out$packed_upper_triangle <- lapply(dense, .study06_pack_triangle)
  out$diagonal <- runtime_diagonal
  out$filter <- list(mode = "runtime_matched_block_csr",
    tau = NA_real_, eta = NA_real_, mu_floor = NA_real_)
  out
}
