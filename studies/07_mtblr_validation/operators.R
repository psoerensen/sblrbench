.study07_source_operator_helpers <- function() {
  needed <- c(".study06_blocks", ".study06_validate_blocks",
    ".study06_inspect_operator", ".study06_dense_blocks",
    ".study06_write_runtime_csr", ".study06_runtime_glist",
    ".study06_inspect_from_csr", ".study06_equivalence_gate",
    ".study06_apply_blocks", ".study06_apply_csr_crossproduct")
  if (!all(vapply(needed, exists, logical(1), inherits = TRUE))) {
    for (f in c("blocks.R", "operators.R", "operator_validation.R"))
      source(file.path("studies", "06_ld_operator", f), local = .GlobalEnv)
  }
  invisible(TRUE)
}

.study07_subset_full_csr <- function(full_csr, marker_count, prefix,
                                     marker_ids) {
  .study07_source_operator_helpers()
  if (marker_count != length(marker_ids) || marker_count > full_csr$nrow)
    stop("Invalid full-CSR subset contract.", call. = FALSE)
  counts <- diff(full_csr$row_ptr)[seq_len(marker_count)]
  rows <- rep.int(seq_len(marker_count), counts)
  starts <- full_csr$row_ptr[seq_len(marker_count)] + 1L
  ends <- full_csr$row_ptr[seq_len(marker_count) + 1L]
  indices <- unlist(Map(function(a, b) if (b < a) integer() else a:b,
    starts, ends), use.names = FALSE)
  cols <- as.integer(full_csr$col_idx[indices]) + 1L
  keep <- cols <= marker_count
  rows <- rows[keep]; cols <- cols[keep]
  values <- as.numeric(full_csr$values[indices][keep])
  per_row <- tabulate(rows, nbins = marker_count)
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  row_file <- paste0(prefix, ".row_ptr.u64.bin")
  col_file <- paste0(prefix, ".col_idx.u32.0based.bin")
  val_file <- paste0(prefix, ".values.f32.bin")
  sblr:::.stblr_write_uint64_file(row_file, c(0, cumsum(per_row)))
  .study06_write_u32(col_file, cols - 1L)
  .study06_write_float32(val_file, values)
  writeLines(c("format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=0", "n_used=0", "n_samples_used=0",
    paste0("n_variants=", marker_count), paste0("nnz=", length(values)),
    "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", row_file), paste0("col_idx_file=", col_file),
    paste0("values_file=", val_file), "row_ptr_type=uint64",
    "col_idx_type=uint32", "values_type=float32", "index_base=0",
    "value=r"), paste0(prefix, ".meta.txt"), useBytes = TRUE)
  writeLines(marker_ids, paste0(prefix, ".markers.txt"), useBytes = TRUE)
  list(prefix = prefix, row_ptr = c(0, cumsum(per_row)),
    col_idx = cols - 1L, values = values, nnz = length(values))
}

.study07_full_csr_glist <- function(Glist, stats, subset_csr, source) {
  out <- Glist
  out$sparseLD <- list(prefix = subset_csr$prefix, chr = stats$chr,
    bed_files = stats$bed_files, cls = stats$cls, af = stats$af,
    marker_names = stats$marker_names, scale = "standardized_genotype",
    source = source, reference_n = stats$n[[1L]], rows = stats$rows,
    orientation_status = "bed_coding_by_construction")
  out$rsidsLD <- lapply(stats$chr, function(cc)
    stats$marker_names[stats$marker_metadata$chromosome_or_file == cc])
  out
}

.study07_operator_bundle <- function(Glist, stats, full_csr, config,
                                     output_dir, effect_matrix = NULL) {
  .study07_source_operator_helpers()
  blocks <- .study06_blocks(stats$marker_names, config$block_size)
  .study06_validate_blocks(blocks, stats$marker_names)
  effect <- if (is.null(effect_matrix)) numeric(length(stats$marker_names))
    else effect_matrix[, 1L]
  unfiltered <- .study06_inspect_operator(Glist, stats, blocks,
    filter = "ridge_fixed", eta = 0, effects = effect)
  runtime <- .study06_write_runtime_csr(unfiltered,
    file.path(output_dir, paste0("runtime_block_", length(stats$marker_names))),
    stats$marker_names)
  runtime_glist <- .study06_runtime_glist(Glist, stats, runtime)
  # The public MT alignment contract recognizes same-BED orientation through
  # the canonical make_sparse_ld source token plus identical BED provenance.
  # This runtime matrix is reconstructed from that same Glist and retains the
  # explicit reference identifier below.
  runtime_glist$sparseLD$source <- "make_sparse_ld"
  runtime_glist$sparseLD$reference_id <-
    "study07_runtime_matched_block_reconstruction"
  runtime_glist$sparseLD$marker_metadata <- stats$marker_metadata
  runtime_inspection <- .study06_inspect_from_csr(runtime$prefix,
    unfiltered$diagonal, blocks, unfiltered)
  scores <- do.call(cbind, lapply(stats$wy, as.numeric))
  if (is.null(effect_matrix)) effect_matrix <- matrix(0,
    nrow(scores), ncol(scores))
  action_vectors <- c(
    setNames(lapply(seq_len(ncol(scores)), function(t) scores[, t]),
      paste0("summary_score_trait", seq_len(ncol(scores)))),
    setNames(lapply(seq_len(ncol(effect_matrix)),
      function(t) effect_matrix[, t]),
      paste0("true_effect_trait", seq_len(ncol(effect_matrix)))))
  equivalence <- .study06_equivalence_gate(runtime_inspection,
    unfiltered, blocks, config, action_vectors = action_vectors)
  dense <- .study06_dense_blocks(unfiltered)
  matrix_action <- do.call(cbind, lapply(seq_len(ncol(scores)),
    function(t) .study06_apply_blocks(dense, scores[, t])))
  direct_action <- do.call(cbind, lapply(seq_len(ncol(scores)),
    function(t) .study06_apply_blocks(
      .study06_dense_blocks(runtime_inspection), scores[, t])))
  matrix_error <- max(abs(matrix_action - direct_action))
  matrix_relative_error <- sqrt(sum((matrix_action - direct_action)^2)) /
    max(sqrt(sum(direct_action^2)), .Machine$double.eps)
  if (matrix_relative_error > config$operator_tolerance$relative)
    stop("MT matrix operator action equivalence failed.", call. = FALSE)
  list(blocks = blocks, unfiltered = unfiltered, runtime = runtime,
    runtime_glist = runtime_glist, runtime_inspection = runtime_inspection,
    equivalence = equivalence, matrix_action_maximum_error = matrix_error,
    matrix_action_relative_error = matrix_relative_error)
}

.study07_operator_summary <- function(bundle, marker_count) {
  x <- bundle$equivalence$summary
  x$marker_count <- marker_count
  x$matrix_matrix_maximum_error <- bundle$matrix_action_maximum_error
  x$matrix_matrix_relative_error <- bundle$matrix_action_relative_error
  x$block_count <- nrow(bundle$blocks)
  x
}
