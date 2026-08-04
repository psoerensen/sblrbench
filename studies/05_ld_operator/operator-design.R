# Study 05-specific LD-operator construction and validation policies.

study05_blocks <- function(marker_ids, block_size = 1000L) {
  if (!is.character(marker_ids) || !length(marker_ids) || anyNA(marker_ids) ||
      anyDuplicated(marker_ids))
    stop("marker_ids must be unique non-missing character identifiers.", call. = FALSE)
  block_size <- as.integer(block_size)
  if (length(block_size) != 1L || is.na(block_size) || block_size < 1L)
    stop("block_size must be one positive integer.", call. = FALSE)
  start <- seq.int(1L, length(marker_ids), by = block_size)
  end <- pmin(start + block_size - 1L, length(marker_ids))
  data.frame(
    block_id = sprintf("block_%04d", seq_along(start)),
    start = start, end = end, size = end - start + 1L,
    first_marker = marker_ids[start], last_marker = marker_ids[end],
    stringsAsFactors = FALSE
  )
}

study05_validate_blocks <- function(blocks, marker_ids) {
  required <- c("block_id", "start", "end", "size", "first_marker", "last_marker")
  if (!is.data.frame(blocks) || !all(required %in% names(blocks)) || !nrow(blocks))
    stop("blocks must contain the complete Study 05 block schema.", call. = FALSE)
  expected <- study05_blocks(marker_ids, max(blocks$size))
  if (!identical(blocks[, required], expected[, required]))
    stop("Blocks do not form the expected contiguous marker partition.", call. = FALSE)
  invisible(blocks)
}

study05_operator_table <- function(spec) {
  eigen <- spec$operators$eigen$proportions
  configurations <- spec$operators$configurations
  data.frame(
    configuration = configurations,
    source_operator = c("packed individual genotypes", "full CSR", "full CSR",
      rep("block CSR", 3L)),
    operator_family = c("BED", "full CSR", "block CSR", rep("block eigen", 3L)),
    representation = c("packed_bed", "csr", "csr", rep("low_rank", 3L)),
    block_size = c(NA_integer_, NA_integer_, rep(spec$operators$block$size, 4L)),
    eigen_policy = c(rep(NA_character_, 3L), rep(spec$operators$eigen$policy, 3L)),
    retained_mass = c(rep(NA_real_, 3L), unname(eigen)),
    residual_rebuild_every = c(rep(NA_integer_, 3L),
      rep(spec$controls$residual_rebuild_every, 3L)),
    stringsAsFactors = FALSE
  )
}

study05_retain_eigensystem <- function(operator, proportion,
    tolerance = 1e-10) {
  if (!is.matrix(operator) || nrow(operator) != ncol(operator) ||
      any(!is.finite(operator)))
    stop("operator must be a finite square matrix.", call. = FALSE)
  if (!is.numeric(proportion) || length(proportion) != 1L ||
      !is.finite(proportion) || proportion <= 0 || proportion >= 1)
    stop("proportion must be strictly between zero and one.", call. = FALSE)
  symmetry_error <- max(abs(operator - t(operator)))
  if (symmetry_error > tolerance)
    stop("operator is not symmetric within tolerance.", call. = FALSE)
  decomposition <- eigen((operator + t(operator)) / 2, symmetric = TRUE)
  positive <- which(decomposition$values > tolerance)
  if (!length(positive)) stop("operator has no positive eigenvalues.", call. = FALSE)
  mass <- cumsum(decomposition$values[positive]) / sum(decomposition$values[positive])
  rank <- match(TRUE, mass >= proportion)
  keep <- positive[seq_len(rank)]
  factor <- sqrt(decomposition$values[keep]) * t(decomposition$vectors[, keep,
    drop = FALSE])
  list(
    factor = factor,
    values = decomposition$values,
    retained_rank = length(keep),
    positive_rank = length(positive),
    retained_mass = mass[rank],
    reconstructed = crossprod(factor)
  )
}

study05_low_rank_identity <- function(factor, transformed_score, beta, yy) {
  if (!is.matrix(factor) || ncol(factor) != length(beta) ||
      length(transformed_score) != nrow(factor))
    stop("Low-rank identity inputs have incompatible dimensions.", call. = FALSE)
  residual <- as.numeric(transformed_score - factor %*% beta)
  projected_score <- as.numeric(crossprod(factor, transformed_score))
  quadratic <- sum((factor %*% beta)^2)
  sse_residual <- yy - sum(transformed_score^2) + sum(residual^2)
  sse_quadratic <- yy - 2 * sum(beta * projected_score) + quadratic
  data.frame(
    projected_sse_residual = sse_residual,
    projected_sse_quadratic = sse_quadratic,
    projected_sse_identity_error = abs(sse_residual - sse_quadratic),
    corrected_score_error = max(abs(
      as.numeric(crossprod(factor, residual)) + colSums(factor^2) * beta -
        (projected_score - as.numeric(crossprod(factor) %*% beta) +
          colSums(factor^2) * beta)
    )),
    stringsAsFactors = FALSE
  )
}

study05_sbayesr_design <- function(spec) {
  x <- spec$components$sbayesr_ld_sensitivity
  data.frame(
    component = x$component,
    marker_count = x$marker_window$count,
    first_source_row = min(x$marker_window$source_rows),
    last_source_row = max(x$marker_window$source_rows),
    simulation_seed = x$simulation_seed,
    retained_mass = x$retained_mass,
    retained_total_rank = x$retained_total_rank,
    stringsAsFactors = FALSE
  )
}

# Resume the completed benchmark through its immutable compact result contract.
# Native fits remain historical local caches; a future explicit scientific rerun
# can inject a coordinate runner through `sblrbench.ld_operator_runner` without
# changing the shared task dispatch or this scientific specification.
study05_reference_operator_runner <- function(spec, profile, coordinates,
    paths, resume = TRUE) {
  if (!isTRUE(resume))
    stop("Study 05 has no implicit refit mode. Supply an explicit scientific ",
      "coordinate runner for a deliberate rerun.", call. = FALSE)
  capsule <- benchmark_spec_path(spec, spec$frozen_capsule)
  read_table <- function(name) {
    path <- file.path(capsule, name)
    if (!file.exists(path)) stop("Frozen Study 05 table is missing: ", path,
      call. = FALSE)
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  status <- read_table("fit_status.csv")
  expected <- coordinates[c("scenario", "replicate", "configuration")]
  names(expected)[1L] <- "architecture"
  key <- function(x) paste(x$architecture, x$replicate, x$configuration,
    sep = "\r")
  selected <- match(key(expected), key(status))
  if (anyNA(selected))
    stop("Frozen Study 05 fit coordinates do not contain the selected profile.",
      call. = FALSE)
  status <- status[selected, , drop = FALSE]
  status$status <- ifelse(status$status == "ok", "ok", status$status)
  status$reused <- TRUE
  runtime <- read_table("computational_summary.csv")
  runtime$study <- spec$study
  runtime$scenario <- runtime$architecture
  runtime$native_method <- runtime$method
  runtime$method <- runtime$configuration
  runtime$status <- "ok"
  runtime$reason <- ""
  runtime$reused <- TRUE
  tables <- list(
    status = status,
    operator_summary = read_table("deterministic_identity_summary.csv"),
    operator_comparisons = read_table("deterministic_identity_summary.csv"),
    eigenvalue_summary = read_table("low_rank_block_diagnostics.csv"),
    convergence = read_table("convergence_validation_summary.csv"),
    recovery_metrics = read_table("parameter_recovery_summary.csv"),
    runtime = runtime,
    sbayesr_fit_status = read_table("sbayesr_fit_status.csv"),
    sbayesr_scheduler = read_table("sbayesr_scheduler_comparison.csv"),
    sbayesr_variance = read_table("sbayesr_variance_summary.csv"),
    sbayesr_conditionals = read_table("sbayesr_conditional_marker_summary.csv"),
    sbayesr_quadratics = read_table("sbayesr_operator_quadratic_checks.csv"),
    sbayesr_recovery = read_table("sbayesr_effect_recovery.csv"),
    sbayesr_eigenvalues = read_table("sbayesr_eigenvalue_summary.csv")
  )
  output_files <- list()
  filenames <- c(status = "fit_status", operator_summary = "operator_summary",
    operator_comparisons = "operator_pairwise_comparison",
    eigenvalue_summary = "eigenvalue_summary", convergence = "convergence",
    recovery_metrics = "recovery_metrics", runtime = "runtime",
    sbayesr_fit_status = "sbayesr_fit_status",
    sbayesr_scheduler = "sbayesr_scheduler_comparison",
    sbayesr_variance = "sbayesr_variance_summary",
    sbayesr_conditionals = "sbayesr_conditional_marker_summary",
    sbayesr_quadratics = "sbayesr_operator_quadratic_checks",
    sbayesr_recovery = "sbayesr_effect_recovery",
    sbayesr_eigenvalues = "sbayesr_eigenvalue_summary")
  for (name in names(tables)) {
    path <- file.path(paths$tables, paste0(filenames[[name]], ".csv"))
    utils::write.csv(tables[[name]], path, row.names = FALSE, na = "")
    output_files[[filenames[[name]]]] <- path
  }
  manifest <- list(schema_version = 1L, study = spec$study, task = spec$task,
    profile = profile, coordinate_count = nrow(coordinates),
    execution = "frozen_reference_reuse", sampler_calls = 0L,
    source_capsule = spec$frozen_capsule,
    checkpoint_schema = "sblrbench-semantic-v2")
  jsonlite::write_json(manifest, paths$manifest, auto_unbox = TRUE, pretty = TRUE)
  writeLines(benchmark_session_information(), paths$session_info)
  c(tables, list(paths = output_files))
}

# Main operator construction -------------------------------------------------
.study05_unpack_triangle <- function(packed, size) {
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

.study05_inspect_operator <- function(Glist, stats, blocks,
                                      filter = "ridge_fixed",
                                      tau = 0, eta = 0,
                                      effects = NULL) {
  study06_validate_blocks(blocks, stats$marker_names)
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

.study05_dense_blocks <- function(inspect) {
  Map(.study05_unpack_triangle, inspect$packed_upper_triangle,
    inspect$block_size)
}

.study05_pack_triangle <- function(A) {
  unlist(lapply(seq_len(nrow(A)), function(i) A[i, i:ncol(A)]),
    use.names = FALSE)
}

.study05_apply_blocks <- function(blocks, x) {
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

.study05_apply_csr_crossproduct <- function(csr, diagonal, x) {
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

.study05_action_comparison <- function(reference_action, candidate_action,
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

.study05_write_u32 <- function(path, x) {
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

.study05_write_float32 <- function(path, x) {
  if (any(!is.finite(x))) stop("Non-finite CSR values.", call. = FALSE)
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.numeric(x), con, size = 4L, endian = "little")
  invisible(path)
}

.study05_write_runtime_csr <- function(inspect, prefix, marker_ids,
                                       zero_tolerance = 0) {
  dense <- .study05_dense_blocks(inspect)
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
  .study05_write_u32(col_file, col_idx)
  .study05_write_float32(val_file, values)
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

.study05_runtime_glist <- function(Glist, stats, csr) {
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

.study05_inspect_from_csr <- function(prefix, runtime_diagonal, blocks,
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
  out$packed_upper_triangle <- lapply(dense, .study05_pack_triangle)
  out$diagonal <- runtime_diagonal
  out$filter <- list(mode = "runtime_matched_block_csr",
    tau = NA_real_, eta = NA_real_, mu_floor = NA_real_)
  out
}


# Operator validation --------------------------------------------------------
.study05_operator_metrics <- function(reference, candidate, blocks,
                                      config, operator_id,
                                      reference_id = "runtime_matched_block_csr",
                                      action_vectors = list()) {
  ref <- .study05_dense_blocks(reference)
  alt <- .study05_dense_blocks(candidate)
  if (length(ref) != length(alt) || length(ref) != nrow(blocks))
    stop("Operator block counts differ.", call. = FALSE)
  set.seed(config$seeds$operator_probe)
  probes <- replicate(config$operator_probe_count,
    stats::rnorm(reference$marker_count), simplify = FALSE)
  rows <- lapply(seq_along(ref), function(k) {
    A <- ref[[k]]; B <- alt[[k]]
    delta <- B - A
    denom <- pmax(abs(A), .Machine$double.eps^0.5)
    eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    eig_alt <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
    positive <- sum(pmax(eig, 0))
    data.frame(
      operator_id = operator_id, reference_operator_id = reference_id,
      block_id = blocks$block_id[k], block_size = nrow(A),
      source_matrix_definition =
        "same installed builder packed runtime cross-product",
      filtering_policy = candidate$filter$mode,
      filtering_parameter = if (candidate$filter$mode == "hard_truncate")
        candidate$filter$tau else candidate$filter$eta,
      original_rank = qr(A)$rank,
      retained_rank = if (!is.null(candidate$diagnostics$n_kept))
        candidate$diagnostics$n_kept[k] else qr(B)$rank,
      retained_rank_proportion = if (!is.null(candidate$diagnostics$n_kept))
        candidate$diagnostics$n_kept[k] / nrow(B) else qr(B)$rank / nrow(B),
      minimum_eigenvalue = min(eig_alt),
      maximum_eigenvalue = max(eig_alt),
      negative_eigenvalue_count = sum(eig_alt < -1e-10),
      positive_eigenvalue_mass = positive,
      retained_positive_eigenvalue_mass = sum(pmax(eig_alt, 0)),
      retained_mass_proportion = if (positive > 0)
        sum(pmax(eig_alt, 0)) / positive else NA_real_,
      trace = sum(diag(B)), trace_difference = sum(diag(B)) - sum(diag(A)),
      frobenius_norm = sqrt(sum(B^2)),
      frobenius_relative_error = sqrt(sum(delta^2)) /
        max(sqrt(sum(A^2)), .Machine$double.eps),
      off_diagonal_frobenius_norm = sqrt(sum((B - diag(diag(B)))^2)),
      off_diagonal_norm_retained = sqrt(sum((B - diag(diag(B)))^2)) /
        max(sqrt(sum((A - diag(diag(A)))^2)), .Machine$double.eps),
      condition_number = {
        pos <- eig_alt[eig_alt > max(eig_alt) * .Machine$double.eps^0.5]
        if (length(pos)) max(pos) / min(pos) else Inf
      },
      runtime_diagonal_minimum = min(diag(B)),
      runtime_diagonal_maximum = max(diag(B)),
      runtime_diagonal_mean = mean(diag(B)),
      reconstruction_maximum_absolute_error = max(abs(delta)),
      reconstruction_relative_error = max(abs(delta) / denom),
      matrix_vector_maximum_error = max(vapply(probes, function(x) {
        idx <- blocks$start[k]:blocks$end[k]
        max(abs(as.numeric(B %*% x[idx]) - as.numeric(A %*% x[idx])))
      }, numeric(1))),
      quadratic_form_maximum_error = max(vapply(probes, function(x) {
        idx <- blocks$start[k]:blocks$end[k]
        abs(drop(crossprod(x[idx], B %*% x[idx])) -
          drop(crossprod(x[idx], A %*% x[idx])))
      }, numeric(1))),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (length(action_vectors)) {
    for (name in names(action_vectors)) {
      x <- as.numeric(action_vectors[[name]])
      if (length(x) != reference$marker_count || any(!is.finite(x)))
        stop("Invalid operator action vector: ", name, call. = FALSE)
      yr <- .study05_apply_blocks(ref, x)
      ya <- .study05_apply_blocks(alt, x)
      rel <- sqrt(sum((ya - yr)^2)) / max(sqrt(sum(yr^2)),
        .Machine$double.eps)
      cosine <- sum(yr * ya) / max(sqrt(sum(yr^2) * sum(ya^2)),
        .Machine$double.eps)
      out[[paste0(name, "_action_relative_error")]] <- rel
      out[[paste0(name, "_action_cosine_similarity")]] <- cosine
    }
  }
  out
}

.study05_equivalence_gate <- function(reference, candidate, blocks, config,
                                      action_vectors = list()) {
  metrics <- .study05_operator_metrics(reference, candidate, blocks, config,
    operator_id = "block_eigen_unfiltered",
    action_vectors = action_vectors)
  tol <- config$operator_tolerance
  mapping_ok <- identical(reference$marker_count, candidate$marker_count) &&
    identical(reference$block_start, candidate$block_start) &&
    identical(reference$block_size, candidate$block_size) &&
    identical(reference$block_of, candidate$block_of) &&
    identical(reference$local_of, candidate$local_of)
  diag_error <- max(abs(reference$diagonal - candidate$diagonal))
  pass <- mapping_ok && all(is.finite(candidate$diagonal)) &&
    all(candidate$diagonal > 0) &&
    diag_error <= tol$absolute &&
    max(metrics$reconstruction_maximum_absolute_error) <= tol$absolute &&
    max(metrics$reconstruction_relative_error) <= tol$relative &&
    max(metrics$matrix_vector_maximum_error) <= tol$product_absolute &&
    max(metrics$quadratic_form_maximum_error) <= tol$quadratic_absolute
  summary <- data.frame(
    gate = "runtime_matched_block_csr_vs_unfiltered_block_eigen",
    pass = pass, marker_mapping_pass = mapping_ok,
    block_coverage_pass = blocks$start[1L] == 1L &&
      blocks$end[nrow(blocks)] == reference$marker_count &&
      identical(blocks$start[-1L],
        blocks$end[-nrow(blocks)] + 1L),
    diagonal_maximum_absolute_error = diag_error,
    reconstruction_maximum_absolute_error =
      max(metrics$reconstruction_maximum_absolute_error),
    reconstruction_maximum_relative_error =
      max(metrics$reconstruction_relative_error),
    matrix_vector_maximum_error =
      max(metrics$matrix_vector_maximum_error),
    quadratic_form_maximum_error =
      max(metrics$quadratic_form_maximum_error),
    matrix_error_block = metrics$block_id[
      which.max(metrics$reconstruction_maximum_absolute_error)],
    relative_error_block = metrics$block_id[
      which.max(metrics$reconstruction_relative_error)],
    matrix_vector_error_block = metrics$block_id[
      which.max(metrics$matrix_vector_maximum_error)],
    quadratic_form_error_block = metrics$block_id[
      which.max(metrics$quadratic_form_maximum_error)],
    float32_reconstruction_contribution =
      max(metrics$reconstruction_maximum_absolute_error),
    absolute_tolerance = tol$absolute,
    relative_tolerance = tol$relative,
    product_absolute_tolerance = tol$product_absolute,
    quadratic_absolute_tolerance = tol$quadratic_absolute,
    rationale = paste("float32 packed runtime storage with",
      "double-precision validation products"),
    stringsAsFactors = FALSE)
  if (!pass) {
    first <- metrics[which.max(metrics$reconstruction_maximum_absolute_error), ]
    stop("Study 05 operator-equivalence gate failed at ",
      first$block_id, "; maximum absolute error = ",
      signif(first$reconstruction_maximum_absolute_error, 6),
      call. = FALSE)
  }
  list(summary = summary, block_metrics = metrics)
}

.study05_filter_dense_matrix <- function(A,
    mode = c("ridge_fixed", "hard_truncate", "ridge_lw"),
    tau = 0, eta = 0, lw_shrinkage = NULL) {
  mode <- match.arg(mode)
  A <- as.matrix(A)
  if (nrow(A) != ncol(A) || any(!is.finite(A)) ||
      max(abs(A - t(A))) > 1e-10 || any(diag(A) <= 0))
    stop("Synthetic Study 05 operator must be finite, symmetric, and positive-diagonal.",
      call. = FALSE)
  d <- diag(A)
  if (mode == "hard_truncate") {
    C <- A / sqrt(outer(d, d))
    e <- eigen((C + t(C)) / 2, symmetric = TRUE)
    effective <- max(tau, 0.01)
    keep <- e$values >= effective
    if (!any(keep)) keep[which.max(e$values)] <- TRUE
    V <- e$vectors[, keep, drop = FALSE]
    B <- (sqrt(d) * V) %*% diag(e$values[keep], nrow = sum(keep)) %*%
      t(sqrt(d) * V)
    return(list(matrix = B, eigenvalues = e$values,
      retained = e$values[keep], retained_rank = sum(keep),
      effective_threshold = effective,
      shrinkage_weight = 1 - sum(keep) / length(keep)))
  }
  a <- if (mode == "ridge_fixed") eta / (1 + eta) else lw_shrinkage
  if (length(a) != 1L || !is.finite(a) || a < 0 || a > 1)
    stop("Synthetic ridge shrinkage must lie in [0, 1].", call. = FALSE)
  B <- (1 - a) * A
  diag(B) <- d
  list(matrix = B, retained_rank = qr(B)$rank,
    effective_threshold = NA_real_, shrinkage_weight = a)
}

.study05_known_spectrum_matrix <- function(values) {
  p <- length(values)
  raw <- outer(seq_len(p), seq_len(p), function(i, j)
    sin(i * (j + 0.5)) + cos((i + 0.25) * j))
  Q <- qr.Q(qr(raw))
  Q %*% diag(values, nrow = p) %*% t(Q)
}

.study05_filter_symmetric_spectrum <- function(A, tau) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  effective <- max(tau, 0.01)
  keep <- e$values >= effective
  if (!any(keep)) keep[which.max(e$values)] <- TRUE
  V <- e$vectors[, keep, drop = FALSE]
  list(matrix = V %*% diag(e$values[keep], nrow = sum(keep)) %*% t(V),
    eigenvalues = e$values, retained = e$values[keep],
    retained_rank = sum(keep), effective_threshold = effective)
}

.study05_synthetic_filter_validation <- function(tolerance = 1e-10) {
  spectra <- list(
    well_conditioned = c(5, 3, 2, 1, 0.5),
    controlled_small = c(5, 2, 1, 0.05, 0.005),
    near_singular = c(5, 2, 1, 0.05, 1e-8))
  rows <- list()
  for (id in names(spectra)) {
    A <- .study05_known_spectrum_matrix(spectra[[id]])
    zero <- .study05_filter_dense_matrix(A, "ridge_fixed", eta = 0)
    for (tau in c(0.001, 0.01, 0.10)) {
      hard <- .study05_filter_symmetric_spectrum(A, tau)
      e <- eigen(A, symmetric = TRUE)
      keep <- e$values >= max(tau, 0.01)
      V <- e$vectors[, keep, drop = FALSE]
      direct <- V %*% diag(e$values[keep], nrow = sum(keep)) %*% t(V)
      x <- seq_len(nrow(A)) / nrow(A)
      rows[[length(rows) + 1L]] <- data.frame(
        matrix_id = id, requested_threshold = tau,
        effective_threshold = max(tau, 0.01),
        original_rank = nrow(A), retained_rank = hard$retained_rank,
        retained_mass_proportion = sum(hard$retained) /
          sum(pmax(hard$eigenvalues, 0)),
        zero_ridge_maximum_error = max(abs(zero$matrix - A)),
        filtered_reconstruction_maximum_error = max(abs(hard$matrix - direct)),
        matrix_vector_maximum_error = max(abs(hard$matrix %*% x - direct %*% x)),
        quadratic_form_error = abs(drop(crossprod(x, hard$matrix %*% x)) -
          drop(crossprod(x, direct %*% x))),
        minimum_diagonal = min(diag(hard$matrix)),
        pass = max(abs(zero$matrix - A)) <= tolerance &&
          max(abs(hard$matrix - direct)) <= tolerance &&
          all(diag(hard$matrix) > 0), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (!all(out$pass)) stop("Synthetic Study 05 filter validation failed.",
    call. = FALSE)
  out
}


# SBayesR LD-sensitivity helpers ---------------------------------------------
# Shared helpers for Study 05 SBayesR spectral/operator diagnostics.

study05_sbayesr_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
study05_sbayesr_hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

study05_sbayesr_atomic_rds <- function(x, path) {
  tmp <- tempfile(".tmp-", dirname(path), ".rds")
  saveRDS(x, tmp, version = 3)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

study05_sbayesr_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".tmp-", dirname(path), ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

study05_sbayesr_unpack_triangle <- function(packed, size) {
  stopifnot(length(packed) == size * (size + 1L) / 2L)
  out <- matrix(0, size, size)
  k <- 1L
  for (i in seq_len(size)) for (j in i:size) {
    out[i, j] <- out[j, i] <- packed[k]
    k <- k + 1L
  }
  out
}

study05_sbayesr_dense_csr <- function(csr) {
  out <- diag(1, csr$nrow)
  ii <- rep.int(seq_len(csr$nrow), diff(as.integer(csr$row_ptr)))
  jj <- as.integer(csr$col_idx)
  out[cbind(ii, jj)] <- csr$values
  out[cbind(jj, ii)] <- csr$values
  out
}

study05_sbayesr_block_definitions <- function(marker_table, target_size = 250L) {
  m <- nrow(marker_table)
  chromosome_starts <- c(1L, which(marker_table$chromosome[-1L] !=
    marker_table$chromosome[-m]) + 1L)
  fixed_starts <- seq.int(1L, m, by = target_size)
  starts <- sort(unique(c(chromosome_starts, fixed_starts)))
  ends <- c(starts[-1L] - 1L, m)
  data.frame(
    block_id = sprintf("block_%02d", seq_along(starts)),
    start = starts, end = ends, size = ends - starts + 1L,
    chromosome_first = marker_table$chromosome[starts],
    chromosome_last = marker_table$chromosome[ends],
    first_marker = marker_table$marker_id[starts],
    last_marker = marker_table$marker_id[ends],
    stringsAsFactors = FALSE)
}

study05_sbayesr_mask_blocks <- function(A, blocks) {
  group <- integer(nrow(A))
  for (i in seq_len(nrow(blocks))) group[blocks$start[i]:blocks$end[i]] <- i
  out <- A
  out[outer(group, group, `!=`)] <- 0
  out
}

study05_sbayesr_retained_operator <- function(A, blocks, eigen_prop = 0.995,
                                   tolerance = 1e-10) {
  out <- matrix(0, nrow(A), ncol(A))
  diagnostics <- vector("list", nrow(blocks))
  for (i in seq_len(nrow(blocks))) {
    idx <- blocks$start[i]:blocks$end[i]
    Ab <- A[idx, idx, drop = FALSE]
    d <- sqrt(diag(Ab))
    C <- Ab / tcrossprod(d)
    C <- (C + t(C)) / 2
    eg <- eigen(C, symmetric = TRUE)
    positive <- which(eg$values > tolerance)
    ord <- positive[order(eg$values[positive], decreasing = TRUE)]
    positive_mass <- sum(eg$values[positive])
    running <- cumsum(eg$values[ord]) / positive_mass
    retained_n <- which(running > eigen_prop)[1L]
    keep <- ord[seq_len(retained_n)]
    Cr <- eg$vectors[, keep, drop = FALSE] %*%
      (eg$values[keep] * t(eg$vectors[, keep, drop = FALSE]))
    out[idx, idx] <- Cr * tcrossprod(d)
    diagnostics[[i]] <- data.frame(
      block_id = blocks$block_id[i], start = blocks$start[i],
      end = blocks$end[i], size = length(idx), positive_rank = length(positive),
      retained_rank = length(keep), discarded_rank = length(idx) - length(keep),
      positive_eigenvalue_mass = positive_mass,
      retained_eigenvalue_mass = sum(eg$values[keep]),
      retained_mass_fraction = sum(eg$values[keep]) / positive_mass,
      minimum_retained_eigenvalue = min(eg$values[keep]),
      maximum_omitted_eigenvalue = if (length(setdiff(positive, keep)))
        max(eg$values[setdiff(positive, keep)]) else 0,
      negative_eigenvalue_count = sum(eg$values < -tolerance),
      negative_eigenvalue_mass = sum(abs(eg$values[eg$values < -tolerance])))
  }
  list(matrix = out, diagnostics = do.call(rbind, diagnostics))
}

study05_sbayesr_operator_metrics <- function(A, reference, id, retained_rank = NA_integer_) {
  ev <- eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values
  err <- A - reference
  data.frame(
    operator = id,
    relative_frobenius_error = sqrt(sum(err^2) / sum(reference^2)),
    maximum_absolute_error = max(abs(err)),
    diagonal_maximum_error = max(abs(diag(err))),
    symmetry_error = max(abs(A - t(A))),
    retained_rank = retained_rank,
    trace = sum(diag(A)), trace_difference = sum(diag(A)) - sum(diag(reference)),
    minimum_eigenvalue = min(ev), maximum_eigenvalue = max(ev),
    negative_eigenvalue_count = sum(ev < -1e-8),
    negative_eigenvalue_mass = sum(abs(ev[ev < -1e-8])),
    stringsAsFactors = FALSE)
}

study05_sbayesr_spectrum <- function(A, id) {
  values <- sort(eigen((A + t(A)) / 2, symmetric = TRUE,
    only.values = TRUE)$values, decreasing = TRUE)
  positive <- values[values > 1e-8]
  abs_sum <- sum(abs(values))
  p <- abs(values) / abs_sum
  effective <- exp(-sum(p[p > 0] * log(p[p > 0])))
  condition <- if (length(positive)) max(positive) / min(positive) else NA_real_
  summary <- data.frame(
    operator = id, minimum = min(values), maximum = max(values),
    below_negative_tolerance = sum(values < -1e-8),
    near_zero = sum(values >= -1e-8 & values <= 1e-8),
    positive = sum(values > 1e-8), sum_negative = sum(values[values < -1e-8]),
    absolute_negative_mass = sum(abs(values[values < -1e-8])),
    negative_component_frobenius = sqrt(sum(values[values < -1e-8]^2)),
    trace = sum(values), effective_rank = effective,
    positive_condition_number = condition,
    q000 = min(values), q001 = unname(quantile(values, .001)),
    q010 = unname(quantile(values, .01)), q050 = unname(quantile(values, .05)),
    q500 = unname(quantile(values, .5)), q950 = unname(quantile(values, .95)),
    q990 = unname(quantile(values, .99)), q999 = unname(quantile(values, .999)),
    q1000 = max(values))
  spectrum <- data.frame(operator = id, rank_descending = seq_along(values),
    eigenvalue = values, cumulative_absolute_proportion = cumsum(abs(values)) / abs_sum)
  list(summary = summary, spectrum = spectrum, values = values)
}

study05_sbayesr_trace_table <- function(fit, id, label, draws = 1000L, chains = 4L) {
  b <- fit$convergence_traces
  q <- b$quantities
  rows <- vector("list", dim(b$values)[3L])
  for (j in seq_len(dim(b$values)[3L])) {
    group <- as.character(q$group[j])
    component <- if ("component_name" %in% names(q))
      as.character(q$component_name[j]) else NA_character_
    quantity <- if (group == "component_pi") paste0("pi_", component) else group
    rows[[j]] <- data.frame(variant = id, label = label,
      iteration = rep(seq_len(draws), chains),
      chain = rep(seq_len(chains), each = draws), quantity = quantity,
      value = as.vector(b$values[, , j]))
  }
  out <- do.call(rbind, rows)
  make_matrix <- function(quantity) matrix(out$value[out$quantity == quantity], draws, chains)
  vg <- make_matrix("vgs")
  ve <- make_matrix("ves")
  out <- rbind(out, data.frame(variant = id, label = label,
    iteration = rep(seq_len(draws), chains), chain = rep(seq_len(chains), each = draws),
    quantity = "heritability", value = as.vector(vg / (vg + ve))))
  p0 <- make_matrix("pi_component_0")
  rbind(out, data.frame(variant = id, label = label,
    iteration = rep(seq_len(draws), chains), chain = rep(seq_len(chains), each = draws),
    quantity = "active_probability", value = as.vector(1 - p0)))
}

study05_sbayesr_summarise_traces <- function(traces) {
  groups <- split(traces, interaction(traces$variant, traces$quantity, drop = TRUE))
  variance <- do.call(rbind, lapply(groups, function(x) data.frame(
    variant = x$variant[1L], label = x$label[1L], quantity = x$quantity[1L],
    mean = mean(x$value), sd = stats::sd(x$value),
    lower_025 = unname(stats::quantile(x$value, .025)),
    upper_975 = unname(stats::quantile(x$value, .975)))))
  convergence <- do.call(rbind, lapply(groups, function(x) {
    m <- matrix(x$value, length(unique(x$iteration)), length(unique(x$chain)))
    mcse <- posterior::mcse_mean(m)
    data.frame(variant = x$variant[1L], label = x$label[1L],
      quantity = x$quantity[1L], rhat = posterior::rhat(m),
      bulk_ess = posterior::ess_bulk(m), tail_ess = posterior::ess_tail(m),
      mcse = mcse, relative_mcse = mcse / stats::sd(x$value))
  }))
  list(variance = variance, convergence = convergence)
}
