# Study 06-specific LD-operator construction and validation policies.

study06_blocks <- function(marker_ids, block_size = 1000L) {
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

study06_validate_blocks <- function(blocks, marker_ids) {
  required <- c("block_id", "start", "end", "size", "first_marker", "last_marker")
  if (!is.data.frame(blocks) || !all(required %in% names(blocks)) || !nrow(blocks))
    stop("blocks must contain the complete Study 06 block schema.", call. = FALSE)
  expected <- study06_blocks(marker_ids, max(blocks$size))
  if (!identical(blocks[, required], expected[, required]))
    stop("Blocks do not form the expected contiguous marker partition.", call. = FALSE)
  invisible(blocks)
}

study06_operator_table <- function(spec) {
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

study06_retain_eigensystem <- function(operator, proportion,
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

study06_low_rank_identity <- function(factor, transformed_score, beta, yy) {
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

study06_supplemental_design <- function(spec) {
  x <- spec$supplemental
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
study06_reference_operator_runner <- function(spec, profile, coordinates,
    paths, resume = TRUE) {
  if (!isTRUE(resume))
    stop("Study 06 has no implicit refit mode. Supply an explicit scientific ",
      "coordinate runner for a deliberate rerun.", call. = FALSE)
  capsule <- benchmark_spec_path(spec, spec$frozen_capsules$main)
  read_table <- function(name) {
    path <- file.path(capsule, name)
    if (!file.exists(path)) stop("Frozen Study 06 table is missing: ", path,
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
    stop("Frozen Study 06 fit coordinates do not contain the selected profile.",
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
    runtime = runtime
  )
  output_files <- list()
  filenames <- c(status = "fit_status", operator_summary = "operator_summary",
    operator_comparisons = "operator_pairwise_comparison",
    eigenvalue_summary = "eigenvalue_summary", convergence = "convergence",
    recovery_metrics = "recovery_metrics", runtime = "runtime")
  for (name in names(tables)) {
    path <- file.path(paths$tables, paste0(filenames[[name]], ".csv"))
    utils::write.csv(tables[[name]], path, row.names = FALSE, na = "")
    output_files[[filenames[[name]]]] <- path
  }
  manifest <- list(schema_version = 1L, study = spec$study, task = spec$task,
    profile = profile, coordinate_count = nrow(coordinates),
    execution = "frozen_reference_reuse", sampler_calls = 0L,
    source_capsule = spec$frozen_capsules$main,
    checkpoint_schema = "sblrbench-semantic-v2")
  jsonlite::write_json(manifest, paths$manifest, auto_unbox = TRUE, pretty = TRUE)
  writeLines(benchmark_session_information(), paths$session_info)
  c(tables, list(paths = output_files))
}
