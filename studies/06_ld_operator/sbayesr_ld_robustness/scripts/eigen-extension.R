# Shared helpers for supplemental Study 06 spectral/operator diagnostics.

s06rob_hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
s06rob_hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

s06rob_atomic_rds <- function(x, path) {
  tmp <- tempfile(".tmp-", dirname(path), ".rds")
  saveRDS(x, tmp, version = 3)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

s06rob_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".tmp-", dirname(path), ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

s06rob_unpack_triangle <- function(packed, size) {
  stopifnot(length(packed) == size * (size + 1L) / 2L)
  out <- matrix(0, size, size)
  k <- 1L
  for (i in seq_len(size)) for (j in i:size) {
    out[i, j] <- out[j, i] <- packed[k]
    k <- k + 1L
  }
  out
}

s06rob_dense_csr <- function(csr) {
  out <- diag(1, csr$nrow)
  ii <- rep.int(seq_len(csr$nrow), diff(as.integer(csr$row_ptr)))
  jj <- as.integer(csr$col_idx)
  out[cbind(ii, jj)] <- csr$values
  out[cbind(jj, ii)] <- csr$values
  out
}

s06rob_block_definitions <- function(marker_table, target_size = 250L) {
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

s06rob_mask_blocks <- function(A, blocks) {
  group <- integer(nrow(A))
  for (i in seq_len(nrow(blocks))) group[blocks$start[i]:blocks$end[i]] <- i
  out <- A
  out[outer(group, group, `!=`)] <- 0
  out
}

s06rob_retained_operator <- function(A, blocks, eigen_prop = 0.995,
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

s06rob_operator_metrics <- function(A, reference, id, retained_rank = NA_integer_) {
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

s06rob_spectrum <- function(A, id) {
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

s06rob_trace_table <- function(fit, id, label, draws = 1000L, chains = 4L) {
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

s06rob_summarise_traces <- function(traces) {
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
