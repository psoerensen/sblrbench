.study04_extract_chain_draws <- function(fit, architecture, method, expected_chains = 4L) {
  bundle <- fit$convergence_traces
  if (is.null(bundle$values) || length(dim(bundle$values)) != 3L) stop("Native convergence trace bundle is unavailable.", call. = FALSE)
  values <- bundle$values; q <- bundle$quantities
  if (dim(values)[2L] != expected_chains) stop("Unexpected number of identifiable chains.", call. = FALSE)
  groups <- as.character(q$group); required <- c("vbs", "vgs", "ves")
  idx <- match(required, groups); if (anyNA(idx)) stop("Required scalar traces are absent.", call. = FALSE)
  out <- expand.grid(raw_iteration = seq_len(dim(values)[1L]), chain = seq_len(dim(values)[2L]))
  out$architecture <- architecture; out$method <- method
  for (j in seq_along(required)) out[[c("effect_variance", "genetic_variance", "residual_variance")[j]]] <- as.vector(values[, , idx[j]])
  out$heritability <- out$genetic_variance / (out$genetic_variance + out$residual_variance)
  out <- out[c("architecture", "method", "chain", "raw_iteration", "effect_variance",
    "genetic_variance", "residual_variance", "heritability")]
  if (any(!is.finite(as.matrix(out[5:8]))) || any(as.matrix(out[5:7]) < 0) ||
      any(out$heritability < 0 | out$heritability > 1) || anyDuplicated(out[c("chain", "raw_iteration")]))
    stop("Invalid scalar chain draws.", call. = FALSE)
  out
}

.study04_window <- function(draws, burnin, retained, expected_chains = 4L) {
  if (anyDuplicated(draws[c("chain", "raw_iteration")]) ||
      !identical(sort(unique(draws$chain)), seq_len(expected_chains)))
    stop("Invalid chain identity.", call. = FALSE)
  z <- draws[draws$raw_iteration > burnin & draws$raw_iteration <= burnin + retained, ]
  if (any(table(z$chain) != retained)) stop("Insufficient chain length for checkpoint.", call. = FALSE)
  z
}
