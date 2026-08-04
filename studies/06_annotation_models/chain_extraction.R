.study06_quantity_id <- function(q) {
  out <- as.character(q$parameter_name)
  missing <- is.na(out) | !nzchar(out)
  if ("group" %in% names(q)) out[missing] <- as.character(q$group[missing])
  alpha <- !is.na(out) & out == "alpha"
  out[alpha] <- paste0("alpha:", q$annotation_name[alpha], ":", q$stick_name[alpha])
  sigma <- !is.na(out) & out == "sigmaSqAlpha"
  out[sigma] <- paste0("sigmaSqAlpha:", q$stick_name[sigma])
  component <- out %in% c("component_pi", "pattern_pi", "pi_active")
  out[component] <- paste0(out[component], ":",
    ifelse(!is.na(q$component_name[component]), q$component_name[component],
      ifelse(!is.na(q$pattern_name[component]), q$pattern_name[component], "global")))
  out
}

.study06_extract_chain_draws <- function(run, A) {
  if (run$status != "ok") stop("Cannot extract a failed fit.", call. = FALSE)
  bundle <- run$fit$convergence_traces
  values <- bundle$values; q <- bundle$quantities
  if (length(dim(values)) != 3L || dim(values)[2L] != 4L ||
      nrow(q) != dim(values)[3L])
    stop("Invalid Study 06 convergence trace bundle.", call. = FALSE)
  quantity_id <- .study06_quantity_id(q)
  base <- expand.grid(iteration = seq_len(dim(values)[1L]),
    chain = seq_len(dim(values)[2L]))
  out <- do.call(rbind, lapply(seq_along(quantity_id), function(j)
    data.frame(base, quantity = quantity_id[j],
      parameter_name = q$parameter_name[j],
      annotation_name = q$annotation_name[j],
      stick_name = q$stick_name[j],
      value = as.vector(values[, , j]), stringsAsFactors = FALSE)))
  vg <- which(q$group == "vgs")
  ve <- which(q$group == "ves")
  if (length(vg) != 1L || length(ve) != 1L)
    stop("Core variance traces are absent.", call. = FALSE)
  h2 <- values[, , vg] / (values[, , vg] + values[, , ve])
  out <- rbind(out, data.frame(base, quantity = "heritability",
    parameter_name = "heritability", annotation_name = NA_character_,
    stick_name = NA_character_, value = as.vector(h2),
    stringsAsFactors = FALSE))
  chain_pis <- lapply(seq_len(4L), function(chain) {
    x <- run$fit$chains[[chain]]$pis
    if (is.null(x) || length(x) != dim(values)[1L]) return(NULL)
    data.frame(iteration = seq_along(x), chain = chain,
      quantity = "global_nonnull_proportion",
      parameter_name = "global_nonnull_proportion",
      annotation_name = NA_character_, stick_name = NA_character_,
      value = as.numeric(x), stringsAsFactors = FALSE)
  })
  if (all(vapply(chain_pis, Negate(is.null), logical(1))))
    out <- rbind(out, do.call(rbind, chain_pis))
  if (isTRUE(run$method$annotation)) {
    alpha_q <- which(q$parameter_name == "alpha" & q$stick_index == 1L)
    alpha_q <- alpha_q[order(q$annotation_index[alpha_q])]
    if (length(alpha_q) != ncol(A))
      stop("First-stick alpha trace dimensions do not match annotations.",
        call. = FALSE)
    enriched <- A[, "enriched_binary"] == 1
    derived <- vector("list", 4L)
    for (chain in seq_len(4L)) {
      coef <- values[, chain, alpha_q, drop = FALSE]
      dim(coef) <- c(dim(values)[1L], length(alpha_q))
      inside <- outside <- numeric(nrow(coef))
      chunks <- split(seq_len(nrow(coef)), ceiling(seq_len(nrow(coef)) / 50L))
      for (idx in chunks) {
        prob <- stats::pnorm(A %*% t(coef[idx, , drop = FALSE]))
        inside[idx] <- colMeans(prob[enriched, , drop = FALSE])
        outside[idx] <- colMeans(prob[!enriched, , drop = FALSE])
      }
      derived[[chain]] <- rbind(
        data.frame(iteration = seq_len(nrow(coef)), chain = chain,
          quantity = "prior_nonnull_mean_enriched",
          parameter_name = "derived_prior_probability",
          annotation_name = "enriched_binary", stick_name = NA_character_,
          value = inside),
        data.frame(iteration = seq_len(nrow(coef)), chain = chain,
          quantity = "prior_nonnull_mean_unannotated",
          parameter_name = "derived_prior_probability",
          annotation_name = "enriched_binary", stick_name = NA_character_,
          value = outside))
    }
    out <- rbind(out, do.call(rbind, derived))
  }
  if (any(!is.finite(out$value)) ||
      anyDuplicated(out[c("chain", "iteration", "quantity")]))
    stop("Invalid Study 06 chain draws.", call. = FALSE)
  out$scenario <- run$scenario
  out$replicate <- run$replicate
  out$method <- run$method$id
  out[c("scenario", "replicate", "method", "chain", "iteration",
    "quantity", "parameter_name", "annotation_name", "stick_name", "value")]
}

.study06_chain_window <- function(draws, burnin, retained) {
  z <- draws[draws$iteration > burnin &
    draws$iteration <= burnin + retained, , drop = FALSE]
  counts <- table(z$quantity, z$chain)
  if (!nrow(counts) || ncol(counts) != 4L || any(counts != retained))
    stop("Candidate window lacks equal retained draws in four chains.",
      call. = FALSE)
  z
}
