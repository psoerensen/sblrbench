# Extraction distinguishes native formatted fits from benchmark result objects
# and never infers unavailable traces or probabilities.

.benchmark_native_fit <- function(fit) {
  if (inherits(fit, "sblrbench_result")) return(fit$native_fit)
  if (is.list(fit) && !is.null(fit$native_fit)) return(fit$native_fit)
  fit
}

extract_marker_effects <- function(fit) {
  value <- if (inherits(fit, "sblrbench_result")) fit$estimates$effects else
    .benchmark_native_fit(fit)$bm
  if (is.null(value)) stop("Posterior mean marker effects are unavailable.",
    call. = FALSE)
  attr(value, "benchmark_quantity") <- "posterior_mean"
  value
}

extract_marker_probabilities <- function(fit) {
  native <- .benchmark_native_fit(fit)
  posterior <- if (inherits(fit, "sblrbench_result")) fit$estimates$pip else
    native$dm
  list(posterior_inclusion = posterior,
    posterior_component = native$comp_prob,
    prior_component = native$pi,
    prior_marker_component = native$pim)
}

extract_variance_components <- function(fit) {
  native <- .benchmark_native_fit(fit)
  list(posterior_means = native[intersect(c("vbs", "vgs", "ves", "vle",
    "vld"), names(native))],
    true_traces = if (is.list(native$chains)) native$chains else NULL,
    final_states = native[intersect(c("vb", "vg", "ve"), names(native))])
}

extract_fit_metadata <- function(fit) {
  native <- .benchmark_native_fit(fit)
  data.frame(
    method = if (inherits(fit, "sblrbench_result")) fit$method_id else NA_character_,
    marker_count = if (!is.null(native$bm)) nrow(native$bm) else NA_integer_,
    trait_count = if (!is.null(native$bm)) ncol(native$bm) else NA_integer_,
    has_posterior_probabilities = !is.null(native$dm),
    has_component_probabilities = !is.null(native$comp_prob),
    has_true_chain_traces = is.list(native$chains),
    stringsAsFactors = FALSE)
}

extract_runtime <- function(fit) {
  value <- if (inherits(fit, "sblrbench_result"))
    fit$computation$elapsed_seconds else fit$elapsed_seconds
  if (is.null(value) || !length(value)) NA_real_ else as.numeric(value[[1L]])
}

extract_chain_information <- function(fit) {
  native <- .benchmark_native_fit(fit)
  diagnostics <- if (inherits(fit, "sblrbench_result"))
    fit$diagnostics$convergence else native$convergence
  list(true_traces = if (is.list(native$chains)) native$chains else NULL,
    compact_summaries = diagnostics,
    final_states = native[intersect(c("b", "d", "vb", "vg", "ve", "pi"),
      names(native))])
}

#' Extract identifiable scalar convergence traces
#'
#' Only the native iteration-by-chain-by-quantity trace bundle is accepted.
#' Posterior means, final states, pooled draws, and compact summaries are not
#' substitutes for true chain traces.
#'
#' @param fit An sblr fit or `sblrbench_result`.
#' @param coordinate A one-row convergence coordinate.
#' @param registry The convergence quantity registry.
#' @param expected_chains Required number of identifiable chains.
#' @return A tidy data frame with iteration, chain, quantity, and value.
#' @export
extract_convergence_traces <- function(fit, coordinate, registry,
                                       expected_chains = 4L) {
  native <- .benchmark_native_fit(fit)
  bundle <- native$convergence_traces
  if (is.null(bundle$values) || length(dim(bundle$values)) != 3L)
    stop("True native convergence traces are unavailable; posterior means, final states, and compact summaries cannot be substituted.",
      call. = FALSE)
  if (dim(bundle$values)[2L] != expected_chains)
    stop("True convergence traces do not contain the required number of identifiable chains.",
      call. = FALSE)
  source <- as.character(bundle$quantities$group)
  required <- registry[registry$required, , drop = FALSE]
  primitive <- required[required$source %in% source, , drop = FALSE]
  missing_primitive <- required[!required$source %in% source &
    required$source != "vgs/(vgs+ves)", , drop = FALSE]
  if (nrow(missing_primitive))
    stop("True convergence traces omit required scalar quantities: ",
      paste(missing_primitive$quantity, collapse = ", "), ".",
      call. = FALSE)
  metadata <- as.list(coordinate[intersect(c("stage", "scenario",
    "replicate", "method"), names(coordinate))])
  index <- match(primitive$source, source)
  out <- benchmark_trace_array_long(bundle$values[, , index, drop = FALSE],
    primitive$quantity, metadata = metadata,
    expected_chains = expected_chains)
  if ("heritability" %in% required$quantity) {
    derived_index <- match(c("vgs", "ves"), source)
    if (anyNA(derived_index))
      stop("Heritability traces require true vgs and ves chain traces.",
        call. = FALSE)
    vg <- bundle$values[, , derived_index[[1L]]]
    ve <- bundle$values[, , derived_index[[2L]]]
    h2 <- benchmark_trace_array_long(array(vg / (vg + ve),
      dim = c(dim(vg), 1L)), "heritability", metadata = metadata,
      expected_chains = expected_chains)
    out <- rbind(out, h2)
  }
  if (any(!is.finite(out$value)))
    stop("True convergence traces contain non-finite values.", call. = FALSE)
  bounds <- required[match(out$quantity, required$quantity), , drop = FALSE]
  if (any(out$value < bounds$lower_bound | out$value > bounds$upper_bound))
    stop("True convergence traces violate registered quantity bounds.",
      call. = FALSE)
  out
}

benchmark_trace_vector <- function(x, name) {
  if (is.null(x)) stop("Missing posterior component: ", name, call.=FALSE)
  if (is.vector(x)) x <- matrix(x,ncol=1L)
  if (!is.matrix(x) || ncol(x)!=1L)
    stop("Expected one-column posterior trace: ",name,call.=FALSE)
  as.numeric(x[,1L])
}

extract_parameter_draws <- function(fit, method, estimands, marker_count,
                                    expected_chains=NULL) {
  native <- .benchmark_native_fit(fit)
  if (is.null(native$input$nit) || is.null(native$input$nburn))
    stop("Fit does not record nit and nburn.",call.=FALSE)
  if (!is.null(expected_chains) && expected_chains > 1L) {
    bundle <- native$convergence_traces
    if (is.null(bundle$values) || length(dim(bundle$values))!=3L)
      stop("Native convergence trace bundle is unavailable.",call.=FALSE)
    values <- bundle$values
    if (dim(values)[2L]!=expected_chains || dim(values)[1L]!=native$input$nit)
      stop("Multichain dimensions disagree with retained draws.",call.=FALSE)
    idx <- match(c("vbs","vgs","ves"),as.character(bundle$quantities$group))
    if (anyNA(idx)) stop("Required scalar traces are absent.",call.=FALSE)
    base <- expand.grid(iteration=seq_len(dim(values)[1L]),
      chain=seq_len(expected_chains))
    base$vbs <- as.vector(values[,,idx[1L]])
    base$vgs <- as.vector(values[,,idx[2L]])
    base$ves <- as.vector(values[,,idx[3L]])
    pi_values <- NULL
    if (is.list(native$chains) && length(native$chains)==expected_chains) {
      pis <- lapply(native$chains,function(ch) {
        z <- if(is.null(ch$pi_trace)) NULL else as.numeric(ch$pi_trace)
        if(is.null(z)) return(NULL)
        if(length(z)==native$input$nit) z else if(length(z)==native$input$nit+
          native$input$nburn) tail(z,native$input$nit) else NULL
      })
      if(all(vapply(pis,length,integer(1))==native$input$nit))
        pi_values <- unlist(pis,use.names=FALSE)
    }
  } else {
    components <- lapply(c("vbs","vgs","ves"),function(n)
      benchmark_trace_vector(native[[n]],n))
    names(components) <- c("vbs","vgs","ves")
    if(!is.null(native$pi_trace)) components$pi_trace <-
      benchmark_trace_vector(native$pi_trace,"pi_trace")
    lengths <- vapply(components,length,integer(1))
    if(length(unique(lengths))!=1L)
      stop("Posterior traces are not jointly aligned.",call.=FALSE)
    nit <- as.integer(native$input$nit); nburn <- as.integer(native$input$nburn)
    nthin <- as.integer(if(is.null(native$input$nthin)) 1L else native$input$nthin)
    if(lengths[1L]==nburn+nit) keep <- seq.int(nburn+1L,nburn+nit,by=nthin)
    else if(lengths[1L]==nit && nburn<nit) keep <- seq.int(nburn+1L,nit,by=nthin)
    else stop("Posterior trace length is inconsistent with fit metadata.",call.=FALSE)
    z <- lapply(components,`[`,keep)
    base <- data.frame(iteration=seq_along(z$vbs),chain=1L,
      vbs=z$vbs,vgs=z$vgs,ves=z$ves)
    pi_values <- z$pi_trace
  }
  if(any(!is.finite(as.matrix(base[c("vbs","vgs","ves")]))) ||
      any(as.matrix(base[c("vbs","vgs","ves")])<0))
    stop("Scalar posterior draws are invalid.",call.=FALSE)
  derived <- list(effect_variance=base$vbs,genetic_variance=base$vgs,
    residual_variance=base$ves,heritability=base$vgs/(base$vgs+base$ves))
  if(!is.null(pi_values)) derived <- c(list(causal_proportion=pi_values,
    total_marker_effect_variance=base$vbs*pi_values*marker_count),derived)
  do.call(rbind,lapply(names(derived),function(id) data.frame(method=method,
    estimand_id=id,chain=base$chain,iteration=base$iteration,
    value=derived[[id]],source_component=estimands$posterior_source[
      match(id,estimands$estimand_id)],status="ok",reason="",
    stringsAsFactors=FALSE)))
}

summarise_parameter_draws <- function(draws) {
  do.call(rbind,lapply(split(draws,draws$estimand_id),function(x) data.frame(
    method=x$method[1L],estimand_id=x$estimand_id[1L],
    posterior_mean=mean(x$value),posterior_sd=stats::sd(x$value),
    posterior_median=stats::median(x$value),
    lower_95=unname(stats::quantile(x$value,.025)),
    upper_95=unname(stats::quantile(x$value,.975)),
    n_posterior_draws=nrow(x),chain_count=length(unique(x$chain)),
    draws_per_chain=unique(as.integer(table(x$chain)))[1L],
    mcse_mean=NA_real_,status="ok",reason="",stringsAsFactors=FALSE)))
}

complete_parameter_summary <- function(summary, estimands, method,
                                       chain_count) {
  missing <- setdiff(estimands$estimand_id,summary$estimand_id)
  if(!length(missing)) return(summary)
  absent <- data.frame(method=method,estimand_id=missing,
    posterior_mean=NA_real_,posterior_sd=NA_real_,posterior_median=NA_real_,
    lower_95=NA_real_,upper_95=NA_real_,n_posterior_draws=0L,
    chain_count=as.integer(chain_count),draws_per_chain=0L,mcse_mean=NA_real_,
    status="unavailable",
    reason="aligned draw-level posterior component is not retained",
    stringsAsFactors=FALSE)
  rbind(summary,absent)
}

#' Extract true annotation coefficient traces
#'
#' This accepts only the identifiable iteration-by-chain-by-quantity trace
#' bundle produced by the annotation-aware samplers. Posterior means and final
#' states are never substituted.
#'
#' @param fit An sblr fit or `sblrbench_result`.
#' @param expected_chains Required number of identifiable chains.
#' @return A tidy data frame. Unavailable traces are represented by one row
#'   with `status = "unavailable"` and an actionable reason.
#' @export
extract_annotation_coefficient_traces <- function(fit,
                                                   expected_chains = 4L) {
  native <- .benchmark_native_fit(fit)
  bundle <- native$convergence_traces
  unavailable <- function(reason) data.frame(iteration = NA_integer_,
    chain = NA_integer_, parameter = NA_character_, annotation = NA_character_,
    stick = NA_character_, value = NA_real_, status = "unavailable",
    reason = reason, stringsAsFactors = FALSE)
  if (is.null(bundle$values) || length(dim(bundle$values)) != 3L)
    return(unavailable(paste("True annotation coefficient traces are absent;",
      "posterior means and final states are not substitutes.")))
  if (dim(bundle$values)[2L] != as.integer(expected_chains))
    return(unavailable("Annotation traces lack the required identifiable chains."))
  q <- bundle$quantities
  required_columns <- c("parameter_name", "annotation_name", "stick_name")
  if (!is.data.frame(q) || !all(required_columns %in% names(q)) ||
      nrow(q) != dim(bundle$values)[3L])
    return(unavailable("Annotation trace descriptors are incomplete."))
  index <- which(q$parameter_name %in% c("alpha", "sigmaSqAlpha"))
  if (!length(index) || !all(c("alpha", "sigmaSqAlpha") %in%
      q$parameter_name[index]))
    return(unavailable("True alpha and sigmaSqAlpha traces are both required."))
  base <- expand.grid(iteration = seq_len(dim(bundle$values)[1L]),
    chain = seq_len(dim(bundle$values)[2L]))
  rows <- lapply(index, function(j) data.frame(base,
    parameter = as.character(q$parameter_name[j]),
    annotation = ifelse(is.na(q$annotation_name[j]), "",
      as.character(q$annotation_name[j])),
    stick = ifelse(is.na(q$stick_name[j]), "",
      as.character(q$stick_name[j])),
    value = as.vector(bundle$values[, , j]), status = "ok", reason = "",
    stringsAsFactors = FALSE))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (any(!is.finite(out$value)) ||
      anyDuplicated(out[c("iteration", "chain", "parameter", "annotation",
        "stick")]))
    return(unavailable("Annotation coefficient traces are non-finite or duplicated."))
  out
}

extract_annotation_fit_components <- function(fit) {
  native <- .benchmark_native_fit(fit)
  list(
    posterior_means = list(alpha = native$alpha,
      sigmaSqAlpha = native$sigmaSqAlpha),
    final_states = list(alpha = native$alpha_final,
      sigmaSqAlpha = native$sigmaSqAlpha_final,
      marker_prior = native$annotation_prior),
    true_draws = extract_annotation_coefficient_traces(fit,
      expected_chains = if (is.list(native$chains)) length(native$chains) else
        4L),
    marker_posterior_components = native$component_probabilities,
    fit_input = native$input)
}

#' Summarise a comparable draw-wise annotation marker prior
#'
#' For each retained chain draw, this reconstructs the complete alpha matrix
#' and evaluates the package's probit-stick transformation. Marker-level output
#' is the posterior mean of those draw-wise probabilities. Summary output
#' retains draw identity for convergence and uncertainty calculations.
#'
#' @param traces Output from `extract_annotation_coefficient_traces()`.
#' @param annotations Ordered marker-by-annotation matrix.
#' @param mixture_var Ordered mixture variances.
#' @param true_marker_prior Optional true marker component probabilities.
#' @param marker_truth Optional data frame containing `true_nonnull`.
#' @param retain_marker_summary Whether to retain posterior marker-component
#'   means. Qualification gates may set this to `FALSE` because their
#'   registered draw-wise non-null summaries depend only on the first stick.
#' @return A list with status, marker summaries, draw summaries, and reason.
#' @export
summarise_drawwise_annotation_prior <- function(traces, annotations,
                                                mixture_var,
                                                true_marker_prior = NULL,
                                                marker_truth = NULL,
                                                retain_marker_summary = TRUE) {
  if (!is.data.frame(traces) || !nrow(traces) ||
      any(traces$status != "ok"))
    return(list(status = "unavailable", marker = NULL, draws = NULL,
      reason = if (is.data.frame(traces) && "reason" %in% names(traces))
        paste(unique(traces$reason[nzchar(traces$reason)]), collapse = "; ") else
          "True annotation coefficient traces are unavailable."))
  alpha <- traces[traces$parameter == "alpha", , drop = FALSE]
  annotations_required <- colnames(annotations)
  sticks_required <- paste0("component_",
    seq_len(length(mixture_var) - 1L) - 1L, "_stick")
  if (!setequal(unique(alpha$annotation), annotations_required) ||
      !setequal(unique(alpha$stick), sticks_required))
    return(list(status = "unavailable", marker = NULL, draws = NULL,
      reason = "Alpha traces do not span every annotation and probit stick."))
  keys <- unique(alpha[c("chain", "iteration")])
  keys <- keys[order(keys$chain, keys$iteration), , drop = FALSE]
  key_id <- paste(keys$chain, keys$iteration, sep = "\r")
  alpha$key_index <- match(paste(alpha$chain, alpha$iteration, sep = "\r"),
    key_id)
  if (anyNA(alpha$key_index) ||
      nrow(alpha) != nrow(keys) * length(annotations_required) *
        length(sticks_required))
    return(list(status = "unavailable", marker = NULL, draws = NULL,
      reason = "At least one retained draw has an incomplete alpha matrix."))
  marker_sum <- if (isTRUE(retain_marker_summary))
    matrix(0, nrow(annotations), length(mixture_var),
      dimnames = list(rownames(annotations),
        paste0("component_", seq_len(length(mixture_var)) - 1L))) else NULL
  enriched <- annotations[, "enriched_binary"] == 1
  causal <- if (is.null(marker_truth)) rep(NA, nrow(annotations)) else
    as.logical(marker_truth$true_nonnull[match(rownames(annotations),
      marker_truth$marker_id)])
  draw_rows <- vector("list", ceiling(nrow(keys) / 64L))
  chunks <- split(seq_len(nrow(keys)),
    ceiling(seq_len(nrow(keys)) / 64L))
  for (i in seq_along(chunks)) {
    draw_index <- chunks[[i]]
    remaining <- matrix(1, nrow(annotations), length(draw_index))
    nonnull <- NULL
    continuous_prior_contrast <- rep(NA_real_, length(draw_index))
    sticks_to_transform <- if (isTRUE(retain_marker_summary))
      seq_along(sticks_required) else 1L
    for (j in sticks_to_transform) {
      z <- alpha[alpha$stick == sticks_required[j] &
        alpha$key_index %in% draw_index, , drop = FALSE]
      coefficient <- matrix(NA_real_, length(annotations_required),
        length(draw_index))
      coefficient[cbind(match(z$annotation, annotations_required),
        match(z$key_index, draw_index))] <- z$value
      if (anyNA(coefficient))
        return(list(status = "unavailable", marker = NULL, draws = NULL,
          reason = "At least one retained draw has an incomplete alpha matrix."))
      stick_probability <- stats::pnorm(annotations %*% coefficient)
      if (any(!is.finite(stick_probability)))
        stop("Draw-wise annotation prior transformation failed.",
          call. = FALSE)
      component_probability <- remaining * (1 - stick_probability)
      if (isTRUE(retain_marker_summary))
        marker_sum[, j] <- marker_sum[, j] +
          rowSums(component_probability)
      remaining <- remaining * stick_probability
      if (j == 1L) {
        nonnull <- stick_probability
        intercept_row <- match("Intercept", annotations_required)
        continuous_row <- match("continuous_signal", annotations_required)
        if (!is.na(intercept_row) && !is.na(continuous_row))
          continuous_prior_contrast <-
            stats::pnorm(coefficient[intercept_row, ] +
              coefficient[continuous_row, ]) -
            stats::pnorm(coefficient[intercept_row, ] -
              coefficient[continuous_row, ])
      }
    }
    if (isTRUE(retain_marker_summary))
      marker_sum[, length(mixture_var)] <-
        marker_sum[, length(mixture_var)] + rowSums(remaining)
    draw_rows[[i]] <- data.frame(chain = keys$chain[draw_index],
      iteration = keys$iteration[draw_index],
      expected_active = colSums(nonnull),
      mean_prior_enriched = colMeans(nonnull[enriched, , drop = FALSE]),
      mean_prior_unannotated = colMeans(nonnull[!enriched, , drop = FALSE]),
      enriched_prior_contrast =
        colMeans(nonnull[enriched, , drop = FALSE]) -
          colMeans(nonnull[!enriched, , drop = FALSE]),
      continuous_prior_contrast = continuous_prior_contrast,
      mean_prior_causal = if (all(is.na(causal))) NA_real_ else
        colMeans(nonnull[causal, , drop = FALSE]),
      mean_prior_noncausal = if (all(is.na(causal))) NA_real_ else
        colMeans(nonnull[!causal, , drop = FALSE]),
      stringsAsFactors = FALSE)
  }
  marker <- NULL
  if (isTRUE(retain_marker_summary)) {
    posterior_mean <- marker_sum / nrow(keys)
    marker <- data.frame(marker_id = rownames(annotations),
      posterior_mean_nonnull_prior = 1 - posterior_mean[, 1L],
      stringsAsFactors = FALSE)
    for (j in seq_len(ncol(posterior_mean)))
      marker[[paste0("posterior_mean_prior_component_", j - 1L)]] <-
        posterior_mean[, j]
    if (!is.null(true_marker_prior)) {
      true_marker_prior <- true_marker_prior[match(marker$marker_id,
        rownames(true_marker_prior)), , drop = FALSE]
      marker$true_nonnull_prior <- 1 - true_marker_prior[, 1L]
    }
  }
  list(status = "ok", marker = marker,
    draws = do.call(rbind, draw_rows), reason = "")
}

prediction_estimate_table <- function(result, scenario, replicate, method) {
  effects <- extract_marker_effects(result)
  data.frame(study = "02_prediction", scenario = scenario,
    replicate = as.integer(replicate), method = method,
    marker = rep(rownames(effects), times = ncol(effects)),
    trait = rep(colnames(effects), each = nrow(effects)),
    estimate = as.numeric(effects), quantity = "posterior_mean_effect",
    stringsAsFactors = FALSE)
}

prediction_marker_table <- function(result, scenario, replicate, method) {
  effects <- extract_marker_effects(result)
  probabilities <- extract_marker_probabilities(result)$posterior_inclusion
  out <- data.frame(study = "02_prediction", scenario = scenario,
    replicate = as.integer(replicate), method = method,
    marker = rep(rownames(effects), times = ncol(effects)),
    trait = rep(colnames(effects), each = nrow(effects)),
    posterior_mean_effect = as.numeric(effects), stringsAsFactors = FALSE)
  if (!is.null(probabilities)) {
    probabilities <- align_traits(align_markers(probabilities,
      rownames(effects)), colnames(effects))
    out$posterior_inclusion_probability <- as.numeric(probabilities)
  }
  out
}
