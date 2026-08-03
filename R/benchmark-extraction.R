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
