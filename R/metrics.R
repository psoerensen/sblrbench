.metric_prepare <- function(simulation,result,kind) { validate_sblrbench_simulation(simulation);validate_sblrbench_result(result); switch(kind,effects=result$estimates$effects,pip=result$estimates$pip,prediction=result$predictions$genetic_value) }
.skip_traits <- function(sim,result,metric,reason) do.call(rbind,lapply(sim$data$trait_names,function(t).metric_row(sim,result,t,metric,status="skipped",reason=reason)))

#' All-marker posterior-mean effect RMSE
#' @param simulation A validated benchmark simulation.
#' @param result A benchmark result.
#' @export
metric_effect_rmse <- function(simulation,result) {
  z<-.metric_prepare(simulation,result,"effects");if(is.null(z))return(.skip_traits(simulation,result,"effect_rmse","estimated effects are unavailable"))
  z<-align_markers(z,simulation$data$marker_ids);z<-align_traits(z,simulation$data$trait_names);truth<-simulation$truth$effects
  do.call(rbind,lapply(simulation$data$trait_names,function(t).metric_row(simulation,result,t,"effect_rmse",sqrt(mean((z[,t]-truth[,t])^2)))))
}

#' Prediction correlation with true genetic value
#' @inheritParams metric_effect_rmse
#' @export
metric_prediction_correlation <- function(simulation,result) {
  z<-.metric_prepare(simulation,result,"prediction");if(is.null(z))return(.skip_traits(simulation,result,"prediction_correlation","genetic-value predictions are unavailable"))
  z<-align_samples(z,simulation$data$sample_ids);z<-align_traits(z,simulation$data$trait_names);truth<-simulation$truth$genetic_values
  do.call(rbind,lapply(simulation$data$trait_names,function(t){if(stats::sd(z[,t])==0||stats::sd(truth[,t])==0)return(.metric_row(simulation,result,t,"prediction_correlation",status="failed",reason="prediction or truth has zero variance"));.metric_row(simulation,result,t,"prediction_correlation",stats::cor(z[,t],truth[,t]))}))
}

#' Prediction MSE against true genetic value
#' @inheritParams metric_effect_rmse
#' @export
metric_prediction_mse <- function(simulation,result) {
  z<-.metric_prepare(simulation,result,"prediction");if(is.null(z))return(.skip_traits(simulation,result,"prediction_mse","genetic-value predictions are unavailable"));z<-align_samples(z,simulation$data$sample_ids);z<-align_traits(z,simulation$data$trait_names);truth<-simulation$truth$genetic_values
  do.call(rbind,lapply(simulation$data$trait_names,function(t).metric_row(simulation,result,t,"prediction_mse",mean((z[,t]-truth[,t])^2))))
}

.causal_matrix <- function(sim) { out<-matrix(0,nrow=length(sim$data$marker_ids),ncol=length(sim$data$trait_names),dimnames=list(sim$data$marker_ids,sim$data$trait_names)); shared<-sim$truth$causal$shared %||% character();if(length(shared))out[shared,]<-1;sp<-sim$truth$causal$specific %||% list();for(t in intersect(names(sp),colnames(out)))if(length(sp[[t]]))out[sp[[t]],t]<-1;out }

#' Trait-specific PIP Brier score
#' @inheritParams metric_effect_rmse
#' @export
metric_pip_brier <- function(simulation,result) {
  z<-.metric_prepare(simulation,result,"pip");if(is.null(z))return(.skip_traits(simulation,result,"pip_brier","PIPs are unavailable"));z<-align_markers(z,simulation$data$marker_ids);z<-align_traits(z,simulation$data$trait_names);truth<-.causal_matrix(simulation)
  do.call(rbind,lapply(simulation$data$trait_names,function(t).metric_row(simulation,result,t,"pip_brier",mean((z[,t]-truth[,t])^2))))
}

.pip_ranking <- function(simulation, result) {
  z <- .metric_prepare(simulation, result, "pip")
  if (is.null(z)) return(NULL)
  z <- align_markers(z, simulation$data$marker_ids)
  z <- align_traits(z, simulation$data$trait_names)
  if (any(!is.finite(z)) || any(z < 0 | z > 1))
    stop("PIPs must be finite and lie in [0, 1].", call. = FALSE)
  z
}

#' Average precision of causal-marker ranking
#'
#' Markers are ordered by decreasing PIP, with ties resolved by first
#' occurrence in canonical marker order. Average precision is the mean of
#' precision at the one-based ranks occupied by causal markers.
#' @inheritParams metric_effect_rmse
#' @export
metric_average_precision <- function(simulation, result) {
  z <- .pip_ranking(simulation, result)
  if (is.null(z)) return(.skip_traits(simulation, result, "average_precision", "PIPs are unavailable"))
  truth <- .causal_matrix(simulation)
  do.call(rbind, lapply(simulation$data$trait_names, function(t) {
    causal <- truth[, t] == 1
    if (!any(causal)) return(.metric_row(simulation, result, t, "average_precision", status = "skipped", reason = "trait has no causal markers"))
    ord <- order(-z[, t], seq_len(nrow(z)))
    ranked <- causal[ord]
    value <- mean(cumsum(ranked)[ranked] / which(ranked))
    .metric_row(simulation, result, t, "average_precision", value)
  }))
}

#' Causal-marker rank summaries
#'
#' Ranks are one-based after decreasing-PIP sorting; ties use first occurrence
#' in canonical marker order. Top-K recall is reported for K equal to the
#' number of causals, twice that number, and 50 (each capped at marker count).
#' @inheritParams metric_effect_rmse
#' @export
metric_causal_ranks <- function(simulation, result) {
  z <- .pip_ranking(simulation, result)
  names0 <- c("causal_rank_mean", "causal_rank_median", "causal_rank_best", "causal_rank_worst")
  if (is.null(z)) return(do.call(rbind, lapply(c(names0, "causal_top_1_recall", "causal_top_2_recall", "causal_top_50_recall"), function(n) .skip_traits(simulation, result, n, "PIPs are unavailable"))))
  truth <- .causal_matrix(simulation)
  do.call(rbind, lapply(simulation$data$trait_names, function(t) {
    causal <- truth[, t] == 1
    if (!any(causal)) return(do.call(rbind, lapply(names0, function(n) .metric_row(simulation, result, t, n, status = "skipped", reason = "trait has no causal markers"))))
    ord <- order(-z[, t], seq_len(nrow(z)))
    ranks <- match(which(causal), ord)
    nc <- sum(causal); ks <- unique(c(nc, 2L * nc, 50L))
    base <- c(mean(ranks), stats::median(ranks), min(ranks), max(ranks))
    rows <- Map(function(n, v) .metric_row(simulation, result, t, n, v), names0, base)
    rows <- c(rows, lapply(ks, function(k) .metric_row(simulation, result, t, paste0("causal_top_", k, "_recall"), mean(ranks <= min(k, nrow(z))))))
    do.call(rbind, rows)
  }))
}

#' Evaluate selected truth-aware metrics
#' @inheritParams metric_effect_rmse
#' @param metrics Metric names.
#' @export
evaluate_metrics <- function(simulation,result,metrics=c("effect_rmse","prediction_correlation","prediction_mse","pip_brier")) {
  funs<-list(effect_rmse=metric_effect_rmse,prediction_correlation=metric_prediction_correlation,prediction_mse=metric_prediction_mse,pip_brier=metric_pip_brier,average_precision=metric_average_precision,causal_ranks=metric_causal_ranks);bad<-setdiff(metrics,names(funs));if(length(bad))stop("Unknown metrics: ",paste(bad,collapse=", "),call.=FALSE);do.call(rbind,lapply(metrics,function(n)funs[[n]](simulation,result)))
}
