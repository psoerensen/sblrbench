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

#' Evaluate selected truth-aware metrics
#' @inheritParams metric_effect_rmse
#' @param metrics Metric names.
#' @export
evaluate_metrics <- function(simulation,result,metrics=c("effect_rmse","prediction_correlation","prediction_mse","pip_brier")) {
  funs<-list(effect_rmse=metric_effect_rmse,prediction_correlation=metric_prediction_correlation,prediction_mse=metric_prediction_mse,pip_brier=metric_pip_brier);bad<-setdiff(metrics,names(funs));if(length(bad))stop("Unknown metrics: ",paste(bad,collapse=", "),call.=FALSE);do.call(rbind,lapply(metrics,function(n)funs[[n]](simulation,result)))
}
