#' Construct a standard benchmark result
#' @param method_id Method identifier.
#' @param effects,pip,genetic_covariance Common estimates or `NULL`.
#' @param genetic_value Prediction matrix or `NULL`.
#' @param convergence,native_diagnostics,warnings Diagnostics.
#' @param elapsed_seconds,memory_estimate Computation summaries.
#' @param provenance Compact provenance.
#' @param native_fit Unmodified native fit.
#' @param keep_native_fit Retain the native fit.
#' @export
new_sblrbench_result <- function(method_id,effects=NULL,pip=NULL,genetic_covariance=NULL,genetic_value=NULL,convergence=NULL,native_diagnostics=NULL,warnings=character(),elapsed_seconds=NA_real_,memory_estimate=NULL,provenance=list(),native_fit=NULL,keep_native_fit=TRUE) {
  out<-list(schema_version=1L,method_id=method_id,estimates=list(effects=effects,pip=pip,genetic_covariance=genetic_covariance),predictions=list(genetic_value=genetic_value),diagnostics=list(convergence=convergence,native=native_diagnostics,warnings=warnings),computation=list(elapsed_seconds=as.numeric(elapsed_seconds),memory_estimate=memory_estimate),provenance=provenance,native_fit=if(keep_native_fit)native_fit else NULL); class(out)<-c("sblrbench_result","list"); validate_sblrbench_result(out); out
}

#' Validate a standard benchmark result
#' @param x Result to validate.
#' @param simulation Optional simulation defining expected identifiers.
#' @export
validate_sblrbench_result <- function(x,simulation=NULL) {
  if(!inherits(x,"sblrbench_result")) stop("x must be an sblrbench_result.",call.=FALSE); .assert_scalar_string(x$method_id,"method_id")
  mats<-list(effects=x$estimates$effects,pip=x$estimates$pip,genetic_value=x$predictions$genetic_value)
  for(nm in names(mats)) if(!is.null(mats[[nm]])){z<-mats[[nm]];if(!is.matrix(z)||!is.numeric(z)||is.null(rownames(z))||is.null(colnames(z))||any(!is.finite(z))) stop(nm," must be a finite numeric matrix with row and column names.",call.=FALSE);.assert_ids(rownames(z),paste(nm,"row IDs"));.assert_ids(colnames(z),paste(nm,"trait names"))}
  if(!is.null(x$estimates$pip)&&any(x$estimates$pip<0|x$estimates$pip>1)) stop("pip values must lie in [0, 1].",call.=FALSE)
  if(!is.null(x$estimates$genetic_covariance)){g<-x$estimates$genetic_covariance;if(!is.matrix(g)||nrow(g)!=ncol(g)||is.null(rownames(g))||!identical(rownames(g),colnames(g))||any(!is.finite(g))) stop("genetic_covariance must be a finite named square matrix.",call.=FALSE)}
  if(!is.null(simulation)){validate_sblrbench_simulation(simulation);m<-simulation$data$marker_ids;s<-simulation$data$sample_ids;t<-simulation$data$trait_names;if(!is.null(mats$effects)) .as_named_matrix(mats$effects,m,t,"effects");if(!is.null(mats$pip)) .as_named_matrix(mats$pip,m,t,"pip");if(!is.null(mats$genetic_value)) .as_named_matrix(mats$genetic_value,s,t,"genetic_value")}
  invisible(x)
}
