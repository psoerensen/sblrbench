.sblr_interfaces <- c("stblr_csr","stblr_csr_annot","stblr_block_eigen","stblr_bed","mtblr_csr","mtblr_block_eigen","mtblr_bed")

#' Create a native public-API sblr method
#' @inheritParams new_sblrbench_method
#' @param interface One of the seven canonical exported fitting interfaces.
#' @param method Optional model argument.
#' @export
new_sblr_native_method <- function(id,label,interface,method=NULL,capabilities,metadata=list()) {
  interface<-match.arg(interface,.sblr_interfaces)
  fit_fun<-function(fit_inputs=list(),controls=list()) {
    if(!is.list(fit_inputs)||is.null(names(fit_inputs))||!is.list(controls)||is.null(names(controls))) stop("fit_inputs and controls must be named lists.",call.=FALSE)
    dup<-intersect(names(fit_inputs),names(controls)); if(length(dup)) stop("Duplicated fit/control arguments: ",paste(dup,collapse=", "),call.=FALSE)
    args<-c(fit_inputs,controls); if(!is.null(method)){if("method"%in%names(args)) stop("method conflicts with the factory-supplied argument.",call.=FALSE);args$method<-method}
    exports<-getNamespaceExports("sblr"); if(!interface%in%exports) stop("Installed sblr does not export ",interface,"().",call.=FALSE)
    fun<-getExportedValue("sblr",interface); warns<-msgs<-character(); start<-proc.time()[["elapsed"]]
    native <- withCallingHandlers(
      do.call(fun, args),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    list(native_fit=native,elapsed_seconds=unname(proc.time()[["elapsed"]]-start),warnings=warns,messages=msgs)
  }
  extract_fun<-function(run) extract_sblr_result(run,method_id=id)
  new_sblrbench_method(id,label,capabilities,fit_fun,extract_fun,metadata=c(metadata,list(interface=interface,method=method)))
}

#' Extract compatible common fields from an sblr run
#' @param run Native run record or native fit list.
#' @param method_id Benchmark method identifier.
#' @param keep_native_fit Retain unmodified fit.
#' @export
extract_sblr_result <- function(run,method_id,keep_native_fit=TRUE) {
  if(!is.list(run)) stop("run must be a list.",call.=FALSE); fit<-run$native_fit %||% run; elapsed<-run$elapsed_seconds %||% NA_real_; warns<-run$warnings %||% character()
  new_sblrbench_result(method_id=method_id,effects=fit$bm %||% NULL,pip=fit$dm %||% NULL,genetic_covariance=fit$cov_g_mean %||% NULL,convergence=fit$convergence %||% NULL,native_diagnostics=fit$diagnostics %||% NULL,warnings=warns,elapsed_seconds=elapsed,memory_estimate=fit$memory_estimate %||% NULL,native_fit=fit,keep_native_fit=keep_native_fit)
}
