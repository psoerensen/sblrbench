#' Construct a benchmark method specification
#' @param id Machine-readable identifier.
#' @param label Human-readable label.
#' @param capabilities Unique character capabilities.
#' @param fit,extract Functions implementing fit and extraction.
#' @param predict Optional prediction function.
#' @param metadata Compact metadata list.
#' @export
new_sblrbench_method <- function(id,label,capabilities=character(),fit,extract,predict=NULL,metadata=list()) {
  out<-list(id=id,label=label,capabilities=unique(capabilities),fit=fit,extract=extract,predict=predict,metadata=metadata); class(out)<-c("sblrbench_method","list"); validate_sblrbench_method(out); out
}

#' Validate a benchmark method specification
#' @param x Object to validate.
#' @export
validate_sblrbench_method <- function(x) {
  if(!inherits(x,"sblrbench_method")) stop("x must be an sblrbench_method.",call.=FALSE)
  .assert_scalar_string(x$id,"id"); if(!grepl("^[a-z][a-z0-9_.-]*$",x$id)) stop("id must be machine-readable and start with a lower-case letter.",call.=FALSE)
  .assert_scalar_string(x$label,"label"); if(!is.character(x$capabilities)||anyNA(x$capabilities)||anyDuplicated(x$capabilities)) stop("capabilities must be a unique character vector.",call.=FALSE)
  if(!is.function(x$fit)) stop("fit must be a function.",call.=FALSE); if(!is.function(x$extract)) stop("extract must be a function.",call.=FALSE); if(!is.null(x$predict)&&!is.function(x$predict)) stop("predict must be NULL or a function.",call.=FALSE); if(!is.list(x$metadata)) stop("metadata must be a list.",call.=FALSE); invisible(x)
}

#' Run a benchmark method
#' @param method A method specification.
#' @param ... Arguments passed to its fit function.
#' @export
run_sblrbench_method <- function(method,...) { validate_sblrbench_method(method); run<-method$fit(...); result<-method$extract(run); if(!is.null(method$predict)) result<-method$predict(result,...); result }
