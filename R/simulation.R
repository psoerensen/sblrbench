#' Convert an sblr simulation to the benchmark contract
#' @param x Current output from [sblr::mtsim()].
#' @param source Simulation source; currently `"mtsim"`.
#' @param study,architecture Scenario labels.
#' @param replicate Non-negative replicate index.
#' @param seed Simulation seed.
#' @param keep_genotypes Retain the genotype matrix.
#' @export
as_sblrbench_simulation <- function(x, source="mtsim", study="unspecified", architecture="unspecified", replicate=0L, seed=NULL, keep_genotypes=TRUE) {
  if (!identical(source,"mtsim") || !is.list(x)) stop("Only a list returned by sblr::mtsim() is currently supported.", call.=FALSE)
  req <- c("y","W","B","G","E","causal","rsids","ids")
  if (length(setdiff(req,names(x)))) stop("mtsim output is missing: ",paste(setdiff(req,names(x)),collapse=", "),call.=FALSE)
  marker_ids <- .assert_ids(as.character(x$rsids),"marker_ids"); sample_ids <- .assert_ids(as.character(x$ids),"sample_ids")
  trait_names <- colnames(x$B); .assert_ids(trait_names,"trait_names")
  parameters <- x[intersect(c("h2_target","h2_observed","rg_target","rg_observed","rb_shared_observed","rb_all_observed","re_target","re_observed","Sigma_e"),names(x))]
  extras <- x[setdiff(names(x),c(req,names(parameters),"shared_idx","specific_idx","causal_rsids"))]
  extras$shared_idx <- x$shared_idx %||% NULL; extras$specific_idx <- x$specific_idx %||% NULL; extras$causal_rsids <- x$causal_rsids %||% NULL
  out <- list(schema_version=1L,
    data=list(marker_ids=marker_ids,sample_ids=sample_ids,trait_names=trait_names,train_ids=NULL,test_ids=NULL,reference_ids=NULL,genotypes=if(keep_genotypes)x$W else NULL),
    truth=list(effects=x$B,genetic_values=x$G,phenotypes=x$y,residuals=x$E,causal=x$causal,parameters=parameters),
    scenario=list(study=study,architecture=architecture,replicate=as.integer(replicate)),
    provenance=list(seed=seed,simulator="sblr::mtsim",transformations=character()),extras=extras)
  class(out) <- c("sblrbench_simulation","list"); validate_sblrbench_simulation(out); out
}

#' Validate a benchmark simulation
#' @param x Object to validate.
#' @export
validate_sblrbench_simulation <- function(x) {
  if (!inherits(x,"sblrbench_simulation") || !is.list(x)) stop("x must be an sblrbench_simulation.",call.=FALSE)
  if (!is.integer(x$schema_version)||length(x$schema_version)!=1L||is.na(x$schema_version)) stop("schema_version must be a scalar integer.",call.=FALSE)
  d <- x$data; m <- .assert_ids(d$marker_ids,"marker_ids"); s <- .assert_ids(d$sample_ids,"sample_ids"); t <- .assert_ids(d$trait_names,"trait_names")
  .as_named_matrix(x$truth$effects,m,t,"truth$effects"); .as_named_matrix(x$truth$genetic_values,s,t,"truth$genetic_values"); .as_named_matrix(x$truth$phenotypes,s,t,"truth$phenotypes"); .as_named_matrix(x$truth$residuals,s,t,"truth$residuals")
  if (!is.null(d$genotypes)) .as_named_matrix(d$genotypes,s,m,"data$genotypes")
  for(nm in c("train_ids","test_ids","reference_ids")) if(!is.null(d[[nm]])){.assert_ids(d[[nm]],nm); if(anyNA(match(d[[nm]],s))) stop(nm," must be a subset of sample_ids.",call.=FALSE)}
  if(length(intersect(d$train_ids %||% character(),d$test_ids %||% character()))) stop("train_ids and test_ids must not overlap.",call.=FALSE)
  r <- x$scenario$replicate; if(!is.integer(r)||length(r)!=1L||is.na(r)||r<0L) stop("replicate must be a non-negative integer.",call.=FALSE)
  causal <- unique(c(x$truth$causal$shared %||% character(),unlist(x$truth$causal$specific %||% list(),use.names=FALSE),x$truth$causal$all %||% character()))
  if(anyNA(match(causal,m))) stop("All causal marker IDs must occur in marker_ids.",call.=FALSE)
  invisible(x)
}

#' Check oracle genetic values
#' @param simulation Optional benchmark simulation.
#' @param genotypes,effects,genetic_values Explicit matrices.
#' @param tolerance Maximum permitted absolute error.
#' @param stop_on_failure Stop when the numeric check fails.
#' @export
check_oracle_genetic_values <- function(simulation=NULL,genotypes=NULL,effects=NULL,genetic_values=NULL,tolerance=sqrt(.Machine$double.eps),stop_on_failure=TRUE) {
  if(!is.null(simulation)){validate_sblrbench_simulation(simulation);genotypes<-simulation$data$genotypes;effects<-simulation$truth$effects;genetic_values<-simulation$truth$genetic_values}
  if(is.null(genotypes)) stop("Genotypes are required for the oracle check.",call.=FALSE)
  if(!is.matrix(genotypes)||!is.matrix(effects)||!is.matrix(genetic_values)) stop("Genotypes, effects, and genetic values must be matrices.",call.=FALSE)
  if(any(!is.finite(genotypes))||any(!is.finite(effects))||any(!is.finite(genetic_values))) stop("Oracle inputs must be finite.",call.=FALSE)
  marker_ids <- .assert_ids(rownames(effects),"effect marker IDs"); traits <- .assert_ids(colnames(effects),"effect trait names"); samples <- .assert_ids(rownames(genetic_values),"genetic-value sample IDs")
  genotypes <- align_markers(genotypes,marker_ids,margin=2L); genotypes <- align_samples(genotypes,samples); genetic_values <- align_traits(genetic_values,traits)
  product <- genotypes %*% effects; if(any(!is.finite(product))) stop("Oracle matrix product contains non-finite values.",call.=FALSE)
  err <- max(abs(product-genetic_values)); out <- list(ok=err<=tolerance,max_abs_error=err,tolerance=tolerance,n_samples=nrow(genotypes),n_markers=ncol(genotypes),n_traits=ncol(effects)); class(out)<-c("sblrbench_oracle_check","list")
  if(!out$ok && isTRUE(stop_on_failure)) stop("Oracle genetic-value check failed; maximum absolute error = ",format(err)," (tolerance = ",format(tolerance),").",call.=FALSE)
  out
}

#' @export
print.sblrbench_oracle_check <- function(x,...) { cat("Oracle genetic-value check: ",if(x$ok)"PASS" else "FAIL"," (max abs error ",format(x$max_abs_error),", tolerance ",format(x$tolerance),")\n",sep=""); invisible(x) }
