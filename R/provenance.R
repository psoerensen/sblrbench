#' Detect a Git commit without changing repository state
#' @param path Repository path.
#' @export
sblrbench_git_commit <- function(path=".") {
  benchmark_git_sha(path, warn = TRUE)
}

#' Construct a compact provenance manifest
#' @param study,scenario Study and scenario labels.
#' @param replicate Non-negative replicate index.
#' @param simulation_seed,chain_seeds Simulation and fitting seeds.
#' @param sblr_commit,sblrbench_commit Source commits.
#' @param genotype_source Compact genotype-source description.
#' @param running_controls Deprecated alias controls.
#' @param fitting_controls Compact fitting controls.
#' @param runtime Compact runtime details.
#' @param convergence_status Compact convergence status.
#' @param package_versions Named compact package versions.
#' @export
new_sblrbench_manifest <- function(study,scenario,replicate,simulation_seed=NULL,chain_seeds=NULL,sblr_commit=NA_character_,sblrbench_commit=sblrbench_git_commit("."),genotype_source,running_controls=list(),fitting_controls=list(),runtime=list(),convergence_status=NA_character_,package_versions=list()) {
  if(length(running_controls)) warning("running_controls is deprecated; use fitting_controls.",call.=FALSE)
  out<-list(schema_version=1L,study=study,scenario=scenario,replicate=as.integer(replicate),simulation_seed=simulation_seed,chain_seeds=chain_seeds,sblr_commit=sblr_commit,sblrbench_commit=sblrbench_commit,genotype_source=genotype_source,fitting_controls=c(fitting_controls,running_controls),runtime=runtime,convergence_status=convergence_status,R_version=as.character(getRversion()),platform=R.version$platform,package_versions=package_versions,created_at=format(Sys.time(),tz="UTC",usetz=TRUE));class(out)<-c("sblrbench_manifest","list");validate_sblrbench_manifest(out);out
}

.compact_value <- function(x) { if(is.null(x))return(TRUE);if(is.matrix(x)||is.data.frame(x)||is.array(x)||inherits(x,"sblrbench_result")||inherits(x,"sblrbench_simulation"))return(FALSE);if(is.list(x))return(all(vapply(x,.compact_value,logical(1))));is.atomic(x)&&length(x)<=1000L }

#' Validate a benchmark provenance manifest
#' @param x Object to validate.
#' @export
validate_sblrbench_manifest <- function(x) {if(!inherits(x,"sblrbench_manifest")||!is.list(x))stop("x must be an sblrbench_manifest.",call.=FALSE);req<-c("schema_version","study","scenario","replicate","simulation_seed","chain_seeds","sblr_commit","sblrbench_commit","genotype_source","fitting_controls","runtime","convergence_status","R_version","platform","package_versions","created_at");if(length(setdiff(req,names(x))))stop("Manifest fields are missing.",call.=FALSE);if(!is.integer(x$schema_version)||length(x$schema_version)!=1L)stop("schema_version must be a scalar integer.",call.=FALSE);if(!.compact_value(unclass(x)))stop("Manifest accepts only compact atomic values and lists; matrices and fitted objects are forbidden.",call.=FALSE);invisible(x)}

#' Write a manifest as pretty JSON
#' @param x Manifest.
#' @param path Output path.
#' @export
write_sblrbench_manifest <- function(x,path) {validate_sblrbench_manifest(x);jsonlite::write_json(unclass(x),path,pretty=TRUE,auto_unbox=TRUE,null="null",digits=NA);invisible(path)}

#' Read a JSON manifest
#' @param path Input path.
#' @export
read_sblrbench_manifest <- function(path) {x<-jsonlite::read_json(path,simplifyVector=FALSE);x$schema_version<-as.integer(x$schema_version);x$replicate<-as.integer(x$replicate);class(x)<-c("sblrbench_manifest","list");validate_sblrbench_manifest(x);x}
