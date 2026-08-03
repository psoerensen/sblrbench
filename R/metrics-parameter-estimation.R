# Study 03 parameter-recovery metrics. Formulas preserve the frozen capsule.

parameter_model_specification <- function(scenario, prior_class) {
  matched <- (scenario=="sparse_homogeneous" & prior_class=="BayesC") |
    (scenario=="sparse_mixture" & prior_class=="BayesR")
  ifelse(matched,"matched","misspecified")
}

parameter_recovery_metrics <- function(summary, truth, method, estimands,
                                       tolerance=1e-10) {
  x <- merge(summary,truth,by="estimand_id",all.x=TRUE,sort=FALSE)
  x$status <- x$status.x
  x$reason <- ifelse(nzchar(x$reason.x),x$reason.x,x$reason.y)
  x$study <- "03_parameter_estimation"; x$method <- method$id
  x$method_label <- method$id; x$prior_class <- method$prior_class
  x$model_specification <- parameter_model_specification(x$architecture,
    method$prior_class)
  x$estimand_label <- estimands$label[match(x$estimand_id,
    estimands$estimand_id)]
  x$bias <- x$posterior_mean-x$truth
  x$absolute_error <- abs(x$bias); x$squared_error <- x$bias^2
  allowed <- estimands$relative_error_allowed[match(x$estimand_id,
    estimands$estimand_id)] & abs(x$truth)>tolerance
  x$relative_error <- ifelse(allowed,x$bias/x$truth,NA_real_)
  x$covered_95 <- x$lower_95<=x$truth & x$truth<=x$upper_95
  x$interval_width_95 <- x$upper_95-x$lower_95
  x[,c("study","architecture","replicate","method","method_label",
    "prior_class","model_specification","estimand_id","estimand_label",
    "truth","truth_type","posterior_mean","posterior_sd",
    "posterior_median","lower_95","upper_95","bias","absolute_error",
    "squared_error","relative_error","covered_95","interval_width_95",
    "n_posterior_draws","chain_count","draws_per_chain","status","reason")]
}

parameter_pair_definitions <- function() data.frame(
  comparison_id=c("csr_vs_bed_bayesc","csr_vs_bed_bayesr",
    "bayesr_vs_bayesc_bed","sbayesr_vs_sbayesc_csr"),
  focal_method=c("st_csr_sbayesc","st_csr_sbayesr","st_bed_bayesr",
    "st_csr_sbayesr"),
  comparison_method=c("st_bed_bayesc","st_bed_bayesr","st_bed_bayesc",
    "st_csr_sbayesc"),
  comparison_type=c("BED-versus-CSR","BED-versus-CSR","prior-class",
    "prior-class"),stringsAsFactors=FALSE)

parameter_paired_differences <- function(x) {
  keys <- c("architecture","replicate","estimand_id")
  defs <- parameter_pair_definitions(); out <- list()
  if(anyDuplicated(x[c(keys,"method")]))
    stop("Duplicate recovery rows prevent pairing.",call.=FALSE)
  for(i in seq_len(nrow(defs))) {
    a <- x[x$method==defs$focal_method[i],]
    b <- x[x$method==defs$comparison_method[i],]
    z <- merge(a,b,by=keys,suffixes=c("_focal","_comparison"),all=FALSE)
    if(!nrow(z)) stop("Incomplete parameter comparison rows.",call.=FALSE)
    complete <- is.finite(z$posterior_mean_focal) &
      is.finite(z$posterior_mean_comparison)
    out[[i]] <- data.frame(z[keys],defs[i,],
      estimate_difference=z$posterior_mean_focal-z$posterior_mean_comparison,
      absolute_error_difference=z$absolute_error_focal-
        z$absolute_error_comparison,
      interval_width_difference=z$interval_width_95_focal-
        z$interval_width_95_comparison,
      coverage_agreement=z$covered_95_focal==z$covered_95_comparison,
      complete_pair=complete,
      reason=ifelse(complete,"","one or both method-specific estimands unavailable"),
      stringsAsFactors=FALSE)
  }
  do.call(rbind,out)
}

parameter_recovery_summary <- function(x) do.call(rbind,lapply(split(x,
  interaction(x$architecture,x$method,x$model_specification,x$estimand_id,
    drop=TRUE)),function(z) {
  ok <- z$status=="ok" & is.finite(z$posterior_mean); v <- z[ok,,drop=FALSE]
  data.frame(architecture=z$architecture[1L],method=z$method[1L],
    model_specification=z$model_specification[1L],estimand_id=z$estimand_id[1L],
    replicate_count=length(unique(z$replicate)),successful_replicates=nrow(v),
    mean_truth=if(nrow(v)) mean(v$truth) else NA_real_,
    mean_posterior_mean=if(nrow(v)) mean(v$posterior_mean) else NA_real_,
    mean_bias=if(nrow(v)) mean(v$bias) else NA_real_,
    rmse=if(nrow(v)) sqrt(mean(v$squared_error)) else NA_real_,
    mae=if(nrow(v)) mean(v$absolute_error) else NA_real_,
    observed_coverage_count=if(nrow(v)) sum(v$covered_95) else 0L,
    observed_coverage_proportion=if(nrow(v)) mean(v$covered_95) else NA_real_,
    coverage_note="coarse observed coverage; each miss changes five-replicate coverage by 20 percentage points",
    empirical_coverage_95=if(nrow(v)>1L) mean(v$covered_95) else NA_real_,
    coverage_observation=if(nrow(v)==1L) v$covered_95[1L] else NA,
    mean_interval_width_95=if(nrow(v)) mean(v$interval_width_95) else NA_real_,
    minimum=if(nrow(v)) min(v$posterior_mean) else NA_real_,
    maximum=if(nrow(v)) max(v$posterior_mean) else NA_real_,
    failure_count=sum(!ok),stringsAsFactors=FALSE)
}))

parameter_paired_summary <- function(x) do.call(rbind,lapply(split(x,
  interaction(x$architecture,x$comparison_id,x$estimand_id,drop=TRUE)),
  function(z) {
    v <- z[z$complete_pair,,drop=FALSE]
    data.frame(architecture=z$architecture[1L],comparison_id=z$comparison_id[1L],
      comparison_type=z$comparison_type[1L],estimand_id=z$estimand_id[1L],
      replicate_count=length(unique(z$replicate)),successful_replicates=nrow(v),
      mean_estimate_difference=if(nrow(v)) mean(v$estimate_difference) else NA_real_,
      sd_estimate_difference=if(nrow(v)>1L) stats::sd(v$estimate_difference) else NA_real_,
      median_estimate_difference=if(nrow(v)) stats::median(v$estimate_difference) else NA_real_,
      minimum_estimate_difference=if(nrow(v)) min(v$estimate_difference) else NA_real_,
      maximum_estimate_difference=if(nrow(v)) max(v$estimate_difference) else NA_real_,
      mean_absolute_error_difference=if(nrow(v)) mean(v$absolute_error_difference) else NA_real_,
      mean_interval_width_difference=if(nrow(v)) mean(v$interval_width_difference) else NA_real_,
      coverage_agreement_observation=if(nrow(v)==1L) v$coverage_agreement else NA,
      stringsAsFactors=FALSE)
  }))
