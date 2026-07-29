config <- source(file.path("studies", "00_contract_smoke", "config.R"), local=TRUE)$value
list(
  targets::tar_target(smoke_simulation, {
    set.seed(config$seed)
    W <- matrix(sample(0:2, config$n*config$m, replace=TRUE), config$n, config$m,
                dimnames=list(paste0("sample",seq_len(config$n)),paste0("marker",seq_len(config$m))))
    raw <- sblr::mtsim(W=W,standardize_W=TRUE,nt=config$nt,n_shared=config$n_shared,n_specific=config$n_specific,seed=config$seed)
    sblrbench::as_sblrbench_simulation(raw,study="00_contract_smoke",architecture="contract",replicate=0L,seed=config$seed)
  }),
  targets::tar_target(smoke_oracle,sblrbench::check_oracle_genetic_values(smoke_simulation)),
  targets::tar_target(smoke_result, {
    causal <- matrix(0,nrow=length(smoke_simulation$data$marker_ids),ncol=length(smoke_simulation$data$trait_names),dimnames=list(smoke_simulation$data$marker_ids,smoke_simulation$data$trait_names))
    causal[smoke_simulation$truth$causal$shared,] <- 1
    for(t in names(smoke_simulation$truth$causal$specific)) causal[smoke_simulation$truth$causal$specific[[t]],t] <- 1
    method <- sblrbench::new_sblrbench_method("oracle","Contract oracle",c("posterior_effects","pip","prediction"),fit=function(...)list(effects=smoke_simulation$truth$effects,pip=causal,prediction=smoke_simulation$truth$genetic_values),extract=function(x)sblrbench::new_sblrbench_result("oracle",effects=x$effects,pip=x$pip,genetic_value=x$prediction))
    sblrbench::run_sblrbench_method(method)
  }),
  targets::tar_target(smoke_metrics,sblrbench::evaluate_metrics(smoke_simulation,smoke_result)),
  targets::tar_target(smoke_outputs, {
    dir.create(config$output_dir,recursive=TRUE,showWarnings=FALSE)
    utils::write.csv(smoke_metrics,file.path(config$output_dir,"metrics.csv"),row.names=FALSE)
    manifest <- sblrbench::new_sblrbench_manifest("00_contract_smoke","contract",0L,config$seed,sblr_commit="fd76b18828bc77756948aa3138de07ae4dc75513",genotype_source="base-R generated matrix",runtime=list(oracle_ok=smoke_oracle$ok),convergence_status="not_applicable",package_versions=list(sblr=as.character(utils::packageVersion("sblr"))))
    sblrbench::write_sblrbench_manifest(manifest,file.path(config$output_dir,"manifest.json"))
    c(file.path(config$output_dir,"metrics.csv"),file.path(config$output_dir,"manifest.json"))
  },format="file")
)
