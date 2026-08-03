# Plain Study 02 method translation. There is intentionally no global registry
# or method class beyond the existing small sblrbench method contract.

read_benchmark_recommendations <- function(spec) {
  path <- benchmark_spec_path(spec,
    spec$controls$benchmark$recommendation_source)
  if (!file.exists(path))
    stop("Frozen Study 04 method recommendations are unavailable: ", path,
      call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  expected <- data.frame(
    method = c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc",
      "st_csr_sbayesr"),
    recommended_nburn = c(250L, 250L, 250L, 250L),
    recommended_post_burnin_draws = c(250L, 2000L, 250L, 2000L),
    recommended_nit_argument = c(250L, 2000L, 250L, 2000L),
    recommended_nthin = rep(1L, 4L), recommended_nchains = rep(4L, 4L),
    recommended_ncores = rep(4L, 4L), stringsAsFactors = FALSE)
  columns <- names(expected)
  if (!all(columns %in% names(x)) || anyDuplicated(x$method) ||
      !setequal(x$method, expected$method))
    stop("Frozen Study 04 recommendation grid is invalid.", call. = FALSE)
  observed <- x[match(expected$method, x$method), columns, drop = FALSE]
  rownames(observed) <- NULL
  if (!isTRUE(all.equal(observed, expected, check.attributes = FALSE)))
    stop("Frozen Study 04 recommendations do not match the validated Study 02 settings.",
      call. = FALSE)
  x[match(expected$method, x$method), , drop = FALSE]
}

resolve_benchmark_methods <- function(spec) {
  validate_benchmark_spec(spec)
  lapply(seq_along(spec$methods), function(i) {
    method <- spec$methods[[i]]
    method$id <- names(spec$methods)[[i]]
    method$method_index <- as.integer(i)
    method
  })
}

benchmark_method_controls <- function(spec, method_id, profile, fit_seed,
                                       chain_seeds) {
  resolve_benchmark_profile(spec, profile)
  if (!method_id %in% names(spec$methods))
    stop("Unknown benchmark method: ", method_id, call. = FALSE)
  method <- spec$methods[[method_id]]
  if (identical(profile, "benchmark")) {
    recommendation <- read_benchmark_recommendations(spec)
    row <- recommendation[recommendation$method == method_id, , drop = FALSE]
    controls <- list(nit = as.integer(row$recommended_nit_argument),
      nburn = as.integer(row$recommended_nburn),
      nthin = as.integer(row$recommended_nthin),
      nchains = as.integer(row$recommended_nchains),
      ncores = as.integer(row$recommended_ncores), convergence = "core",
      keep_chains = TRUE,
      convergence_control = list(warn = FALSE, keep_traces = TRUE))
  } else {
    controls <- spec$controls$workshop[c("nit", "nburn", "nthin",
      "nchains", "ncores", "convergence", "keep_chains",
      "convergence_control")]
  }
  controls$seed <- as.integer(fit_seed)
  controls$chain_seeds <- as.integer(chain_seeds)
  if (length(controls$chain_seeds) != controls$nchains ||
      anyDuplicated(controls$chain_seeds))
    stop("Resolved chain seeds do not match the configured chain count.",
      call. = FALSE)
  controls$verbose <- FALSE
  controls$h2 <- spec$controls$priors$h2
  if (identical(method$prior_class, "BayesR")) {
    active <- spec$controls$priors$bayesr_active_probability
    controls$pi <- c(1 - active, rep(active / 3, 3L))
    controls$mixture_var <- spec$controls$priors$bayesr_mixture_var
  } else {
    controls$pi_init <- spec$controls$priors$bayesc_inclusion_probability
  }
  controls
}

fit_prediction_method <- function(method, controls, simulation, stats, glist,
                                  split) {
  capabilities <- c("posterior_effects", "pip", "scalar_trait",
    if (identical(method$representation, "BED")) "individual_level" else
      "summary_statistics")
  method_spec <- new_sblr_native_method(method$id, method$label,
    method$interface, method$native_method, capabilities = capabilities,
    metadata = list(task = "single_trait_prediction"))
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  inputs <- if (identical(method$representation, "BED"))
    list(y = y, Glist = glist, rows = split$train_rows) else
    list(stats = stats, Glist = glist)
  result <- run_sblrbench_method(method_spec, fit_inputs = inputs,
    controls = controls)
  validate_sblrbench_result(result, simulation)
  result
}

predict_prediction_result <- function(result, simulation, test_simulation,
                                      scaled_test) {
  effects <- align_traits(align_markers(extract_marker_effects(result),
    simulation$data$marker_ids), simulation$data$trait_names)
  prediction <- scaled_test %*% effects
  add_sblrbench_predictions(result, prediction, test_simulation)
}

fit_parameter_estimation_method <- function(method, controls, simulation,
                                             stats, glist) {
  capabilities <- c("posterior_effects","pip","scalar_trait",
    if(identical(method$representation,"BED")) "individual_level" else
      "summary_statistics")
  method_spec <- new_sblr_native_method(method$id,method$label,
    method$interface,method$native_method,capabilities=capabilities,
    metadata=list(task="parameter_estimation"))
  inputs <- if(identical(method$representation,"BED"))
    list(y=simulation$truth$phenotypes,Glist=glist,
      rows=seq_len(nrow(simulation$truth$phenotypes))) else
    list(stats=stats,Glist=glist)
  result <- run_sblrbench_method(method_spec,fit_inputs=inputs,
    controls=controls)
  validate_sblrbench_result(result,simulation)
  result
}
