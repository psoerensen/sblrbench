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

read_annotation_baseline_recommendations <- function(spec) {
  path <- benchmark_spec_path(spec,
    spec$controls$baseline$recommendation_source)
  if (!file.exists(path))
    stop("Frozen Study 04 baseline recommendations are unavailable: ", path,
      call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  ids <- c("st_bed_bayesr", "st_csr_sbayesr")
  x <- x[match(ids, x$method), , drop = FALSE]
  required <- c("method", "recommendation_status", "recommended_nburn",
    "recommended_nit_argument", "recommended_nthin", "recommended_nchains",
    "recommended_ncores")
  if (!all(required %in% names(x)) || anyNA(x$method) ||
      any(x$recommendation_status != "available") ||
      any(x$recommended_nburn != 250L) ||
      !identical(as.integer(x$recommended_nit_argument), c(2000L, 2000L)) ||
      any(x$recommended_nthin != 1L) || any(x$recommended_nchains != 4L) ||
      any(x$recommended_ncores != 4L))
    stop("Study 04 BayesR/SBayesR recommendations differ from the audited Study 06 baseline controls.",
      call. = FALSE)
  x
}

annotation_method_controls <- function(spec, coordinate,
                                       profile = "benchmark",
                                       mode = c("qualification", "final"),
                                       qualification_decision = NULL) {
  mode <- match.arg(mode)
  method_id <- as.character(coordinate$method)
  if (!method_id %in% names(spec$methods))
    stop("Unknown Study 06 method: ", method_id, call. = FALSE)
  method <- spec$methods[[method_id]]
  if (identical(mode, "qualification") && !isTRUE(method$annotation_aware))
    stop("Qualification mode accepts annotation-aware methods only.",
      call. = FALSE)
  if (identical(profile, "workshop") && identical(mode, "final")) {
    controls <- spec$controls$workshop
  } else if (identical(mode, "qualification")) {
    controls <- list(nit = as.integer(spec$qualification$maximum_history),
      nburn = as.integer(spec$qualification$nburn),
      nthin = as.integer(spec$qualification$nthin),
      nchains = as.integer(spec$qualification$nchains),
      ncores = as.integer(spec$qualification$ncores),
      convergence = "extended", keep_chains = TRUE,
      convergence_control = list(warn = FALSE,
        extended_groups = c("annotations", "probability"),
        keep_traces = TRUE, max_trace_gb = 2,
        allow_large_traces = FALSE))
  } else if (isTRUE(method$annotation_aware)) {
    if (is.null(qualification_decision) ||
        !identical(qualification_decision$overall_decision, "passed"))
      stop("Final annotation controls require a passing qualification decision.",
        call. = FALSE)
    entries <- qualification_decision$entries
    hit <- vapply(entries, function(x)
      identical(x$scenario, as.character(coordinate$scenario)) &&
        identical(x$method, method_id), logical(1))
    if (sum(hit) != 1L)
      stop("Qualification decision has no unique compatible entry for ",
        coordinate$scenario, " / ", method_id, ".", call. = FALSE)
    selected <- entries[[which(hit)]]
    controls <- list(nit = as.integer(selected$selected_retained),
      nburn = as.integer(selected$selected_burnin), nthin = 1L,
      nchains = 4L, ncores = 4L, convergence = "extended",
      keep_chains = TRUE,
      convergence_control = list(warn = FALSE,
        extended_groups = c("annotations", "probability"),
        keep_traces = TRUE, max_trace_gb = 2,
        allow_large_traces = FALSE))
  } else {
    recommendation <- read_annotation_baseline_recommendations(spec)
    row <- recommendation[recommendation$method == method_id, , drop = FALSE]
    controls <- list(nit = as.integer(row$recommended_nit_argument),
      nburn = as.integer(row$recommended_nburn),
      nthin = as.integer(row$recommended_nthin),
      nchains = as.integer(row$recommended_nchains),
      ncores = as.integer(row$recommended_ncores), convergence = "core",
      keep_chains = TRUE,
      convergence_control = list(warn = FALSE, keep_traces = TRUE))
  }
  controls$seed <- as.integer(coordinate$fit_seed)
  controls$chain_seeds <- as.integer(coordinate$chain_seeds[[1L]])
  if (length(controls$chain_seeds) != controls$nchains ||
      anyDuplicated(controls$chain_seeds))
    stop("Study 06 chain seeds do not match the configured chain count.",
      call. = FALSE)
  controls$verbose <- FALSE
  controls$h2 <- spec$controls$priors$h2
  controls$mixture_var <- spec$controls$priors$bayesr_mixture_var
  if (isTRUE(method$annotation_aware)) {
    controls$add_intercept <- FALSE
    controls$standardize_annotations <- FALSE
    controls$center_binary_annotations <- FALSE
    controls$sigmaSqAlpha_init <- spec$controls$priors$sigmaSqAlpha_init
    controls$intercept_flat <- spec$controls$priors$intercept_flat
    controls$sigmaSqAlpha_a <- spec$controls$priors$sigmaSqAlpha_a
    controls$sigmaSqAlpha_b <- spec$controls$priors$sigmaSqAlpha_b
    controls$pi_floor <- spec$controls$priors$pi_floor
    controls$alpha_update_every <-
      as.integer(spec$controls$priors$alpha_update_every)
    controls$updateAlpha <- spec$controls$priors$updateAlpha
    controls$updateB <- spec$controls$priors$updateB
    controls$updateE <- spec$controls$priors$updateE
  }
  controls
}

fit_annotation_method <- function(method, controls, simulation, stats, glist,
                                  split, annotations, annotation_truth) {
  if (!identical(rownames(annotations), simulation$data$marker_ids))
    stop("Annotation rows are not aligned to fitted marker order.",
      call. = FALSE)
  if (isTRUE(method$annotation_aware)) {
    controls$alpha_init <- annotation_truth$uninformative_annotations
  } else {
    controls$pi <- annotation_truth$marginal_component_probability
  }
  y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
  if (identical(method$representation, "BED")) {
    inputs <- list(y = y, Glist = glist, rows = split$train_rows)
    if (isTRUE(method$annotation_aware)) {
      inputs$annotation <- as.data.frame(annotations)
      controls$annot_alpha_init <- controls$alpha_init
      controls$annot_sigma_sq_alpha_init <- controls$sigmaSqAlpha_init
      controls$annot_alpha_update_every <- controls$alpha_update_every
      controls$alpha_init <- NULL
      controls$sigmaSqAlpha_init <- NULL
      controls$alpha_update_every <- NULL
    }
  } else if (isTRUE(method$annotation_aware)) {
    inputs <- list(stats = stats, Glist = glist, annotations = annotations,
      annotation_model = "annotation_probit_stick")
  } else {
    inputs <- list(stats = stats, Glist = glist)
  }
  capabilities <- c("posterior_effects", "component_probabilities",
    "scalar_trait", "multichain",
    if (isTRUE(method$annotation_aware)) "annotation_coefficients" else NULL,
    if (identical(method$representation, "BED")) "individual_level" else
      "summary_statistics")
  method_spec <- new_sblr_native_method(method$id, method$label,
    method$interface, method$native_method, capabilities = capabilities,
    metadata = list(task = "annotation_models",
      annotation_aware = isTRUE(method$annotation_aware)))
  result <- run_sblrbench_method(method_spec, fit_inputs = inputs,
    controls = controls)
  validate_sblrbench_result(result, simulation)
  if (isTRUE(method$annotation_aware)) {
    trace <- extract_annotation_coefficient_traces(result,
      expected_chains = controls$nchains)
    if (!identical(trace$status[[1L]], "ok"))
      stop(trace$reason[[1L]], call. = FALSE)
  }
  result
}

benchmark_method_controls <- function(spec, method_id, profile, fit_seed,
                                       chain_seeds) {
  resolve_benchmark_profile(spec, profile)
  if (!method_id %in% names(spec$methods))
    stop("Unknown benchmark method: ", method_id, call. = FALSE)
  method <- spec$methods[[method_id]]
  if (identical(spec$task, "finemapping")) {
    controls <- spec$controls[[profile]][c("nit", "nburn", "nthin",
      "nchains", "ncores", "convergence", "keep_chains",
      "convergence_control")]
    controls$seed <- as.integer(fit_seed)
    controls$chain_seeds <- as.integer(chain_seeds)
    controls$verbose <- FALSE
    return(controls)
  }
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

fit_finemapping_method <- function(method, controls, simulation, stats, glist) {
  capabilities <- c("posterior_effects", "pip", "scalar_trait",
    if (identical(method$representation, "BED")) "individual_level" else
      "summary_statistics")
  method_spec <- new_sblr_native_method(method$id,method$label,
    method$interface,method$native_method,capabilities=capabilities,
    metadata=list(task="finemapping",development_settings=TRUE))
  inputs <- if(identical(method$representation,"BED"))
    list(y=simulation$truth$phenotypes,Glist=glist) else
    list(stats=stats,Glist=glist)
  result <- run_sblrbench_method(method_spec,fit_inputs=inputs,controls=controls)
  validate_sblrbench_result(result,simulation)
  sblr::check_stblr_consistency(result$native_fit,verbose=FALSE)
  result
}

convergence_method_controls <- function(spec, parameter_spec, coordinate,
                                        recommendations = NULL) {
  method_id <- as.character(coordinate$method)
  if (!method_id %in% names(parameter_spec$methods))
    stop("Unknown matched convergence method: ", method_id, call. = FALSE)
  if (identical(as.character(coordinate$stage), "selection")) {
    controls <- spec$controls$selection
  } else {
    if (is.null(recommendations))
      stop("Validation controls require selected Study 04 recommendations.",
        call. = FALSE)
    row <- recommendations[recommendations$method == method_id, , drop = FALSE]
    if (nrow(row) != 1L || row$recommendation_status != "available")
      stop("No unique available convergence recommendation for method `",
        method_id, "`.", call. = FALSE)
    controls <- spec$controls$validation[c("nthin", "nchains", "ncores",
      "convergence", "keep_chains", "convergence_control")]
    controls$nit <- as.integer(row$recommended_nit_argument)
    controls$nburn <- as.integer(row$recommended_nburn)
  }
  controls$seed <- as.integer(coordinate$fit_seed)
  controls$chain_seeds <- as.integer(coordinate$chain_seeds[[1L]])
  controls$verbose <- FALSE
  controls$h2 <- parameter_spec$controls$priors$h2
  method <- parameter_spec$methods[[method_id]]
  if (identical(method$prior_class, "BayesR")) {
    active <- parameter_spec$controls$priors$bayesr_active_probability
    controls$pi <- c(1 - active, rep(active / 3, 3L))
    controls$mixture_var <- parameter_spec$controls$priors$bayesr_mixture_var
  } else {
    controls$pi_init <-
      parameter_spec$controls$priors$bayesc_inclusion_probability
  }
  controls
}
