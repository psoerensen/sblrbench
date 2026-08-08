study07_joint_controls <- function(spec, annotation = FALSE,
                                   selection = FALSE) {
  groups <- if (annotation) c("annotations", "probability") else "probability"
  list(
    nit = spec$joint$nit,
    nburn = spec$joint$nburn,
    nthin = spec$joint$nthin,
    nchains = spec$joint$nchains,
    ncores = spec$joint$ncores,
    seed = if (selection) spec$joint$selection_seed else if (annotation)
      spec$joint$seed else spec$joint$baseline_seed,
    chain_seeds = if (selection) spec$joint$selection_chain_seeds else
      if (annotation) spec$joint$chain_seeds else
        spec$joint$baseline_chain_seeds,
    convergence = "extended",
    keep_chains = TRUE,
    convergence_control = list(warn = FALSE, extended_groups = groups,
      keep_traces = TRUE, aggregate_component_states = TRUE,
      max_trace_gb = 2, allow_large_traces = FALSE))
}

study07_public_block_call <- function(method, spec, inputs) {
  annotation <- identical(method, "SBayesRC")
  controls <- study07_joint_controls(spec, annotation = annotation)
  args <- c(list(
    stats = inputs$stats,
    Glist = inputs$Glist,
    block_start = inputs$block_start,
    method = if (annotation) "sbayesrc" else "sbayesr",
    representation = spec$operator$representation,
    eigen_policy = spec$operator$eigen_policy,
    eigen_prop = spec$operator$eigen_prop,
    residual_policy = spec$operator$residual_policy,
    block_ve_mode = spec$operator$block_ve_mode,
    block_ve_keep_history = TRUE,
    h2 = spec$model$h2,
    updateB = spec$operator$updateB,
    updateE = spec$operator$updateE,
    verbose = FALSE), controls)
  if (!annotation) {
    args$mixture_var <- spec$model$gamma
    args$pi <- spec$model$sbayesr_pi
  }
  if (annotation) {
    args$gamma <- spec$model$gamma
    args$annotation <- inputs$truth$annotations
    args$add_intercept <- FALSE
    args$standardize_annotations <- FALSE
    args$center_binary_annotations <- FALSE
    args$alpha_init <- inputs$alpha_init
    args$sigmaSqAlpha_init <- spec$model$sigmaSqAlpha_fixed_em
    args$sigmaSqAlpha_a <- 2
    args$sigmaSqAlpha_b <- 2
    args$pi_floor <- 1e-12
    args$alpha_update_every <- 1L
    args$updateAlpha <- TRUE
    args$annotation_intercept_prior <- inputs$intercept_spec
  }
  do.call(sblr::stblr_block_eigen, args)
}

study07_joint_selection_call <- function(spec, inputs) {
  controls <- study07_joint_controls(spec, annotation = TRUE,
    selection = TRUE)
  intercept <- rbind(inputs$intercept_native,
    update_sigmaSqAlpha = rep(1, 3L),
    allocation_updates_per_cycle = rep(1, 3L),
    annotation_updates_per_cycle = rep(1, 3L))
  m <- spec$source$marker_count
  stats <- inputs$stats
  sblr:::.st_sbayesrc_selection_csr(
    stats$wy, stats$ww, as.numeric(stats$yy),
    list(rep(0, m)), list(rep(0L, m)), TRUE,
    lapply(stats$wy, as.numeric), TRUE, inputs$csr_prefix,
    matrix(spec$model$initial_B, 1L, 1L),
    matrix(spec$model$initial_E, 1L, 1L),
    list(spec$model$ssb_prior), list(spec$model$sse_prior),
    inputs$truth$annotations, spec$model$gamma, inputs$alpha_init,
    as.integer(spec$model$selection$joint_delta_init),
    spec$model$selection$joint_pi_A_init,
    spec$model$selection$joint_tau2_init,
    spec$model$selection$joint_pi_prior[["a"]],
    spec$model$selection$joint_pi_prior[["b"]],
    spec$model$selection$joint_tau_prior[["a"]],
    spec$model$selection$joint_tau_prior[["b"]],
    integer(), TRUE, TRUE, TRUE, intercept, 1e-12, 4, 4,
    spec$operator$updateB, spec$operator$updateE, 0.9,
    spec$source$training_count,
    controls$nit, controls$nburn, controls$nthin, controls$ncores,
    controls$seed, controls$nchains, controls$chain_seeds, TRUE, FALSE)
}

study07_continuous_em_call <- function(start, spec, inputs) {
  alpha <- study07_alpha_starts(inputs)[[start]]
  sblr:::.stblr_mcem_sbayesrc_block_eigen(
    stats = inputs$stats, Glist = inputs$Glist,
    annotation = inputs$truth$annotations,
    block_start = inputs$block_start,
    B = matrix(spec$model$initial_B, 1L, 1L),
    E = matrix(spec$model$initial_E, 1L, 1L),
    ssb_prior = spec$model$ssb_prior,
    sse_prior = spec$model$sse_prior,
    gamma = spec$model$gamma,
    alpha_init = alpha,
    sigmaSqAlpha_init = spec$model$sigmaSqAlpha_fixed_em,
    intercept_prior_resolved = inputs$intercept_native,
    representation = spec$operator$representation,
    eigen_prop = spec$operator$eigen_prop,
    residual_policy = spec$operator$residual_policy,
    block_ve_mode = spec$operator$block_ve_mode,
    updateB = spec$operator$updateB,
    updateE = spec$operator$updateE,
    inner_sweeps = spec$em$inner_sweeps,
    inner_burn = spec$em$inner_burn,
    final_sweeps = spec$em$final_sweeps,
    final_burn = spec$em$final_burn,
    damping = spec$em$damping,
    tol_alpha = spec$em$tol_alpha,
    tol_prior = spec$em$tol_prior,
    min_outer = spec$em$min_outer,
    max_outer = spec$em$max_outer,
    ncores = 1L,
    seed = unname(spec$em$seeds[[start]]),
    return_responsibilities = TRUE,
    verbose = FALSE)
}

study07_selection_em_call <- function(start, spec, inputs) {
  alpha <- study07_alpha_starts(inputs)[[start]]
  delta <- study07_delta_starts()[[start]]
  sblr:::.stblr_mcem_sbayesrc_s_block_eigen(
    stats = inputs$stats, Glist = inputs$Glist,
    annotation = inputs$truth$annotations,
    block_start = inputs$block_start,
    B = matrix(spec$model$initial_B, 1L, 1L),
    E = matrix(spec$model$initial_E, 1L, 1L),
    ssb_prior = spec$model$ssb_prior,
    sse_prior = spec$model$sse_prior,
    gamma = spec$model$gamma,
    alpha_init = alpha,
    delta_init = delta,
    pi_a = spec$model$selection$pi_A_fixed_em,
    tau2 = spec$model$selection$tau2_fixed_em,
    intercept_prior_resolved = inputs$intercept_native,
    representation = spec$operator$representation,
    eigen_prop = spec$operator$eigen_prop,
    residual_policy = spec$operator$residual_policy,
    block_ve_mode = spec$operator$block_ve_mode,
    updateB = spec$operator$updateB,
    updateE = spec$operator$updateE,
    inner_sweeps = spec$em$inner_sweeps,
    inner_burn = spec$em$inner_burn,
    selection_sweeps = spec$em$selection_sweeps,
    selection_burn = spec$em$selection_burn,
    final_sweeps = spec$em$final_sweeps,
    final_burn = spec$em$final_burn,
    damping = spec$em$damping,
    tol_alpha = spec$em$tol_alpha,
    tol_prior = spec$em$tol_prior,
    min_outer = spec$em$min_outer,
    max_outer = spec$em$max_outer,
    ncores = 1L,
    seed = unname(spec$em$selection_seeds[[start]]),
    return_responsibilities = TRUE,
    verbose = FALSE)
}

study07_checkpoint_path <- function(spec, method, start = "primary") {
  id <- gsub("[^A-Za-z0-9]+", "_", tolower(method))
  file.path(spec$output$local_dir, "checkpoints",
    paste0(id, "--", start, ".rds"))
}

study07_fit_identity <- function(method, start, spec, validation) {
  controls <- switch(method,
    SBayesR = study07_joint_controls(spec),
    SBayesRC = study07_joint_controls(spec, TRUE),
    `SBayesRC-S` = c(study07_joint_controls(spec, TRUE, TRUE),
      list(delta_init = spec$model$selection$joint_delta_init)),
    `SBayesRC-EM` = c(spec$em, list(seed = spec$em$seeds[[start]])),
    `SBayesRC-S-EM` = c(spec$em,
      list(seed = spec$em$selection_seeds[[start]],
        delta_init = study07_delta_starts()[[start]])))
  benchmark_semantic_checkpoint_identity("study07-joint-em-fit", list(
    study = spec$study,
    study_version = spec$study_version,
    source_study = spec$source$study,
    source_scenario = spec$source$scenario,
    source_replicate = spec$source$replicate,
    method = method,
    start = start,
    design_hash = validation$design_hash,
    input_hash = validation$input_hash,
    controls = controls,
    sblr_sha = spec$packages$sblr$sha,
    installed_tree_sha256 = spec$packages$sblr$installed_tree_sha256))
}

study07_run_fit <- function(method, spec, inputs, validation,
                            start = "primary", resume = TRUE) {
  identity <- study07_fit_identity(method, start, spec, validation)
  hash <- benchmark_semantic_checkpoint_hash(identity)
  path <- study07_checkpoint_path(spec, method, start)
  if (isTRUE(resume) && file.exists(path)) {
    loaded <- benchmark_load_semantic_checkpoint(path, hash,
      validator = function(x) !is.null(x$result) && is.finite(x$elapsed_seconds))
    return(c(loaded$value, list(reused = TRUE, checkpoint = path)))
  }
  elapsed <- system.time({
    fit <- switch(method,
      SBayesR = study07_public_block_call(method, spec, inputs),
      SBayesRC = study07_public_block_call(method, spec, inputs),
      `SBayesRC-S` = study07_joint_selection_call(spec, inputs),
      `SBayesRC-EM` = study07_continuous_em_call(start, spec, inputs),
      `SBayesRC-S-EM` = study07_selection_em_call(start, spec, inputs),
      stop("Unknown Study 07 method: ", method, call. = FALSE))
  })[["elapsed"]]
  payload <- list(checkpoint_schema = "sblrbench-semantic-v2",
    identity_payload = identity, semantic_hash = hash,
    elapsed_seconds = as.numeric(elapsed), result = fit)
  benchmark_atomic_save_rds(payload, path, compress = FALSE,
    temporary_prefix = ".study07-fit-")
  c(payload, list(reused = FALSE, checkpoint = path))
}
