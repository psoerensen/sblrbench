.study06_method_specs <- function(config) {
  map <- list(
    st_bed_bayesr = list(interface = "stblr_bed", native = "bayesr",
      representation = "BED", annotation = FALSE),
    st_bed_bayesrc = list(interface = "stblr_bed", native = "bayesrc",
      representation = "BED", annotation = TRUE),
    st_csr_sbayesr = list(interface = "stblr_csr", native = "sbayesr",
      representation = "CSR", annotation = FALSE),
    st_csr_sbayesrc = list(interface = "stblr_csr_annot", native = "sbayesrc",
      representation = "CSR", annotation = TRUE))
  if (!identical(config$methods, names(map)))
    stop("Study 06 method grid differs from the frozen four-method contract.", call. = FALSE)
  lapply(seq_along(map), function(i) c(list(id = names(map)[i],
    label = unname(config$method_labels[names(map)[i]]), method_index = i), map[[i]]))
}

.study06_fit_seeds <- function(scenario, replicate, method, config) {
  seed <- config$seeds$fit_base + match(scenario, config$scenarios) * 100000L +
    as.integer(replicate) * 1000L + match(method, config$methods) * 10L
  chain <- as.integer(seed + seq_len(4L) * config$seeds$chain_stride)
  c(fit_seed = as.integer(seed), stats::setNames(chain, paste0("chain_", 1:4)))
}

.study06_baseline_recommendations <- function(config) {
  path <- config$baseline_recommendation_source
  if (!file.exists(path)) stop("Frozen Study 04 recommendations are absent.", call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  ids <- c("st_bed_bayesr", "st_csr_sbayesr")
  x <- x[x$method %in% ids, , drop = FALSE]
  x <- x[match(ids, x$method), , drop = FALSE]
  if (nrow(x) != 2L || anyNA(x$method) ||
      any(x$recommendation_status != "available") ||
      any(x$recommended_nchains != 4L) || any(x$recommended_ncores != 4L) ||
      any(x$recommended_nthin != 1L) ||
      any(x$recommended_nit_argument != x$recommended_post_burnin_draws))
    stop("Study 04 baseline recommendation contract failed.", call. = FALSE)
  data.frame(method = x$method, recommendation_status = x$recommendation_status,
    nchains = x$recommended_nchains, ncores = x$recommended_ncores,
    nburn = x$recommended_nburn, nit = x$recommended_nit_argument,
    nthin = x$recommended_nthin, limiting_estimand = x$limiting_estimand,
    limiting_diagnostic = x$limiting_diagnostic,
    source = normalizePath(path, winslash = "/"), stringsAsFactors = FALSE)
}

.study06_annotation_recommendations <- function(config) {
  path <- file.path(config$convergence_capsule, "method_recommendations.csv")
  if (!file.exists(path)) stop("Study 06 annotation recommendations are absent.", call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  ids <- c("st_bed_bayesrc", "st_csr_sbayesrc")
  x <- x[match(ids, x$method), , drop = FALSE]
  if (nrow(x) != 2L || anyNA(x$method) ||
      any(x$recommendation_status != "available") ||
      any(x$nchains != 4L) || any(x$ncores != 4L) || any(x$nthin != 1L))
    stop("Study 06 annotation recommendation contract failed.", call. = FALSE)
  x
}

.study06_mcmc_for <- function(method, config, phase = c("benchmark", "convergence")) {
  phase <- match.arg(phase)
  if (phase == "convergence") return(config$convergence_pilot[
    c("nit", "nburn", "nthin", "nchains", "ncores")])
  table <- rbind(.study06_baseline_recommendations(config),
    .study06_annotation_recommendations(config)[names(.study06_baseline_recommendations(config))])
  z <- table[table$method == method, , drop = FALSE]
  if (nrow(z) != 1L) stop("No unique MCMC recommendation for ", method, call. = FALSE)
  as.list(z[1L, c("nit", "nburn", "nthin", "nchains", "ncores")])
}

.study06_fit <- function(method, simulation, stats, Glist, split, A,
                         alpha_bundle, config,
                         phase = c("benchmark", "convergence")) {
  phase <- match.arg(phase)
  mcmc <- .study06_mcmc_for(method$id, config, phase)
  seeds <- .study06_fit_seeds(simulation$scenario, simulation$replicate,
    method$id, config)
  controls <- c(mcmc, list(seed = seeds[["fit_seed"]],
    chain_seeds = unname(seeds[paste0("chain_", 1:4)]),
    keep_chains = TRUE, verbose = FALSE, h2 = config$simulation$h2,
    mixture_var = config$mixture_var))
  marginal <- alpha_bundle$marginal_component_probability
  if (!isTRUE(method$annotation)) {
    controls$pi <- marginal
    controls$convergence <- "core"
    controls$convergence_control <- list(warn = FALSE, keep_traces = TRUE)
  } else {
    controls$add_intercept <- FALSE
    controls$standardize_annotations <- FALSE
    controls$center_binary_annotations <- FALSE
    controls$alpha_init <- alpha_bundle$uninformative_annotations
    controls$sigmaSqAlpha_init <- config$alpha_prior$sigmaSqAlpha_init
    controls$intercept_flat <- config$alpha_prior$intercept_flat
    controls$sigmaSqAlpha_a <- config$alpha_prior$sigmaSqAlpha_a
    controls$sigmaSqAlpha_b <- config$alpha_prior$sigmaSqAlpha_b
    controls$pi_floor <- config$alpha_prior$pi_floor
    controls$alpha_update_every <- config$alpha_prior$alpha_update_every
    controls$updateAlpha <- TRUE
    controls$convergence <- "extended"
    controls$convergence_control <- list(warn = FALSE,
      extended_groups = c("annotations", "probability"), keep_traces = TRUE,
      max_trace_gb = 2, allow_large_traces = FALSE)
  }
  y <- simulation$phenotype[split$train_ids, , drop = FALSE]
  if (method$representation == "BED") {
    input <- list(y = y, Glist = Glist, rows = split$train_rows)
    if (isTRUE(method$annotation)) {
      input$annotation <- as.data.frame(A)
      controls$annot_alpha_init <- controls$alpha_init
      controls$annot_sigma_sq_alpha_init <- controls$sigmaSqAlpha_init
      controls$annot_alpha_update_every <- controls$alpha_update_every
      controls$alpha_init <- NULL
      controls$sigmaSqAlpha_init <- NULL
      controls$alpha_update_every <- NULL
    }
  } else if (isTRUE(method$annotation)) {
    input <- list(stats = stats, Glist = Glist, annotations = A,
      annotation_model = "annotation_probit_stick")
  } else {
    input <- list(stats = stats, Glist = Glist)
  }
  native <- sblrbench::new_sblr_native_method(method$id, method$label,
    method$interface, method$native,
    capabilities = c("scalar_trait", "multichain", "posterior_effects",
      "component_probabilities"))
  started <- Sys.time(); elapsed <- proc.time()[["elapsed"]]
  tryCatch({
    result <- sblrbench::run_sblrbench_method(native,
      fit_inputs = input, controls = controls)
    fit <- result$native_fit
    if (is.null(fit$convergence_traces$values) ||
        dim(fit$convergence_traces$values)[2L] != 4L ||
        is.null(fit$chains) || length(fit$chains) != 4L)
      stop("Four identifiable chains and convergence traces were not retained.",
        call. = FALSE)
    if (isTRUE(method$annotation)) {
      q <- fit$convergence_traces$quantities
      if (!all(c("alpha", "sigmaSqAlpha") %in% q$parameter_name))
        stop("Chain-level annotation coefficient draws are unavailable.",
          call. = FALSE)
    }
    list(status = "ok", reason = "", method = method,
      scenario = simulation$scenario, replicate = simulation$replicate,
      result = result,
      fit = fit, controls = controls, seeds = seeds, started_at = started,
      finished_at = Sys.time(), runtime = proc.time()[["elapsed"]] - elapsed)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e),
    scenario = simulation$scenario, replicate = simulation$replicate,
    method = method, result = NULL, fit = NULL, controls = controls,
    seeds = seeds, started_at = started, finished_at = Sys.time(),
    runtime = proc.time()[["elapsed"]] - elapsed))
}
