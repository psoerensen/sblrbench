# Study 06 v2 paired power and allocation-isolation diagnostic.
# This file deliberately does not alter the formal qualification registry.

study06_power_profile <- function() list(
  id = "v2_paired_power_isolation",
  shuffle_seed = 6201L,
  output = file.path("results", "local", "06_annotation_models",
    "v2_paired_power_isolation"),
  trace_policy = "selected_component_trace_unavailable_native_abort",
  primary_metric = "pip_auprc")

study06_power_registry <- function(spec) {
  condition <- rep(c("baseline", "learned_informative",
    "learned_shuffled", "fixed_true_alpha"), each = 2L)
  route <- rep(c("bed", "block_eigen"), 4L)
  method <- c("st_bed_bayesr", "st_block_eigen_sbayesr",
    rep(c("st_bed_bayesrc", "st_block_eigen_sbayesrc"), 3L))
  out <- data.frame(
    fit_id = paste(condition, route, sep = "--"), condition = condition,
    route = route, method = method,
    annotation_treatment = c("none", "none", "informative", "informative",
      "shuffled", "shuffled", "informative", "informative"),
    update_alpha = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    fixed_true_alpha = c(rep(FALSE, 6L), TRUE, TRUE),
    stringsAsFactors = FALSE)
  if (nrow(out) != 8L || anyDuplicated(out$fit_id) ||
      any(!out$method %in% names(spec$methods)))
    stop("The paired power diagnostic registry is invalid.", call. = FALSE)
  out
}

study06_shuffle_annotations <- function(annotations, seed) {
  if (!is.matrix(annotations) || ncol(annotations) < 2L ||
      !identical(as.numeric(annotations[, 1L]), rep(1, nrow(annotations))))
    stop("The informative annotation matrix has no exact intercept.",
      call. = FALSE)
  set.seed(as.integer(seed))
  permutation <- sample.int(nrow(annotations), replace = FALSE)
  shuffled <- annotations
  shuffled[, -1L] <- annotations[permutation, -1L, drop = FALSE]
  if (!identical(shuffled[, 1L], annotations[, 1L]) ||
      !all(vapply(seq.int(2L, ncol(annotations)), function(j)
        identical(unname(sort(shuffled[, j])),
          unname(sort(annotations[, j]))), logical(1))) ||
      max(abs(stats::cor(shuffled[, -1L, drop = FALSE]) -
        stats::cor(annotations[, -1L, drop = FALSE]))) > 1e-14)
    stop("The shuffled annotation audit failed.", call. = FALSE)
  list(annotations = shuffled, permutation = permutation,
    permutation_hash = benchmark_hash_object(permutation),
    annotation_hash = benchmark_hash_object(shuffled))
}

study06_power_controls <- function(spec, marker_ids, fit_seed, chain_seeds,
                                   annotation_aware, update_alpha) {
  groups <- if (isTRUE(annotation_aware))
    c("annotations", "probability") else "probability"
  convergence_control <- list(warn = FALSE, extended_groups = groups,
    keep_traces = TRUE, max_trace_gb = 2, allow_large_traces = FALSE)
  if (length(marker_ids)) {
    convergence_control$selected_markers <- marker_ids
    convergence_control$selected_marker_quantities <- "component"
  }
  controls <- list(
    nit = as.integer(spec$qualification$maximum_history),
    nburn = as.integer(spec$qualification$nburn),
    nthin = as.integer(spec$qualification$nthin),
    nchains = as.integer(spec$qualification$nchains),
    ncores = as.integer(spec$qualification$ncores),
    convergence = "extended", keep_chains = TRUE,
    convergence_control = convergence_control,
    seed = as.integer(fit_seed), chain_seeds = as.integer(chain_seeds),
    verbose = FALSE, h2 = spec$controls$priors$h2,
    mixture_var = spec$controls$priors$bayesr_mixture_var)
  if (isTRUE(annotation_aware)) {
    controls$add_intercept <- FALSE
    controls$standardize_annotations <- FALSE
    controls$center_binary_annotations <- FALSE
    controls$sigmaSqAlpha_init <- spec$controls$priors$sigmaSqAlpha_init
    controls$sigmaSqAlpha_a <- spec$controls$priors$sigmaSqAlpha_a
    controls$sigmaSqAlpha_b <- spec$controls$priors$sigmaSqAlpha_b
    controls$pi_floor <- spec$controls$priors$pi_floor
    controls$alpha_update_every <-
      as.integer(spec$controls$priors$alpha_update_every)
    controls$updateAlpha <- isTRUE(update_alpha)
    controls$updateB <- spec$controls$priors$updateB
    controls$updateE <- spec$controls$priors$updateE
  }
  controls$diagnostic_mode <- if (isTRUE(annotation_aware) &&
      !isTRUE(update_alpha)) "fixed_true_annotation_coefficients" else
    "standard"
  controls
}

study06_trace_marker_set <- function(marker_truth, budget_gib = 1) {
  bytes_per_marker <- 4 * 4 * 9000 + 8 * 4 * 9000 * 6 +
    8 * 4 * 9000 + 256 + 8 * 48
  maximum <- floor(budget_gib * 1024^3 / bytes_per_marker)
  causal <- which(marker_truth$true_nonnull)
  noncausal <- which(!marker_truth$true_nonnull)
  if (maximum < length(causal))
    stop("The trace budget cannot retain every causal marker.", call. = FALSE)
  extra <- maximum - length(causal)
  selected_noncausal <- if (extra) noncausal[unique(pmax(1L,
    pmin(length(noncausal), round(seq(1, length(noncausal),
      length.out = extra)))))] else integer()
  if (length(selected_noncausal) != extra)
    stop("Deterministic noncausal trace selection was not unique.",
      call. = FALSE)
  index <- sort(c(causal, selected_noncausal))
  list(index = index, marker_ids = marker_truth$marker_id[index],
    marker_count = length(index), causal_count = length(causal),
    noncausal_count = length(selected_noncausal),
    budget_gib = budget_gib, estimated_extended_gib =
      length(index) * bytes_per_marker / 1024^3,
    complete_genomewide_occupancy = length(index) == nrow(marker_truth))
}

study06_power_truth_identity <- function(data, bundle) {
  simulation <- bundle$simulation
  list(
    marker_ids = simulation$data$marker_ids,
    sample_ids = simulation$data$sample_ids,
    train_ids = data$split$train_ids, test_ids = data$split$test_ids,
    component_index = bundle$marker_truth$component_index,
    causal = bundle$marker_truth$true_nonnull,
    effects = simulation$truth$effects,
    phenotype = simulation$truth$phenotypes,
    genetic_values = simulation$truth$genetic_values,
    summary_statistics = bundle$stats,
    block_start = data$block_start,
    block_panel = data$marker_panel,
    mixture_var = bundle$spec$controls$simulation$mixture_var,
    realized_heritability = simulation$extras$realized_heritability)
}

study06_component_trace <- function(result, marker_ids) {
  native <- .benchmark_native_fit(result)
  traces <- native$convergence_traces
  q <- traces$quantities
  index <- which(q$group == "selected_component")
  if (length(index) != length(marker_ids) ||
      !identical(as.character(q$marker_id[index]), marker_ids))
    stop("Complete ordered marker-component traces are unavailable.",
      call. = FALSE)
  values <- traces$values[, , index, drop = FALSE]
  dim(values) <- c(dim(traces$values)[1L], dim(traces$values)[2L],
    length(index))
  values
}

study06_occupancy_summary <- function(component, marker_truth, fit_id) {
  nit <- dim(component)[1L]
  nchains <- dim(component)[2L]
  causal <- which(marker_truth$true_nonnull)
  rows <- vector("list", nchains)
  for (chain in seq_len(nchains)) {
    state <- component[, chain, , drop = TRUE]
    counts <- vapply(0:3, function(k) rowSums(state == k), numeric(nit))
    causal_counts <- vapply(0:3, function(k)
      rowSums(state[, causal, drop = FALSE] == k), numeric(nit))
    rows[[chain]] <- data.frame(fit_id = fit_id, iteration = seq_len(nit),
      chain = chain, traced_component_0 = counts[, 1L],
      traced_component_1 = counts[, 2L], traced_component_2 = counts[, 3L],
      traced_component_3 = counts[, 4L],
      traced_active_count = rowSums(counts[, -1L, drop = FALSE]),
      causal_active_count = rowSums(causal_counts[, -1L, drop = FALSE]),
      causal_component_1 = causal_counts[, 2L],
      causal_component_2 = causal_counts[, 3L],
      causal_component_3 = causal_counts[, 4L], stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

study06_long_trace <- function(value, quantity, fit_id) {
  base <- expand.grid(iteration = seq_len(nrow(value)),
    chain = seq_len(ncol(value)))
  data.frame(base, quantity = quantity, value = as.vector(value),
    scenario = "shared_informative_truth", replicate = 1L,
    method = fit_id, stringsAsFactors = FALSE)
}

study06_power_required_traces <- function(result, component, occupancy,
                                          annotations, marker_truth,
                                          fit_row, spec) {
  native <- .benchmark_native_fit(result)
  values <- native$convergence_traces$values
  q <- native$convergence_traces$quantities
  index <- match(c("vbs", "vgs", "ves"), q$group)
  if (anyNA(index)) stop("Core variance traces are unavailable.", call. = FALSE)
  out <- list(
    study06_long_trace(values[, , index[1L]], "effect_variance", fit_row$fit_id),
    study06_long_trace(values[, , index[2L]], "genetic_variance", fit_row$fit_id),
    study06_long_trace(values[, , index[3L]], "residual_variance", fit_row$fit_id),
    study06_long_trace(values[, , index[2L]] /
      (values[, , index[2L]] + values[, , index[3L]]), "heritability",
      fit_row$fit_id))
  probability <- which(q$group == "component_pi")
  for (j in probability) {
    component_name <- q$component_name[j]
    if (is.na(component_name) || !nzchar(component_name))
      component_name <- as.character(j)
    out[[length(out) + 1L]] <- study06_long_trace(values[, , j],
      paste0("global_component_probability:", component_name), fit_row$fit_id)
  }
  if (length(probability)) {
    component_names <- as.character(q$component_name[probability])
    zero <- probability[match("component_0", component_names)]
    if (is.na(zero)) zero <- probability[1L]
    out[[length(out) + 1L]] <- study06_long_trace(
      nrow(annotations) * (1 - values[, , zero]),
      "expected_active_count", fit_row$fit_id)
  }
  for (name in c("traced_component_0", "traced_component_1",
      "traced_component_2", "traced_component_3", "traced_active_count",
      "causal_active_count", "causal_component_1",
      "causal_component_2", "causal_component_3")) {
    if (!is.null(occupancy)) {
      matrix_value <- matrix(occupancy[[name]], ncol = 4L)
      out[[length(out) + 1L]] <- study06_long_trace(matrix_value,
        paste0("occupancy_", name), fit_row$fit_id)
    }
  }
  if (isTRUE(fit_row$update_alpha)) {
    alpha <- extract_annotation_coefficient_traces(result,
      expected_chains = spec$qualification$nchains)
    if (any(alpha$status != "ok")) stop(alpha$reason[[1L]], call. = FALSE)
    alpha$quantity <- ifelse(alpha$parameter == "alpha",
      paste("alpha", alpha$annotation, alpha$stick, sep = ":"),
      paste("sigmaSqAlpha", alpha$stick, sep = ":"))
    z <- alpha[c("iteration", "chain", "quantity", "value")]
    z$scenario <- "shared_informative_truth"; z$replicate <- 1L
    z$method <- fit_row$fit_id
    out[[length(out) + 1L]] <- z
    prior <- summarise_drawwise_annotation_prior(alpha, annotations,
      spec$controls$simulation$mixture_var, marker_truth = marker_truth,
      retain_marker_summary = FALSE)
    if (!identical(prior$status, "ok")) stop(prior$reason, call. = FALSE)
    for (name in c("expected_active", "mean_prior_enriched",
        "mean_prior_unannotated", "enriched_prior_contrast",
        "continuous_prior_contrast", "mean_prior_causal",
        "mean_prior_noncausal")) {
      z <- prior$draws[c("iteration", "chain", name)]
      names(z)[3L] <- "value"
      z$quantity <- paste0("prior_", name)
      z$scenario <- "shared_informative_truth"; z$replicate <- 1L
      z$method <- fit_row$fit_id
      out[[length(out) + 1L]] <- z[c("iteration", "chain", "quantity",
        "value", "scenario", "replicate", "method")]
    }
  }
  do.call(rbind, out)
}

study06_selected_diagnostics <- function(draws, spec) {
  diagnostics <- .annotation_candidate_diagnostics(draws, spec)
  choice <- .annotation_select_candidate(diagnostics)
  selected <- choice$selected
  rows <- diagnostics[diagnostics$burnin == selected$burnin &
    diagnostics$retained == selected$retained, , drop = FALSE]
  list(selected = selected, rows = rows, all = diagnostics)
}

study06_window_pip <- function(component, burnin, retained) {
  keep <- seq.int(burnin + 1L, burnin + retained)
  total <- numeric(dim(component)[3L])
  chain <- matrix(NA_real_, dim(component)[3L], dim(component)[2L])
  for (i in seq_len(dim(component)[2L])) {
    value <- colMeans(component[keep, i, , drop = TRUE] > 0)
    total <- total + value
    chain[, i] <- value
  }
  list(pip = total / dim(component)[2L], chain_pip = chain)
}

study06_bayesian_fdr <- function(pip, causal, level) {
  order <- order(-pip, seq_along(pip))
  expected_fdp <- cumsum(1 - pip[order]) / seq_along(order)
  eligible <- which(expected_fdp <= level)
  selected <- if (length(eligible)) order[seq_len(max(eligible))] else integer()
  c(selected = length(selected), true_causal = sum(causal[selected]))
}

study06_power_metrics <- function(pip, effect, marker_truth, prediction,
                                  genetic_truth, phenotype, fit_id) {
  causal <- marker_truth$true_nonnull
  ordering <- order(-pip, seq_along(pip))
  ranked <- causal[ordering]
  rank_value <- rank(-pip, ties.method = "average")
  average_precision <- mean(cumsum(ranked)[ranked] / which(ranked))
  r <- rank(pip, ties.method = "average")
  auroc <- (sum(r[causal]) - sum(seq_len(sum(causal)))) /
    (sum(causal) * sum(!causal))
  values <- c(pip_auprc = average_precision, pip_auroc = auroc,
    causal_rank_median = stats::median(rank_value[causal]),
    causal_rank_mean = mean(rank_value[causal]),
    causal_mean_pip = mean(pip[causal]),
    noncausal_mean_pip = mean(pip[!causal]),
    causal_noncausal_pip_ratio = mean(pip[causal]) / mean(pip[!causal]),
    posterior_effect_truth_correlation = stats::cor(effect,
      marker_truth$effect),
    validation_genetic_value_correlation = stats::cor(prediction,
      genetic_truth), phenotype_prediction_correlation = stats::cor(prediction,
      phenotype))
  for (k in c(10L, 25L, 50L, 100L)) {
    selected <- ordering[seq_len(k)]
    values[paste0("recall_", k)] <- sum(causal[selected]) / sum(causal)
    values[paste0("precision_", k)] <- sum(causal[selected]) / k
  }
  for (level in c(.05, .10)) {
    fdr <- study06_bayesian_fdr(pip, causal, level)
    suffix <- if (level == .05) "05" else "10"
    values[paste0("bayes_fdr_", suffix, "_selected")] <- fdr[["selected"]]
    values[paste0("bayes_fdr_", suffix, "_true_causal")] <-
      fdr[["true_causal"]]
  }
  data.frame(fit_id = fit_id, metric = names(values),
    value = as.numeric(values), stringsAsFactors = FALSE)
}

study06_component_recovery <- function(pip, marker_truth, fit_id) {
  do.call(rbind, lapply(1:3, function(component) {
    causal <- marker_truth$component_index == component
    ordering <- order(-pip, seq_along(pip))
    ranked <- causal[ordering]
    data.frame(fit_id = fit_id, true_component = component,
      causal_count = sum(causal), mean_pip = mean(pip[causal]),
      auprc = mean(cumsum(ranked)[ranked] / which(ranked)),
      recall_50 = sum(causal[ordering[1:50]]) / sum(causal),
      recall_100 = sum(causal[ordering[1:100]]) / sum(causal),
      stringsAsFactors = FALSE)
  }))
}

study06_fixed_alpha_audit <- function(result, true_alpha, true_prior,
                                      mixture_var, fit_id) {
  native <- .benchmark_native_fit(result)
  input <- native$input
  initialized <- input$annot_alpha_init %||% input$alpha_init
  alpha <- extract_annotation_coefficient_traces(result, expected_chains = 4L)
  alpha <- alpha[alpha$parameter == "alpha", , drop = FALSE]
  expected <- true_alpha[cbind(match(alpha$annotation, rownames(true_alpha)),
    match(alpha$stick, paste0("component_", 0:2, "_stick")))]
  final <- native$alpha_final[[1L]]
  native_prior <- native$annotation_prior[[1L]]
  implied_prior <- sblr::sbayesrc_marker_pi(native$annotation, final,
    gamma = mixture_var)
  if (!identical(dim(implied_prior), dim(true_prior)) ||
      any(!is.finite(implied_prior)) ||
      any(abs(rowSums(implied_prior) - 1) > 1e-12))
    stop("Fixed-alpha marker probabilities are invalid.", call. = FALSE)
  prior_error <- max(abs(implied_prior - true_prior))
  native_prior_error <- if (is.null(native_prior)) NA_real_ else
    max(abs(native_prior - true_prior))
  numerical_tolerance <- 1e-12
  data.frame(fit_id = fit_id, update_alpha = input$updateAlpha,
    initialization_max_abs_error = max(abs(initialized - true_alpha)),
    retained_trace_max_abs_error = max(abs(alpha$value - expected)),
    final_alpha_max_abs_error = max(abs(final - true_alpha)),
    marker_prior_max_abs_error = prior_error,
    native_marker_prior_available = !is.null(native_prior),
    native_marker_prior_max_abs_error = native_prior_error,
    marker_prior_verification = "recomputed_with_public_sbayesrc_marker_pi",
    pass = identical(input$updateAlpha, FALSE) &&
      max(abs(initialized - true_alpha)) == 0 &&
      max(abs(alpha$value - expected)) == 0 &&
      max(abs(final - true_alpha)) <= numerical_tolerance &&
      prior_error <= numerical_tolerance &&
      (is.null(native_prior) || native_prior_error <= numerical_tolerance),
    stringsAsFactors = FALSE)
}

study06_metric_contrasts <- function(metrics, registry, convergence_status) {
  comparisons <- data.frame(
    comparison = rep(c("informative_vs_baseline",
      "informative_vs_shuffled", "fixed_vs_learned"), each = 2L),
    route = rep(c("bed", "block_eigen"), 3L),
    focal_condition = rep(c("learned_informative", "learned_informative",
      "fixed_true_alpha"), each = 2L),
    reference_condition = rep(c("baseline", "learned_shuffled",
      "learned_informative"), each = 2L), stringsAsFactors = FALSE)
  wanted <- c("pip_auprc", "pip_auroc", paste0("recall_", c(10,25,50,100)),
    "bayes_fdr_05_selected", "bayes_fdr_10_selected",
    "phenotype_prediction_correlation")
  rows <- list()
  for (i in seq_len(nrow(comparisons))) {
    a <- registry$fit_id[registry$route == comparisons$route[i] &
      registry$condition == comparisons$focal_condition[i]]
    b <- registry$fit_id[registry$route == comparisons$route[i] &
      registry$condition == comparisons$reference_condition[i]]
    x <- merge(metrics[metrics$fit_id == a & metrics$metric %in% wanted, ],
      metrics[metrics$fit_id == b & metrics$metric %in% wanted, ],
      by = "metric", suffixes = c("_focal", "_reference"))
    rows[[i]] <- data.frame(comparisons[i, ], metric = x$metric,
      focal_fit = a, reference_fit = b,
      delta = x$value_focal - x$value_reference,
      interpretation = if (convergence_status[[a]] && convergence_status[[b]])
        "credible" else "descriptive under non-converged chains",
      stringsAsFactors = FALSE, row.names = NULL)
  }
  do.call(rbind, rows)
}

study06_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}
