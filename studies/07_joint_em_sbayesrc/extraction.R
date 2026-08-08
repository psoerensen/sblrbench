study07_chain_records <- function(fit) {
  chains <- fit$chains
  if (is.null(chains)) return(list())
  out <- list()
  for (item in chains) {
    if (is.list(item) && !is.null(item$chain_index)) {
      out[[length(out) + 1L]] <- item
    } else if (is.list(item)) {
      for (candidate in item)
        if (is.list(candidate)) out[[length(out) + 1L]] <- candidate
    }
  }
  out
}

study07_genomic_fit <- function(result) {
  if (!is.null(result$genomic)) result$genomic else result
}

study07_marker_vector <- function(x, field, fallback = NULL) {
  value <- x[[field]]
  if (is.null(value) && !is.null(x$marker)) value <- x$marker[[field]]
  if (is.null(value)) value <- fallback
  if (is.list(value) && length(value) == 1L) value <- value[[1L]]
  if (is.matrix(value)) value <- value[, 1L]
  as.numeric(value)
}

study07_component_probability <- function(fit) {
  value <- fit$component_probabilities
  if (is.null(value) && !is.null(fit$component$prob))
    value <- fit$component$prob[[1L]]
  if (is.list(value) && length(value) == 1L) value <- value[[1L]]
  if (is.null(value)) return(NULL)
  as.matrix(value)
}

study07_variance_draws <- function(fit, field) {
  records <- study07_chain_records(fit)
  values <- unlist(lapply(records, function(chain) {
    value <- chain[[field]]
    if (is.null(value) && !is.null(chain$trace)) value <- chain$trace[[field]]
    value <- as.numeric(value)
    if (length(value) && !is.null(fit$meta$nit))
      value <- tail(value, fit$meta$nit)
    value
  }), use.names = FALSE)
  if (length(values)) return(values)
  value <- fit[[field]]
  if (is.null(value) && !is.null(fit$trace)) value <- fit$trace[[field]]
  as.numeric(value)
}

study07_genomic_values <- function(payload, method, inputs) {
  fit <- study07_genomic_fit(payload$result)
  beta <- study07_marker_vector(fit, "bm")
  pip <- study07_marker_vector(fit, "dm")
  m <- length(inputs$truth$effects)
  study07_assert(length(beta) == m && length(pip) == m &&
    all(is.finite(beta)) && all(is.finite(pip)),
    paste(method, "does not expose aligned posterior beta/PIP summaries."))
  component_probability <- study07_component_probability(fit)
  if (!is.null(component_probability))
    study07_assert(identical(dim(component_probability), c(m, 4L)),
      paste(method, "component-probability dimensions changed."))
  list(beta = beta, pip = pip, component_probability = component_probability,
    fit = fit)
}

study07_genomic_metrics <- function(payload, method, inputs,
                                    start = "primary") {
  value <- study07_genomic_values(payload, method, inputs)
  truth <- as.numeric(inputs$truth$effects)
  causal <- as.logical(inputs$marker_truth$true_nonnull)
  component <- as.integer(inputs$marker_truth$component_index) + 1L
  probability <- value$component_probability
  true_component_probability <- component_accuracy <- NA_real_
  if (!is.null(probability)) {
    true_component_probability <- mean(probability[cbind(seq_len(nrow(probability)),
      component)])
    component_accuracy <- mean(max.col(probability, ties.method = "first") ==
      component)
  }
  trace_bundle <- value$fit$convergence_traces
  trace_group <- function(group) {
    if (is.null(trace_bundle$values)) return(numeric())
    index <- match(group, trace_bundle$quantities$group)
    if (is.na(index)) numeric() else as.numeric(trace_bundle$values[, , index])
  }
  vb <- trace_group("vbs"); vg <- trace_group("vgs"); ve <- trace_group("ves")
  if (!length(vb)) vb <- study07_variance_draws(value$fit, "vbs")
  if (!length(vg)) vg <- study07_variance_draws(value$fit, "vgs")
  if (!length(ve)) ve <- study07_variance_draws(value$fit, "ves")
  if (!length(vb) && !is.null(payload$result$mcem))
    vb <- as.numeric(payload$result$mcem$genomic_hyperparameters$B_final)
  if (!length(ve) && !is.null(payload$result$mcem))
    ve <- as.numeric(payload$result$mcem$genomic_hyperparameters$E_final)
  h2 <- if (length(vg) && length(ve)) mean(vg / (vg + ve)) else NA_real_
  validation_genetic <- as.numeric(inputs$truth$validation_x %*% value$beta)
  data.frame(
    method = method, start = start,
    beta_cor = stats::cor(value$beta, truth),
    beta_rmse = sqrt(mean((value$beta - truth)^2)),
    pip_brier = mean((value$pip - causal)^2),
    mean_pip_causal = mean(value$pip[causal]),
    mean_pip_null = mean(value$pip[!causal]),
    true_component_probability = true_component_probability,
    component_accuracy = component_accuracy,
    B = if (length(vb)) mean(vb) else NA_real_,
    E = if (length(ve)) mean(ve) else NA_real_,
    heritability = h2,
    expected_active_count = sum(value$pip),
    validation_genetic_cor = stats::cor(validation_genetic,
      inputs$truth$validation_genetic),
    validation_prediction_cor = stats::cor(validation_genetic,
      inputs$truth$validation_phenotype),
    runtime_seconds = payload$elapsed_seconds,
    stringsAsFactors = FALSE)
}

study07_alpha_trace <- function(fit) {
  bundle <- fit$convergence_traces
  if (is.null(bundle$values) || is.null(bundle$quantities)) return(NULL)
  index <- which(bundle$quantities$group %in% c("vbs", "vgs", "ves",
      "realized_active_count") |
    bundle$quantities$parameter_name %in% c("alpha", "sigmaSqAlpha"))
  if (!length(index)) return(NULL)
  values <- bundle$values[, , index, drop = FALSE]
  ids <- bundle$quantities$quantity[index]
  benchmark_trace_array_long(values, ids,
    expected_chains = dim(values)[2L])
}

study07_joint_diagnostics <- function(payload, method, spec) {
  trace <- study07_alpha_trace(payload$result)
  if (is.null(trace) && identical(method, "SBayesRC-S"))
    trace <- study07_selection_joint_trace(payload$result)
  if (is.null(trace)) return(data.frame(method = method,
    quantity = "unavailable", rhat = NA_real_, ess_bulk = NA_real_,
    ess_tail = NA_real_, mcse_mean = NA_real_, relative_mcse = NA_real_,
    stringsAsFactors = FALSE))
  rows <- lapply(unique(trace$quantity), function(quantity) {
    x <- trace[trace$quantity == quantity, , drop = FALSE]
    diagnostic <- tryCatch(benchmark_scalar_diagnostics(x,
      spec$joint$convergence_thresholds), error = function(e) data.frame(
        rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
        mcse_mean = NA_real_, relative_mcse = NA_real_))
    cbind(method = method, quantity = quantity,
      diagnostic[c("rhat", "ess_bulk", "ess_tail", "mcse_mean",
        "relative_mcse")], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

study07_selection_joint_trace <- function(fit) {
  records <- study07_chain_records(fit)
  if (length(records) != 4L) return(NULL)
  pieces <- list()
  add <- function(field, matrices) {
    if (any(vapply(matrices, is.null, logical(1)))) return(NULL)
    matrices <- lapply(matrices, as.matrix)
    if (length(unique(vapply(matrices, nrow, integer(1)))) != 1L ||
        length(unique(vapply(matrices, ncol, integer(1)))) != 1L)
      return(NULL)
    values <- array(NA_real_, c(nrow(matrices[[1L]]), length(matrices),
      ncol(matrices[[1L]])))
    for (chain in seq_along(matrices)) values[, chain, ] <- matrices[[chain]]
    ids <- if (ncol(matrices[[1L]]) == 1L) field else
      paste0(field, "[", seq_len(ncol(matrices[[1L]])), "]")
    benchmark_trace_array_long(values, ids,
      expected_chains = length(matrices))
  }
  for (field in c("alpha", "sigmaSqAlpha", "annotation_delta",
      "annotation_pi_A", "annotation_included_count",
      "realized_active_count")) {
    value <- add(field, lapply(records, function(x)
      x$convergence_trace[[field]]))
    if (!is.null(value)) pieces[[field]] <- value
  }
  for (field in c("vbs", "vgs", "ves")) {
    value <- add(field, lapply(records, function(x) {
      z <- as.numeric(x$trace[[field]])
      matrix(tail(z, fit$meta$nit), ncol = 1L)
    }))
    if (!is.null(value)) pieces[[field]] <- value
  }
  if (!length(pieces)) NULL else do.call(rbind, pieces)
}

study07_prior_metrics <- function(prior, inputs, method, start) {
  prior <- as.matrix(prior)
  truth <- inputs$prior_truth
  study07_assert(identical(dim(prior), dim(truth)) &&
    max(abs(rowSums(prior) - 1)) < 1e-10,
    paste(method, "induced component priors are invalid."))
  data.frame(method = method, start = start,
    prior_rmse = sqrt(mean((prior - truth)^2)),
    prior_cor = stats::cor(as.numeric(prior), as.numeric(truth)),
    active_prior_cor = stats::cor(1 - prior[, 1L], 1 - truth[, 1L]),
    active_prior_rmse = sqrt(mean(((1 - prior[, 1L]) -
      (1 - truth[, 1L]))^2)),
    stringsAsFactors = FALSE)
}

study07_joint_annotation_metrics <- function(payload, inputs,
                                               prior_summary = NULL) {
  fit <- payload$result
  trace <- fit$convergence_traces
  index <- which(trace$quantities$parameter_name == "alpha")
  values <- trace$values[, , index, drop = FALSE]
  alpha_mean <- matrix(colMeans(matrix(values,
    ncol = length(index))), nrow = 4L, ncol = 3L)
  dimnames(alpha_mean) <- dimnames(inputs$alpha_truth)
  if (is.null(prior_summary)) prior_summary <-
    study07_joint_prior_summaries(fit, inputs$truth$annotations)
  prior <- prior_summary$pooled
  alpha_truth <- as.numeric(inputs$alpha_truth)
  cbind(data.frame(method = "SBayesRC", start = "primary",
    alpha_rmse = sqrt(mean((as.numeric(alpha_mean) - alpha_truth)^2)),
    alpha_cor = stats::cor(as.numeric(alpha_mean), alpha_truth),
    stringsAsFactors = FALSE),
    study07_prior_metrics(prior, inputs, "SBayesRC", "primary")[-c(1, 2)])
}

study07_joint_prior_summaries <- function(fit, annotation) {
  trace <- fit$convergence_traces
  descriptor <- trace$quantities
  index <- which(descriptor$parameter_name == "alpha")
  study07_assert(length(index) == ncol(annotation) * 3L,
    "Joint SBayesRC lacks complete retained alpha histories.")
  descriptor <- descriptor[index, , drop = FALSE]
  values <- trace$values[, , index, drop = FALSE]
  chains <- lapply(seq_len(dim(values)[2L]), function(chain) {
    q <- lapply(seq_len(3L), function(stick) {
      columns <- which(descriptor$stick_index == stick)
      columns <- columns[order(descriptor$annotation_index[columns])]
      alpha <- values[, chain, columns, drop = FALSE]
      dim(alpha) <- c(dim(values)[1L], length(columns))
      stats::pnorm(annotation %*% t(alpha))
    })
    cbind(
      component_0 = rowMeans(1 - q[[1L]]),
      component_1 = rowMeans(q[[1L]] * (1 - q[[2L]])),
      component_2 = rowMeans(q[[1L]] * q[[2L]] * (1 - q[[3L]])),
      component_3 = rowMeans(q[[1L]] * q[[2L]] * q[[3L]]))
  })
  pooled <- Reduce(`+`, chains) / length(chains)
  study07_assert(max(abs(rowSums(pooled) - 1)) < 1e-10,
    "Retained-alpha marker priors do not normalize.")
  list(chains = chains, pooled = pooled,
    source = "genuine retained alpha histories")
}

study07_joint_prior_stability <- function(fit, annotation,
                                           prior_summary = NULL) {
  if (is.null(prior_summary)) prior_summary <-
    study07_joint_prior_summaries(fit, annotation)
  prior <- prior_summary$chains
  pairs <- utils::combn(seq_along(prior), 2L, simplify = FALSE)
  do.call(rbind, lapply(c("component_prior", "active_prior"), function(target) {
    metric <- do.call(rbind, lapply(pairs, function(pair) {
      x <- if (target == "component_prior") as.numeric(prior[[pair[[1L]]]]) else
        1 - prior[[pair[[1L]]]][, 1L]
      y <- if (target == "component_prior") as.numeric(prior[[pair[[2L]]]]) else
        1 - prior[[pair[[2L]]]][, 1L]
      data.frame(pearson = stats::cor(x, y),
        spearman = stats::cor(x, y, method = "spearman"),
        rmse = sqrt(mean((x - y)^2)), max_abs = max(abs(x - y)))
    }))
    data.frame(method = "SBayesRC", quantity = target,
      minimum_pearson = min(metric$pearson),
      minimum_spearman = min(metric$spearman),
      maximum_rmse = max(metric$rmse),
      maximum_absolute_difference = max(metric$max_abs),
      source = "genuine retained alpha histories")
  }))
}

study07_em_annotation_metrics <- function(payload, inputs, start) {
  mcem <- payload$result$mcem
  alpha_truth <- as.numeric(inputs$alpha_truth)
  cbind(data.frame(method = "SBayesRC-EM", start = start,
    alpha_rmse = sqrt(mean((as.numeric(mcem$alpha_map) - alpha_truth)^2)),
    alpha_cor = stats::cor(as.numeric(mcem$alpha_map), alpha_truth),
    stringsAsFactors = FALSE),
    study07_prior_metrics(mcem$component_prior, inputs,
      "SBayesRC-EM", start)[-c(1, 2)])
}

study07_pairwise_max <- function(values) {
  if (length(values) < 2L) return(0)
  max(vapply(utils::combn(seq_along(values), 2L, simplify = FALSE),
    function(i) max(abs(values[[i[[1L]]]] - values[[i[[2L]]]])), numeric(1)))
}

study07_pairwise_min_cor <- function(values, method = "pearson") {
  if (length(values) < 2L || length(values[[1L]]) < 2L) return(NA_real_)
  correlations <- vapply(utils::combn(seq_along(values), 2L,
      simplify = FALSE), function(i) suppressWarnings(stats::cor(
    as.numeric(values[[i[[1L]]]]), as.numeric(values[[i[[2L]]]]),
    method = method)), numeric(1))
  if (all(!is.finite(correlations))) NA_real_ else
    min(correlations[is.finite(correlations)])
}

study07_em_stability <- function(payloads, method, inputs) {
  mcem <- lapply(payloads, function(x) x$result$mcem)
  genomic <- lapply(payloads, function(x) study07_genomic_values(x, method,
    inputs))
  values <- list(
    alpha_map = lapply(mcem, `[[`, "alpha_map"),
    component_prior = lapply(mcem, `[[`, "component_prior"),
    active_prior = lapply(mcem, function(x) 1 - x$component_prior[, 1L]),
    genomic_pip = lapply(genomic, `[[`, "pip"),
    genomic_beta = lapply(genomic, `[[`, "beta"),
    B = lapply(mcem, function(x) as.numeric(
      x$genomic_hyperparameters$B_final)),
    E = lapply(mcem, function(x) as.numeric(
      x$genomic_hyperparameters$E_final)))
  if (identical(method, "SBayesRC-S-EM")) {
    values$annotation_pip_eb <- lapply(mcem, `[[`, "annotation_pip_eb")
    values$delta_map <- lapply(mcem, `[[`, "delta_map")
    values$alpha_model_average <- lapply(mcem, `[[`,
      "alpha_model_average")
  }
  terminal <- lapply(mcem, function(x) tail(x$history$summary, 1L))
  do.call(rbind, lapply(names(values), function(quantity) data.frame(
    method = method, quantity = quantity,
    between_start_max_abs = study07_pairwise_max(values[[quantity]]),
    between_start_min_pearson = study07_pairwise_min_cor(values[[quantity]]),
    between_start_min_spearman = study07_pairwise_min_cor(values[[quantity]],
      "spearman"),
    all_converged = all(vapply(mcem, function(x) isTRUE(x$converged),
      logical(1))),
    outer_iterations_min = min(vapply(mcem, `[[`, integer(1), "n_outer")),
    outer_iterations_max = max(vapply(mcem, `[[`, integer(1), "n_outer")),
    terminal_max_delta_alpha = max(vapply(terminal,
      function(x) x$max_delta_alpha, numeric(1))),
    terminal_max_delta_prior = max(vapply(terminal,
      function(x) x$max_delta_prior, numeric(1))),
    stringsAsFactors = FALSE)))
}

study07_selection_table <- function(joint_payload, em_payloads, spec, inputs) {
  records <- study07_chain_records(joint_payload$result)
  joint_chain_pip <- do.call(rbind, lapply(records, function(chain)
    as.numeric(chain$annotation$annotation_pip)))
  em_pip <- do.call(rbind, lapply(em_payloads, function(x)
    as.numeric(x$result$mcem$annotation_pip_eb)))
  em_delta <- do.call(rbind, lapply(em_payloads, function(x)
    as.integer(x$result$mcem$delta_map)))
  primary <- em_payloads[["baseline"]]$result$mcem
  annotations <- names(spec$model$selection_truth)
  data.frame(annotation = annotations,
    true_role = as.integer(spec$model$selection_truth),
    joint_annotation_pip = colMeans(joint_chain_pip),
    joint_annotation_pip_min = apply(joint_chain_pip, 2L, min),
    joint_annotation_pip_max = apply(joint_chain_pip, 2L, max),
    em_annotation_pip_eb = as.numeric(primary$annotation_pip_eb),
    em_annotation_pip_eb_min = apply(em_pip, 2L, min),
    em_annotation_pip_eb_max = apply(em_pip, 2L, max),
    em_delta_map = as.integer(primary$delta_map),
    em_delta_map_patterns = apply(em_delta, 2L,
      function(x) paste(x, collapse = "/")),
    em_alpha_model_average_stick1 =
      as.numeric(primary$alpha_model_average[-1L, 1L]),
    stringsAsFactors = FALSE)
}

study07_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

study07_read_payload <- function(spec, method, start = "primary")
  readRDS(study07_checkpoint_path(spec, method, start))

study07_summarize <- function(spec, inputs) {
  starts <- spec$em$starts
  payload <- list(
    SBayesR = study07_read_payload(spec, "SBayesR"),
    SBayesRC = study07_read_payload(spec, "SBayesRC"),
    `SBayesRC-S` = study07_read_payload(spec, "SBayesRC-S"),
    `SBayesRC-EM` = setNames(lapply(starts, function(x)
      study07_read_payload(spec, "SBayesRC-EM", x)), starts),
    `SBayesRC-S-EM` = setNames(lapply(starts, function(x)
      study07_read_payload(spec, "SBayesRC-S-EM", x)), starts))
  genomic <- rbind(
    study07_genomic_metrics(payload$SBayesR, "SBayesR", inputs),
    study07_genomic_metrics(payload$SBayesRC, "SBayesRC", inputs),
    study07_genomic_metrics(payload[["SBayesRC-S"]], "SBayesRC-S", inputs),
    do.call(rbind, lapply(starts, function(start)
      study07_genomic_metrics(payload[["SBayesRC-EM"]][[start]],
        "SBayesRC-EM", inputs, start))),
    do.call(rbind, lapply(starts, function(start)
      study07_genomic_metrics(payload[["SBayesRC-S-EM"]][[start]],
        "SBayesRC-S-EM", inputs, start))))
  joint_prior <- study07_joint_prior_summaries(payload$SBayesRC$result,
    inputs$truth$annotations)
  annotation <- rbind(study07_joint_annotation_metrics(payload$SBayesRC,
    inputs, joint_prior), do.call(rbind, lapply(starts, function(start)
      study07_em_annotation_metrics(payload[["SBayesRC-EM"]][[start]],
        inputs, start))))
  joint_stability <- rbind(
    study07_joint_diagnostics(payload$SBayesR, "SBayesR", spec),
    study07_joint_diagnostics(payload$SBayesRC, "SBayesRC", spec),
    study07_joint_diagnostics(payload[["SBayesRC-S"]], "SBayesRC-S", spec))
  prior_stability <- study07_joint_prior_stability(payload$SBayesRC$result,
    inputs$truth$annotations, joint_prior)
  em_stability <- rbind(
    study07_em_stability(payload[["SBayesRC-EM"]], "SBayesRC-EM", inputs),
    study07_em_stability(payload[["SBayesRC-S-EM"]],
      "SBayesRC-S-EM", inputs))
  selection <- study07_selection_table(payload[["SBayesRC-S"]],
    payload[["SBayesRC-S-EM"]], spec, inputs)
  runtime <- rbind(
    data.frame(method = c("SBayesR", "SBayesRC", "SBayesRC-S"),
      start = "primary", runtime_seconds = c(payload$SBayesR$elapsed_seconds,
        payload$SBayesRC$elapsed_seconds,
        payload[["SBayesRC-S"]]$elapsed_seconds)),
    do.call(rbind, lapply(c("SBayesRC-EM", "SBayesRC-S-EM"), function(method)
      data.frame(method = method, start = starts,
        runtime_seconds = vapply(payload[[method]], `[[`, numeric(1),
          "elapsed_seconds")))))
  tables <- file.path(spec$output$local_dir, "tables")
  files <- c(
    genomic_summary = study07_write_csv(genomic,
      file.path(tables, "genomic_summary.csv")),
    annotation_summary = study07_write_csv(annotation,
      file.path(tables, "annotation_summary.csv")),
    joint_stability = study07_write_csv(joint_stability,
      file.path(tables, "joint_stability.csv")),
    prior_stability = study07_write_csv(prior_stability,
      file.path(tables, "joint_prior_stability.csv")),
    em_stability = study07_write_csv(em_stability,
      file.path(tables, "em_stability.csv")),
    annotation_selection = study07_write_csv(selection,
      file.path(tables, "annotation_selection.csv")),
    runtime = study07_write_csv(runtime, file.path(tables, "runtime.csv")))
  list(payload = payload, joint_prior = joint_prior, genomic = genomic,
    joint_stability = joint_stability, prior_stability = prior_stability,
    em_stability = em_stability,
    selection = selection, runtime = runtime, files = files)
}
