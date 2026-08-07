# Study 06 large information-scale feasibility analysis.
#
# This file reads completed, identity-checked checkpoints. It never launches or
# resumes MCMC.

study06_large_trace_quantity <- function(fit, group, quantity = NULL) {
  meta <- fit$convergence_traces$quantities
  take <- meta$group == group
  if (!is.null(quantity)) take <- take & meta$quantity == quantity
  idx <- which(take)
  if (!length(idx)) return(NULL)
  fit$convergence_traces$values[, , idx, drop = FALSE]
}

study06_large_draw_summary <- function(x) {
  z <- as.numeric(x)
  c(mean = mean(z), sd = stats::sd(z), median = stats::median(z),
    q025 = unname(stats::quantile(z, .025)),
    q975 = unname(stats::quantile(z, .975)))
}

study06_large_prediction_matrix <- function(bundle, effects) {
  out <- matrix(0, nrow = length(bundle$sample_ids), ncol = ncol(effects),
    dimnames = list(bundle$sample_ids, colnames(effects)))
  for (block in unique(bundle$blocks$panel$block_id)) {
    idx <- which(bundle$blocks$panel$block_id == block)
    X <- benchmark_extract_scaled_genotypes(bundle$glist,
      bundle$spec$source$chromosome, bundle$sample_ids,
      bundle$markers$marker_ids[idx])
    out <- out + X %*% effects[idx, , drop = FALSE]
  }
  out
}

study06_large_bfdr_set <- function(pip, level) {
  ordering <- order(-pip, seq_along(pip))
  eligible <- which(cumsum(1 - pip[ordering]) / seq_along(ordering) <= level)
  if (length(eligible)) ordering[seq_len(max(eligible))] else integer()
}

study06_large_analyze_completed <- function(checkpoint_dir, output_dir,
                                             bundle,
                                             spec = study06_large_spec()) {
  ids <- c("E0", "B0", "E2", "B2", "E1", "B1")
  paths <- file.path(checkpoint_dir, paste0("gctb_block_", ids, "_fit.rds"))
  if (!all(file.exists(paths))) stop("All six completed checkpoints are required.")
  saved <- stats::setNames(lapply(paths, readRDS), ids)
  fits <- lapply(saved, `[[`, "fit")
  for (id in ids) {
    row <- bundle$registry[bundle$registry$fit_id == id, , drop = FALSE]
    if (!identical(saved[[id]]$identity$truth_hash, bundle$identities$truth_hash) ||
        !identical(saved[[id]]$identity$marker_hash,
          bundle$identities$marker_order_hash) ||
        !identical(saved[[id]]$identity$controls$chain_seeds,
          spec$seeds$chain)) stop("Completed checkpoint identity mismatch: ", id)
    if (row$route == "block_eigen") {
      contract <- study06_large_block_contract()
      if (!identical(saved[[id]]$identity$block_contract, contract) ||
          !identical(fits[[id]]$input$residual_policy,
            contract$residual_policy)) stop("Block contract mismatch: ", id)
    }
  }

  effects <- do.call(cbind, lapply(fits, function(x) as.numeric(x$bm)))
  pips <- do.call(cbind, lapply(fits, function(x) as.numeric(x$dm)))
  colnames(effects) <- colnames(pips) <- ids
  predictions <- study06_large_prediction_matrix(bundle, effects)

  convergence <- do.call(rbind, lapply(ids, function(id) {
    x <- fits[[id]]$convergence$summary
    x$fit_id <- id
    x$threshold_failure <- (x$rhat_flag %in% TRUE) |
      (x$ess_bulk_flag %in% TRUE) | (x$ess_tail_flag %in% TRUE) |
      (x$mcse_flag %in% TRUE)
    x$diagnostic_unavailable <- x$status %in% c("constant_chain_mismatch",
      "computed_partial")
    x
  }))
  rownames(convergence) <- NULL
  selected <- grepl("^selected_", convergence$group)
  computed <- convergence$status %in% c("computed", "computed_partial")
  overview <- do.call(rbind, lapply(ids, function(id) {
    x <- convergence[convergence$fit_id == id, , drop = FALSE]
    core <- !grepl("^selected_", x$group)
    valid <- core & computed[convergence$fit_id == id] &
      is.finite(x$rhat) & is.finite(x$ess_bulk) &
      is.finite(x$ess_tail) & is.finite(x$mcse_mean_over_sd)
    data.frame(fit_id = id,
      core_threshold_failures = sum(x$threshold_failure & core),
      selected_threshold_failures = sum(x$threshold_failure & !core),
      core_unavailable = sum(x$diagnostic_unavailable & core),
      selected_unavailable = sum(x$diagnostic_unavailable & !core),
      max_rhat_core = if (any(valid)) max(x$rhat[valid]) else NA_real_,
      min_bulk_ess_core = if (any(valid)) min(x$ess_bulk[valid]) else NA_real_,
      min_tail_ess_core = if (any(valid)) min(x$ess_tail[valid]) else NA_real_,
      max_relative_mcse_core = if (any(valid))
        max(x$mcse_mean_over_sd[valid]) else NA_real_,
      core_contract_pass = !any(x$threshold_failure & core),
      stringsAsFactors = FALSE)
  }))

  metric_rows <- vector("list", length(ids))
  model_rows <- vector("list", length(ids))
  component_rows <- list()
  stick_rows <- list()
  prior_rows <- list()
  marker_prior_matrices <- list()
  block_rows <- list()
  phenotype_variance <- stats::var(bundle$truth$phenotype)
  for (ii in seq_along(ids)) {
    id <- ids[[ii]]; fit <- fits[[id]]
    metric_rows[[ii]] <- study06_large_metrics(pips[, id], effects[, id],
      bundle$truth$marker_truth, predictions[, id],
      bundle$truth$genetic_value, id)
    vg <- study06_large_trace_quantity(fit, "vgs")
    ve <- study06_large_trace_quantity(fit, "ves")
    vb <- study06_large_trace_quantity(fit, "vbs")
    h2 <- if (startsWith(id, "B")) vg / phenotype_variance else vg / (vg + ve)
    cc <- study06_large_trace_quantity(fit, "component_count")
    active <- study06_large_trace_quantity(fit, "realized_active_count")
    cp <- fit$component_probabilities[[1L]]
    posterior_expected_occupancy <- sum(1 - cp[, 1L])
    if (grepl("1$|2$", id) && !is.null(fit$annotation_effects)) {
      alpha_plugin <- as.matrix(fit$annotation_effects[[1L]])
      colnames(alpha_plugin) <- colnames(bundle$calibrated$alpha)
      marker_prior <- study06_large_component_probabilities(bundle$annotations,
        alpha_plugin)
    } else {
      marker_prior <- matrix(rep(as.numeric(fit$pi_mean),
        each = nrow(bundle$annotations)), ncol = 4L)
    }
    marker_prior_matrices[[id]] <- marker_prior
    prior_expected_active <- sum(1 - marker_prior[, 1L])
    ss_vg <- study06_large_draw_summary(vg)
    ss_ve <- study06_large_draw_summary(ve)
    ss_vb <- study06_large_draw_summary(vb)
    ss_h2 <- study06_large_draw_summary(h2)
    ss_active <- study06_large_draw_summary(active)
    model_rows[[ii]] <- data.frame(fit_id = id,
      runtime_seconds = saved[[id]]$elapsed_seconds,
      checkpoint_bytes = file.info(paths[[ii]])$size,
      effect_variance_mean = ss_vb[["mean"]],
      genetic_variance_mean = ss_vg[["mean"]],
      genetic_variance_q025 = ss_vg[["q025"]],
      genetic_variance_q975 = ss_vg[["q975"]],
      residual_variance_mean = ss_ve[["mean"]],
      residual_variance_q025 = ss_ve[["q025"]],
      residual_variance_q975 = ss_ve[["q975"]],
      heritability_mean = ss_h2[["mean"]],
      heritability_q025 = ss_h2[["q025"]],
      heritability_q975 = ss_h2[["q975"]],
      realized_active_mean = ss_active[["mean"]],
      realized_active_q025 = ss_active[["q025"]],
      realized_active_q975 = ss_active[["q975"]],
      posterior_expected_active_from_component_probabilities =
        posterior_expected_occupancy,
      annotation_prior_expected_active_plugin = prior_expected_active,
      genetic_value_truth_correlation = stats::cor(predictions[, id],
        bundle$truth$genetic_value),
      phenotype_correlation = stats::cor(predictions[, id],
        bundle$truth$phenotype), stringsAsFactors = FALSE)
    cc_meta <- fit$convergence_traces$quantities[
      fit$convergence_traces$quantities$group == "component_count", ]
    for (j in seq_len(dim(cc)[3L])) {
      ss <- study06_large_draw_summary(cc[, , j])
      component_rows[[length(component_rows) + 1L]] <- data.frame(
        fit_id = id, component = j - 1L,
        truth_count = bundle$truth$component_count[[j]],
        posterior_mean = ss[["mean"]], posterior_sd = ss[["sd"]],
        q025 = ss[["q025"]], q975 = ss[["q975"]],
        absolute_error = ss[["mean"]] - bundle$truth$component_count[[j]],
        relative_error = (ss[["mean"]] - bundle$truth$component_count[[j]]) /
          bundle$truth$component_count[[j]], stringsAsFactors = FALSE)
    }
    for (group in c("stick_eligible_count", "stick_continue_count",
                    "stick_stop_count")) {
      z <- study06_large_trace_quantity(fit, group)
      for (j in seq_len(dim(z)[3L])) {
        ss <- study06_large_draw_summary(z[, , j])
        stick_rows[[length(stick_rows) + 1L]] <- data.frame(fit_id = id,
          quantity = group, stick = j, posterior_mean = ss[["mean"]],
          posterior_sd = ss[["sd"]], q025 = ss[["q025"]],
          q975 = ss[["q975"]], stringsAsFactors = FALSE)
      }
    }
    true_cp <- bundle$calibrated$probability
    prior_active <- 1 - marker_prior[, 1L]
    true_active <- 1 - true_cp[, 1L]
    enriched <- bundle$annotations[, "enriched_binary"] == 1
    prior_rows[[ii]] <- data.frame(fit_id = id,
      marker_prior_pearson = suppressWarnings(stats::cor(prior_active,
        true_active)), marker_prior_spearman = suppressWarnings(stats::cor(
        prior_active, true_active, method = "spearman")),
      expected_active_plugin = sum(prior_active),
      enriched_mean_nonnull = mean(prior_active[enriched]),
      unenriched_mean_nonnull = mean(prior_active[!enriched]),
      enriched_contrast = mean(prior_active[enriched]) -
        mean(prior_active[!enriched]), stringsAsFactors = FALSE)
    if (startsWith(id, "B")) {
      bve <- fit$block_ve
      all_pm <- as.numeric(bve$posterior_mean_per_chain_block)
      block_rows[[length(block_rows) + 1L]] <- data.frame(fit_id = id,
        residual_policy = bve$residual_policy,
        block_ve_mode = bve$block_ve_mode,
        mean_block_ve = mean(all_pm), median_block_ve = stats::median(all_pm),
        min_block_ve = min(all_pm), max_block_ve = max(all_pm),
        resampled = sum(bve$resampled_per_chain_block),
        minimum_ratio_resets = sum(bve$minimum_ratio_resets_per_chain_block),
        history_available = !is.null(bve$history), stringsAsFactors = FALSE)
    }
  }
  metrics <- do.call(rbind, metric_rows)
  model <- do.call(rbind, model_rows)
  component <- do.call(rbind, component_rows)
  stick <- do.call(rbind, stick_rows)
  prior <- do.call(rbind, prior_rows)
  block_ve <- do.call(rbind, block_rows)

  alpha_rows <- list(); sigma_rows <- list()
  annotation_names <- rownames(bundle$calibrated$alpha)
  stick_names <- colnames(bundle$calibrated$alpha)
  for (id in c("E1", "B1")) {
    fit <- fits[[id]]; meta <- fit$convergence_traces$quantities
    for (j in which(meta$group == "annotations" & grepl("^alpha\\[",
        meta$quantity))) {
      z <- fit$convergence_traces$values[, , j]
      ss <- study06_large_draw_summary(z)
      annotation <- meta$annotation_name[[j]]
      stick_name <- meta$stick_name[[j]]
      truth <- bundle$calibrated$alpha[annotation, stick_name]
      diag <- fit$convergence$summary[fit$convergence$summary$quantity ==
        meta$quantity[[j]], , drop = FALSE]
      alpha_rows[[length(alpha_rows) + 1L]] <- data.frame(fit_id = id,
        stick = stick_name, annotation = annotation, truth = truth,
        posterior_mean = ss[["mean"]], posterior_sd = ss[["sd"]],
        posterior_median = ss[["median"]], q025 = ss[["q025"]],
        q975 = ss[["q975"]], bias = ss[["mean"]] - truth,
        absolute_error = abs(ss[["mean"]] - truth),
        standardized_error = (ss[["mean"]] - truth) / ss[["sd"]],
        coverage = ss[["q025"]] <= truth && ss[["q975"]] >= truth,
        sign_recovery = truth == 0 || sign(ss[["mean"]]) == sign(truth),
        rhat = diag$rhat, bulk_ess = diag$ess_bulk,
        tail_ess = diag$ess_tail,
        relative_mcse = diag$mcse_mean_over_sd,
        pass = !diag$rhat_flag && !diag$ess_bulk_flag &&
          !diag$ess_tail_flag && !diag$mcse_flag,
        stringsAsFactors = FALSE)
    }
    for (j in which(meta$group == "annotations" &
        grepl("^sigmaSqAlpha\\[", meta$quantity))) {
      z <- fit$convergence_traces$values[, , j]
      ss <- study06_large_draw_summary(z)
      diag <- fit$convergence$summary[fit$convergence$summary$quantity ==
        meta$quantity[[j]], , drop = FALSE]
      sigma_rows[[length(sigma_rows) + 1L]] <- data.frame(fit_id = id,
        stick = meta$stick_name[[j]], registered_reference = 1,
        posterior_mean = ss[["mean"]], posterior_sd = ss[["sd"]],
        posterior_median = ss[["median"]], q025 = ss[["q025"]],
        q975 = ss[["q975"]], rhat = diag$rhat,
        bulk_ess = diag$ess_bulk, tail_ess = diag$ess_tail,
        relative_mcse = diag$mcse_mean_over_sd,
        pass = !diag$rhat_flag && !diag$ess_bulk_flag &&
          !diag$ess_tail_flag && !diag$mcse_flag,
        stringsAsFactors = FALSE)
    }
  }
  alpha <- do.call(rbind, alpha_rows)
  sigma <- do.call(rbind, sigma_rows)

  pair_rows <- list()
  for (pair in list(c("E0", "B0"), c("E2", "B2"), c("E1", "B1"))) {
    a <- pair[[1L]]; b <- pair[[2L]]
    top50a <- order(-pips[, a])[1:50]; top50b <- order(-pips[, b])[1:50]
    top100a <- order(-pips[, a])[1:100]; top100b <- order(-pips[, b])[1:100]
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(pair = paste(a, b,
      sep = "_vs_"), pip_pearson = stats::cor(pips[, a], pips[, b]),
      pip_spearman = stats::cor(pips[, a], pips[, b], method = "spearman"),
      effect_pearson = stats::cor(effects[, a], effects[, b]),
      validation_g_pearson = stats::cor(predictions[, a], predictions[, b]),
      top50_overlap = length(intersect(top50a, top50b)),
      top100_overlap = length(intersect(top100a, top100b)),
      marker_prior_pearson = suppressWarnings(stats::cor(
        1 - marker_prior_matrices[[a]][, 1L],
        1 - marker_prior_matrices[[b]][, 1L])),
      heritability_difference_bed_minus_block =
        model$heritability_mean[model$fit_id == a] -
          model$heritability_mean[model$fit_id == b],
      genetic_variance_difference_bed_minus_block =
        model$genetic_variance_mean[model$fit_id == a] -
          model$genetic_variance_mean[model$fit_id == b],
      stringsAsFactors = FALSE)
  }
  route <- do.call(rbind, pair_rows)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_one <- function(x, name) utils::write.csv(x,
    file.path(output_dir, name), row.names = FALSE, na = "")
  write_one(convergence, "convergence_all.csv")
  write_one(convergence[convergence$threshold_failure |
    convergence$diagnostic_unavailable, ], "convergence_failures.csv")
  write_one(overview, "convergence_overview.csv")
  write_one(model, "model_summary.csv")
  write_one(component, "component_recovery.csv")
  write_one(stick, "stick_count_summary.csv")
  write_one(prior, "marker_prior_recovery.csv")
  write_one(metrics, "snp_metrics.csv")
  write_one(alpha, "alpha_recovery.csv")
  write_one(sigma, "sigmaSqAlpha_summary.csv")
  write_one(route, "route_comparison.csv")
  write_one(block_ve, "block_ve_summary.csv")
  summary_path <- file.path(output_dir, "posterior_scientific_summaries.rds")
  saveRDS(list(effects = effects, pips = pips, genetic_values = predictions),
    summary_path, compress = FALSE)
  analysis_files <- setdiff(list.files(output_dir, full.names = TRUE),
    file.path(output_dir, "manifest.json"))
  manifest <- list(
    schema = "sblrbench-study06-large-feasibility-completed-v1",
    decision = "LARGE-G2",
    secondary_flags = c("LARGE-G3", "LARGE-G4"),
    identities = bundle$identities,
    sblrbench_sha = trimws(system2("git", "rev-parse HEAD", stdout = TRUE)),
    sblr_sha = saved[[1L]]$identity$sblr_sha,
    installed_package = saved[[1L]]$identity$installed_package,
    chain_seeds = as.list(spec$seeds$chain),
    checkpoint_semantic_hashes = stats::setNames(lapply(saved,
      `[[`, "semantic_hash"), ids),
    checkpoint_sha256 = stats::setNames(lapply(paths, digest::digest,
      file = TRUE, algo = "sha256"), basename(paths)),
    analysis_sha256 = stats::setNames(lapply(analysis_files, digest::digest,
      file = TRUE, algo = "sha256"), basename(analysis_files)),
    runtime_seconds = stats::setNames(as.list(vapply(saved,
      `[[`, numeric(1), "elapsed_seconds")), ids),
    additional_replicates = 0L,
    changed_seed_retries = 0L,
    formal_qualification_changed = FALSE,
    final_benchmark_authorized = FALSE,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(manifest, file.path(output_dir, "manifest.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = 16)
  list(saved = saved, fits = fits, effects = effects, pips = pips,
    predictions = predictions, convergence = convergence, overview = overview,
    model = model, component = component, stick = stick, prior = prior,
    metrics = metrics, alpha = alpha, sigma = sigma, route = route,
    block_ve = block_ve, paths = paths)
}
