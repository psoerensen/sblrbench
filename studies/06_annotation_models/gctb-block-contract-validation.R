# Study 06 pinned GCTB-compatible block residual-contract validation.

study06_gctb_block_constants <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  list(
    schema = "sblrbench-study06-gctb-block-contract-v1",
    specification_hash = "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56",
    truth_hash = "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb",
    official_sha = "b95d3fcbad8ff358290922a58fff893439296138",
    gctb_sha = "cc7fa7d765c83a89c6375946cf77fe50ba1a317e",
    sblr_sha = "0c89234273389e14112ba0e08ef9d11d3e1819dc",
    sblr_version = "0.2.0",
    installed_tree_sha256 = "e723528e7d5d570a31b5b1d1c90551896ac48f86ab05261c181c8109af971fd0",
    output = file.path(root, "results", "local", "06_annotation_models",
      "gctb_block_contract_validation"),
    official_output = file.path(root, "results", "local",
      "06_annotation_models", "gctb_single_trajectory"),
    export = file.path(root, "results", "local", "06_annotation_models",
      "gctb_parity", "export"),
    niter = 9000L, burn = 3000L, retained = 6000L,
    seeds = c(S0_new = 711121L, S1_new = 721121L),
    residual_policy = "gctb_block", block_ve_mode = "allMixVe",
    resam_thresh = 1.1, minimum_ve_ratio = 0.7,
    gates = list(
      S0_new = c(pip_pearson = .99, effect_pearson = .995,
        validation_g_pearson = .995, mean_block_ve_abs_diff = .005,
        heritability_abs_diff = .05),
      S1_new = c(pip_pearson = .95, effect_pearson = .995,
        validation_g_pearson = .995, mean_block_ve_abs_diff = .005,
        heritability_abs_diff = .05))
  )
}

study06_gctb_block_registry <- function(cfg = study06_gctb_block_constants()) {
  data.frame(
    fit_id = c("S0_new", "S1_new"), official_id = c("D0", "D1"),
    method_id = c("st_block_eigen_sbayesr", "st_block_eigen_sbayesrc"),
    annotation_aware = c(FALSE, TRUE), update_alpha = c(FALSE, TRUE),
    seed = unname(cfg$seeds), niter = cfg$niter, burn = cfg$burn,
    residual_policy = cfg$residual_policy,
    block_ve_mode = cfg$block_ve_mode, resam_thresh = cfg$resam_thresh,
    minimum_ve_ratio = cfg$minimum_ve_ratio, stringsAsFactors = FALSE)
}

study06_gctb_block_controls <- function(spec, row, annotation_truth) {
  controls <- study06_power_controls(spec, character(), row$seed, row$seed,
    row$annotation_aware, row$update_alpha)
  # sblr's public `nit` is the post-burn retained count; official niter is the
  # total trajectory length. Preserve 9,000 total = 3,000 burn + 6,000 kept.
  controls$nit <- row$niter - row$burn
  controls$nburn <- row$burn
  controls$nthin <- 1L
  controls$nchains <- 1L
  controls$ncores <- 1L
  controls$chain_seeds <- as.integer(row$seed)
  controls$residual_policy <- row$residual_policy
  controls$block_ve_mode <- row$block_ve_mode
  controls$resam_thresh <- row$resam_thresh
  controls$minimum_ve_ratio <- row$minimum_ve_ratio
  controls$block_ve_keep_history <- TRUE
  controls$convergence_control$aggregate_component_states <- TRUE
  controls
}

study06_gctb_block_fit <- function(method, controls, simulation, stats, glist,
                                   split, annotations, annotation_truth,
                                   block_start) {
  if (isTRUE(method$annotation_aware)) {
    controls$alpha_init <- annotation_truth$uninformative_annotations
  } else {
    controls$pi <- annotation_truth$marginal_component_probability
  }
  controls$diagnostic_mode <- NULL
  controls <- .annotation_public_api_controls(method, controls)
  inputs <- list(stats = stats, Glist = glist,
    block_start = as.integer(block_start),
    representation = method$operator_representation,
    eigen_policy = method$eigen_policy, eigen_prop = method$eigen_prop)
  if (isTRUE(method$annotation_aware)) inputs$annotation <- annotations
  capabilities <- c("posterior_effects", "component_probabilities",
    "scalar_trait", "summary_statistics",
    if (isTRUE(method$annotation_aware)) "annotation_coefficients")
  method_spec <- new_sblr_native_method(method$id, method$label,
    method$interface, method$native_method, capabilities = capabilities,
    metadata = list(task = "annotation_models",
      study_version = simulation$extras$study_version,
      annotation_aware = isTRUE(method$annotation_aware),
      representation = method$representation))
  result <- run_sblrbench_method(method_spec, fit_inputs = inputs,
    controls = controls)
  validate_sblrbench_result(result, simulation)
  result
}

study06_gctb_block_extract <- function(result, marker_ids, truth) {
  native <- result$native_fit
  effect <- extract_marker_effects(result)[, 1L]
  effect <- as.numeric(effect[match(marker_ids, names(effect))])
  pip <- extract_marker_probabilities(result)$posterior_inclusion
  if (is.matrix(pip)) pip <- pip[, 1L]
  pip <- as.numeric(pip[match(marker_ids, names(pip))])
  if (any(!is.finite(effect)) || any(!is.finite(pip)))
    stop("SNP outputs are incomplete after immutable marker alignment.",
      call. = FALSE)
  ve_history <- native$block_ve$history
  if (is.list(ve_history)) ve_history <- ve_history[[1L]]
  if (!is.null(ve_history) && nrow(as.matrix(ve_history)) > 6000L)
    ve_history <- tail(as.matrix(ve_history), 6000L)
  mean_ve <- if (!is.null(ve_history)) mean(ve_history) else
    mean(native$block_ve$posterior_mean)
  median_ve <- if (!is.null(ve_history)) stats::median(ve_history) else
    stats::median(native$block_ve$posterior_mean)
  ve_range <- if (!is.null(ve_history)) range(ve_history) else
    range(native$block_ve$posterior_mean)
  h2 <- tail(as.numeric(native$heritability_summary), 6000L)
  chain_trace <- native$chains[[1L]]$convergence_trace
  active <- if (!is.null(chain_trace$realized_active_count))
    mean(chain_trace$realized_active_count) else sum(native$dm)
  component <- if (!is.null(chain_trace$component_count))
    colMeans(chain_trace$component_count) else as.numeric(native$ncomp)
  list(result = result, native = native, effect = effect, pip = pip,
    prediction = drop(truth$validation_x %*% effect),
    mean_block_ve = mean_ve, median_block_ve = median_ve,
    block_ve_min = ve_range[[1L]], block_ve_max = ve_range[[2L]],
    heritability = mean(h2), active_count = active,
    component_count = component,
    resamples = sum(native$block_ve$resampled_per_chain_block),
    resets = sum(native$block_ve$minimum_ratio_resets_per_chain_block))
}

study06_gctb_block_official <- function(id, cfg, truth) {
  prefix <- file.path(cfg$official_output, "runs", id, id)
  fit <- readRDS(paste0(prefix, ".rds"))
  snp <- data.table::fread(paste0(prefix, ".txt"), data.table = FALSE)
  index <- match(truth$marker_ids, snp$SNP)
  if (anyNA(index) || !identical(as.character(snp$SNP[index]),
      truth$marker_ids)) stop("Official SNP alignment failed.", call. = FALSE)
  vy <- stats::var(truth$training_y)
  effect <- as.numeric(snp$BETA[index]) * sqrt(vy)
  keep <- (cfg$burn + 1L):cfg$niter
  block_ve <- as.matrix(fit$vare_hist) * vy
  h2 <- as.numeric(fit$hsq_hist[keep])
  list(effect = effect, pip = as.numeric(snp$PIP[index]),
    prediction = drop(truth$validation_x %*% effect),
    mean_block_ve = mean(block_ve), median_block_ve = median(block_ve),
    block_ve_min = min(block_ve), block_ve_max = max(block_ve),
    heritability = mean(h2),
    active_count = mean(rowSums(fit$n_comp_hist[keep, -1L, drop = FALSE])),
    component_count = colMeans(fit$n_comp_hist[keep, , drop = FALSE]),
    phenotype_scale = vy)
}

study06_gctb_block_compare <- function(sblr_fit, official, truth, fit_id) {
  causal <- truth$marker_truth$true_nonnull
  auc_s <- study06_gctb_auc(sblr_fit$pip, causal)
  auc_o <- study06_gctb_auc(official$pip, causal)
  top <- function(x, k) order(x, decreasing = TRUE, method = "radix")[seq_len(k)]
  metrics <- c(
    pip_pearson = stats::cor(sblr_fit$pip, official$pip),
    pip_spearman = stats::cor(sblr_fit$pip, official$pip,
      method = "spearman"),
    effect_pearson = stats::cor(sblr_fit$effect, official$effect),
    effect_rmse = sqrt(mean((sblr_fit$effect - official$effect)^2)),
    validation_g_pearson = stats::cor(sblr_fit$prediction,
      official$prediction),
    top50_overlap = length(intersect(top(sblr_fit$pip, 50L),
      top(official$pip, 50L))),
    top100_overlap = length(intersect(top(sblr_fit$pip, 100L),
      top(official$pip, 100L))),
    sblr_auprc = auc_s[["auprc"]], official_auprc = auc_o[["auprc"]],
    sblr_auroc = auc_s[["auroc"]], official_auroc = auc_o[["auroc"]],
    sblr_active_count = sblr_fit$active_count,
    official_active_count = official$active_count,
    sblr_mean_block_ve = sblr_fit$mean_block_ve,
    official_mean_block_ve = official$mean_block_ve,
    mean_block_ve_abs_diff = abs(sblr_fit$mean_block_ve -
      official$mean_block_ve),
    sblr_heritability = sblr_fit$heritability,
    official_heritability = official$heritability,
    heritability_abs_diff = abs(sblr_fit$heritability -
      official$heritability))
  data.frame(fit_id = fit_id, metric = names(metrics),
    value = unname(metrics), stringsAsFactors = FALSE)
}

study06_gctb_block_gate <- function(metrics, cfg) {
  pass <- vapply(names(cfg$gates), function(id) {
    x <- stats::setNames(metrics$value[metrics$fit_id == id],
      metrics$metric[metrics$fit_id == id])
    g <- cfg$gates[[id]]
    x[["pip_pearson"]] >= g[["pip_pearson"]] &&
      x[["effect_pearson"]] >= g[["effect_pearson"]] &&
      x[["validation_g_pearson"]] >= g[["validation_g_pearson"]] &&
      x[["mean_block_ve_abs_diff"]] <= g[["mean_block_ve_abs_diff"]] &&
      x[["heritability_abs_diff"]] <= g[["heritability_abs_diff"]]
  }, logical(1))
  list(pass = pass, decision = if (all(pass)) "SMALL-G1" else {
    snp <- all(vapply(names(pass), function(id) {
      x <- stats::setNames(metrics$value[metrics$fit_id == id],
        metrics$metric[metrics$fit_id == id])
      g <- cfg$gates[[id]]
      x[["pip_pearson"]] >= g[["pip_pearson"]] &&
        x[["effect_pearson"]] >= g[["effect_pearson"]] &&
        x[["validation_g_pearson"]] >= g[["validation_g_pearson"]]
    }, logical(1)))
    residual <- all(vapply(names(pass), function(id) {
      x <- stats::setNames(metrics$value[metrics$fit_id == id],
        metrics$metric[metrics$fit_id == id])
      g <- cfg$gates[[id]]
      x[["mean_block_ve_abs_diff"]] <= g[["mean_block_ve_abs_diff"]] &&
        x[["heritability_abs_diff"]] <= g[["heritability_abs_diff"]]
    }, logical(1)))
    if (snp && !residual) "SMALL-G2" else if (!snp && residual)
      "SMALL-G3" else "SMALL-G4"
  })
}

study06_residual_semantic_crosswalk <- function() data.frame(
  route = c("historical block-eigen", "new block-eigen", "BED", "CSR"),
  policy = c("global_projected_legacy", "gctb_block",
    "global individual-level residual", "global operator residual"),
  ves = c("global projected residual variance",
    "retained mean block residual variance",
    "global individual-level residual variance",
    "global operator residual variance"),
  heritability = c("Vg/(Vg+projected Ve)", "sum(block Vg)/phenotype variance",
    "Vg/(Vg+Ve)", "Vg/(Vg+Ve) under the operator contract"),
  historical_relabeling_allowed = c(FALSE, NA, NA, NA),
  stringsAsFactors = FALSE)
