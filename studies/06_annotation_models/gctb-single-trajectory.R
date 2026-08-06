study06_gctb_single_constants <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  list(
    schema = "sblrbench-study06-gctb-single-trajectory-v1",
    profile = "gctb_single_trajectory",
    specification_hash = "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56",
    truth_hash = "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb",
    marker_order_hash = "135c3604e0c4395349475b8126e0957db265c4b348a2e339daa3e7ddf2316a29",
    gwas_hash = "1bd0abb220c4e7f9ca58ed11f2c2913e6e142868ba17820b988670cd01cd4610",
    annotation_hash = "b2442b5c074b1cd6fa0eb047fcedae6d3296a15d2b5f985c2ac288b534a3b156",
    ld_block_hash = "369df6bbed513da7913f960b9f79e35d96cdeac321cf3105b303c932d8e2c0a4",
    official_version = "0.2.6",
    official_sha = "b95d3fcbad8ff358290922a58fff893439296138",
    sblr_sha = "a165fb0635afcb8a712e8658175dfbb19896b3c3",
    export_dir = file.path(root, "results", "local", "06_annotation_models", "gctb_parity", "export"),
    official_lib = file.path(root, "results", "local", "06_annotation_models", "gctb_parity", "rlib"),
    official_source = file.path(root, "results", "local", "06_annotation_models", "gctb_parity", "official-source"),
    smoke_dir = file.path(root, "results", "local", "06_annotation_models", "gctb_parity", "smoke"),
    output_dir = file.path(root, "results", "local", "06_annotation_models", "gctb_single_trajectory"),
    reference_dir = file.path(root, "results", "local", "06_annotation_models", "v2_paired_power_isolation"),
    niter = 9000L,
    burn = 3000L,
    out_freq = 1L,
    requested_seeds = c(D0 = 711121L, D1 = 721121L, D2 = 731121L),
    matched_gamma = c(0, 0.01, 0.1, 1),
    matched_start_pi = c(0.88, 0.06, 0.036, 0.024),
    native_gamma = c(0, 0.001, 0.01, 0.1, 1),
    native_start_pi = c(0.990, 0.005, 0.003, 0.001, 0.001)
  )
}

study06_gctb_sha256_real <- function(path) {
  if (!requireNamespace("openssl", quietly = TRUE)) stop("openssl is required")
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  paste0(openssl::sha256(raw))
}

study06_gctb_hash_object <- function(x) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf, version = 3)
  study06_gctb_sha256_real(tf)
}

study06_gctb_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

study06_gctb_validate_export <- function(cfg, verify_files = TRUE) {
  manifest_path <- file.path(cfg$export_dir, "manifest.json")
  study06_gctb_assert(file.exists(manifest_path), "validated GCTB export manifest is missing")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  checks <- c(
    specification_hash = identical(manifest$specification_hash, cfg$specification_hash),
    truth_hash = identical(manifest$truth_hash, cfg$truth_hash),
    marker_order_hash = identical(manifest$marker_order_hash, cfg$marker_order_hash),
    gwas_hash = identical(manifest$gwas_hash, cfg$gwas_hash),
    annotation_hash = identical(manifest$annotation_hash, cfg$annotation_hash),
    ld_block_hash = identical(manifest$ld_block_hash, cfg$ld_block_hash),
    package_sha = identical(manifest$official_package_sha, cfg$official_sha),
    sample_size = identical(as.integer(manifest$n), 1400L),
    markers = identical(as.integer(manifest$marker_count), 1500L),
    blocks = identical(as.integer(manifest$block_count), 15L),
    intercept_slot = isTRUE(manifest$annotations_include_explicit_official_intercept_slot)
  )
  study06_gctb_assert(all(checks), paste("export identity failure:", paste(names(checks)[!checks], collapse = ", ")))
  if (verify_files) {
    file_rows <- as.data.frame(manifest$files, stringsAsFactors = FALSE)
    got <- vapply(file_rows$path, function(p) study06_gctb_sha256_real(file.path(cfg$export_dir, p)), character(1))
    study06_gctb_assert(identical(unname(got), as.character(file_rows$sha256)), "one or more validated export file hashes changed")
  }
  truth <- readRDS(file.path(cfg$export_dir, "truth.rds"))
  ma <- data.table::fread(file.path(cfg$export_dir, "study06_informative.ma"), data.table = FALSE)
  annot <- data.table::fread(file.path(cfg$export_dir, "study06_informative.annot"), data.table = FALSE)
  study06_gctb_assert(nrow(ma) == 1500L && all(ma$N == 1400L), "GWAS N or marker count changed")
  study06_gctb_assert(identical(as.character(ma$SNP), truth$marker_ids), "GWAS marker order changed")
  study06_gctb_assert(identical(as.character(annot$SNP), truth$marker_ids), "annotation marker order changed")
  study06_gctb_assert(ncol(annot) == 5L && all(annot[[2L]] == 1), "official intercept slot or annotation columns changed")
  audit <- data.table::fread(file.path(cfg$export_dir, "ld_audit.csv"), data.table = FALSE)
  study06_gctb_assert(nrow(audit) == 15L && all(audit$marker_count == 100L), "LD block identity changed")
  study06_gctb_assert(all(audit$rank == 100L & audit$reader_k == audit$rank & audit$retained_mass == 1), "matched export no longer retains all positive modes")
  list(manifest = manifest, truth = truth, ma = ma, annotation = annot, ld_audit = audit, checks = checks)
}

study06_gctb_validate_official <- function(cfg) {
  .libPaths(c(cfg$official_lib, .libPaths()))
  study06_gctb_assert(requireNamespace("SBayesRC", quietly = TRUE), "official SBayesRC is unavailable in the isolated library")
  desc <- utils::packageDescription("SBayesRC")
  sha_field <- desc$RemoteSha %||% desc$GithubSHA1 %||% NA_character_
  path <- normalizePath(find.package("SBayesRC"), winslash = "/", mustWork = TRUE)
  study06_gctb_assert(as.character(utils::packageVersion("SBayesRC")) == cfg$official_version, "official SBayesRC version changed")
  study06_gctb_assert(startsWith(tolower(path), tolower(normalizePath(cfg$official_lib, winslash = "/", mustWork = TRUE))), "official SBayesRC loaded outside the isolated library")
  # v0.2.6 does not retain RemoteSha in its installed DESCRIPTION. The exact
  # SHA is therefore established by the clean pinned source checkout plus the
  # existing validated export/smoke provenance checked by the runner.
  if (!is.na(sha_field)) study06_gctb_assert(identical(sha_field, cfg$official_sha), "official SBayesRC RemoteSha changed")
  list(version = cfg$official_version, sha = cfg$official_sha, sha_field = sha_field,
    sha_provenance = if (is.na(sha_field)) "clean pinned source checkout and committed export/smoke records" else "installed DESCRIPTION RemoteSha",
    path = path, description = desc)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[[1L]])) y else x

study06_gctb_condition <- function(id, cfg) {
  switch(id,
    D0 = list(id = id, label = "matched official SBayesR", gamma = cfg$matched_gamma,
      start_pi = cfg$matched_start_pi, annot = "", b_tune = FALSE, thresh = 1,
      tune_iter = 150L, tune_burn = 100L, tune_step = c(0.995, 0.99, 0.95, 0.9), b_tune_prior = FALSE),
    D1 = list(id = id, label = "matched official SBayesRC", gamma = cfg$matched_gamma,
      start_pi = cfg$matched_start_pi, annot = file.path(cfg$export_dir, "study06_informative.annot"), b_tune = FALSE, thresh = 1,
      tune_iter = 150L, tune_burn = 100L, tune_step = c(0.995, 0.99, 0.95, 0.9), b_tune_prior = FALSE),
    D2 = list(id = id, label = "native official SBayesRC", gamma = cfg$native_gamma,
      start_pi = cfg$native_start_pi, annot = file.path(cfg$export_dir, "study06_informative.annot"), b_tune = TRUE, thresh = 0.995,
      tune_iter = 150L, tune_burn = 100L, tune_step = c(0.995, 0.99, 0.95, 0.9), b_tune_prior = FALSE),
    stop("unknown condition: ", id)
  )
}

study06_gctb_run_condition <- function(id, cfg) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLAS_NUM_THREADS = "1")
  export <- study06_gctb_validate_export(cfg)
  official <- study06_gctb_validate_official(cfg)
  cond <- study06_gctb_condition(id, cfg)
  out_dir <- file.path(cfg$output_dir, "runs", id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(out_dir, id)
  record_path <- file.path(out_dir, "run-record.json")
  if (file.exists(record_path)) {
    old <- jsonlite::read_json(record_path, simplifyVector = TRUE)
    study06_gctb_assert(identical(old$status, "complete"), paste(id, "has an incomplete existing record"))
    study06_gctb_assert(identical(old$config_hash, study06_gctb_hash_object(cond)), paste(id, "existing configuration differs"))
    message(id, " already complete with the matching identity; reusing without rerun")
    return(invisible(old))
  }
  config_hash <- study06_gctb_hash_object(cond)
  required <- c(paste0(prefix, ".rds"), paste0(prefix, ".txt"), paste0(prefix, ".par"), paste0(prefix, ".log"))
  existing <- list.files(out_dir, all.files = TRUE, no.. = TRUE)
  recover_complete <- length(existing) > 0L && all(file.exists(required))
  study06_gctb_assert(!length(existing) || recover_complete, paste(id, "output directory contains incomplete unregistered files"))
  warnings <- character()
  if (!recover_complete) {
    start <- Sys.time()
    console <- capture.output(
      withCallingHandlers(
        SBayesRC::sbayesrc(
        mafile = file.path(cfg$export_dir, "study06_informative.ma"),
        LDdir = file.path(cfg$export_dir, "ld"),
        outPrefix = prefix,
        annot = cond$annot,
        niter = cfg$niter,
        burn = cfg$burn,
        gamma = cond$gamma,
        startPi = cond$start_pi,
        starth2 = 0.5,
        thresh = cond$thresh,
        bTune = cond$b_tune,
        tuneIter = cond$tune_iter,
        tuneBurn = cond$tune_burn,
        tuneStep = cond$tune_step,
        bTunePrior = cond$b_tune_prior,
        sSamVe = "allMixVe",
        twopq = "nbsq",
        annoSigmaScale = 1,
        seed = unname(cfg$requested_seeds[[id]]),
        bOutDetail = TRUE,
        outFreq = cfg$out_freq,
        bOutBeta = TRUE,
        log2file = TRUE
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ), type = "output"
    )
    end <- Sys.time()
    writeLines(console, file.path(out_dir, "wrapper-console.txt"), useBytes = TRUE)
  } else {
    start <- file.info(required)$ctime[1L]
    end <- max(file.info(required)$mtime)
    warnings <- "runner recovered complete native outputs after post-fit manifest serialization failure; fit was not rerun"
  }
  study06_gctb_assert(all(file.exists(required)), paste(id, "did not produce required official outputs"))
  fit <- readRDS(paste0(prefix, ".rds"))
  tab <- data.table::fread(paste0(prefix, ".txt"), data.table = FALSE)
  study06_gctb_assert(identical(as.character(tab$SNP), export$truth$marker_ids), paste(id, "marker order mismatch"))
  study06_gctb_assert(identical(as.character(tab$A1), as.character(export$ma$A1)), paste(id, "effect allele mismatch"))
  study06_gctb_assert(all(is.finite(tab$BETA)) && all(is.finite(tab$PIP)) && all(tab$PIP >= 0 & tab$PIP <= 1), paste(id, "non-finite or invalid SNP output"))
  all_files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  hashes <- data.frame(
    path = substring(normalizePath(all_files, winslash = "/"), nchar(normalizePath(out_dir, winslash = "/")) + 2L),
    bytes = file.info(all_files)$size,
    sha256 = vapply(all_files, study06_gctb_sha256_real, character(1)),
    stringsAsFactors = FALSE
  )
  record <- list(
    schema = cfg$schema,
    status = "complete",
    condition = id,
    label = cond$label,
    specification_hash = cfg$specification_hash,
    truth_hash = cfg$truth_hash,
    config_hash = config_hash,
    official_version = official$version,
    official_sha = official$sha,
    official_path = official$path,
    requested_seed = unname(cfg$requested_seeds[[id]]),
    native_seed_contract = "public seed argument is not connected to native RNG streams; deterministic default-stream trajectory",
    niter = cfg$niter,
    burn = cfg$burn,
    retained_draws = cfg$niter - cfg$burn,
    start_utc = format(start, tz = "UTC", usetz = TRUE),
    end_utc = format(end, tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(end, start, units = "secs")),
    warnings = unique(warnings),
    recovered_without_rerun = recover_complete,
    files = hashes
  )
  jsonlite::write_json(record, record_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(record)
}

study06_gctb_read_trace_file <- function(path, burn = 3000L, niter = 9000L) {
  x <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (nrow(x) == niter - burn) return(x)
  study06_gctb_assert(nrow(x) == niter, paste("unexpected annotation trace length:", path))
  x[(burn + 1L):niter, , drop = FALSE]
}

study06_gctb_official_result <- function(id, cfg) {
  prefix <- file.path(cfg$output_dir, "runs", id, id)
  fit <- readRDS(paste0(prefix, ".rds"))
  snp <- data.table::fread(paste0(prefix, ".txt"), data.table = FALSE)
  cond <- study06_gctb_condition(id, cfg)
  keep <- (cfg$burn + 1L):cfg$niter
  histories <- list(
    hsq = fit$hsq_hist[keep],
    effect_variance = fit$ssq_hist[keep],
    pi = fit$pi_hist[keep, , drop = FALSE],
    component_count = fit$n_comp_hist[keep, , drop = FALSE],
    component_vg_fraction = fit$vg_comp_hist[keep, , drop = FALSE],
    active_count = rowSums(fit$n_comp_hist[keep, -1L, drop = FALSE]),
    block_residual_variance = fit$vare_hist,
    block_effect_variance = fit$ssq_block_hist,
    block_genetic_variance = fit$hsq_block_hist
  )
  histories$total_genetic_variance <- rowSums(histories$block_genetic_variance)
  histories$residual_variance <- rowMeans(histories$block_residual_variance)
  alpha <- conditional <- joint <- enrichment <- list()
  if (id != "D0") {
    trace_dir <- paste0(prefix, ".mcmcsamples")
    n_sticks <- length(cond$gamma) - 1L
    for (s in seq_len(n_sticks)) {
      alpha[[paste0("p", s)]] <- study06_gctb_read_trace_file(file.path(trace_dir, paste0(id, ".mcmcsamples.AnnoEffects_p", s)), cfg$burn)
      conditional[[paste0("p", s)]] <- study06_gctb_read_trace_file(file.path(trace_dir, paste0(id, ".mcmcsamples.AnnoCondProb_p", s)), cfg$burn)
    }
    for (k in seq_along(cond$gamma)) {
      joint[[paste0("pi", k)]] <- study06_gctb_read_trace_file(file.path(trace_dir, paste0(id, ".mcmcsamples.AnnoJointProb_pi", k)), cfg$burn)
    }
    enrichment <- study06_gctb_read_trace_file(file.path(trace_dir, paste0(id, ".mcmcsamples.AnnoPerSnpHsqEnrichment")), cfg$burn)
  }
  list(id = id, condition = cond, fit = fit, snp = snp, histories = histories,
    alpha = alpha, conditional = conditional, joint = joint, enrichment = enrichment,
    record = jsonlite::read_json(file.path(cfg$output_dir, "runs", id, "run-record.json"), simplifyVector = TRUE))
}

study06_gctb_drift_one <- function(x, condition, quantity, group) {
  x <- as.numeric(x)
  n <- length(x)
  study06_gctb_assert(n == 6000L && all(is.finite(x)), paste("invalid retained history:", condition, quantity))
  thirds <- split(seq_len(n), rep(1:3, each = n / 3L))
  first <- mean(x[thirds[[1L]]]); middle <- mean(x[thirds[[2L]]]); final <- mean(x[thirds[[3L]]])
  full_mean <- mean(x); full_sd <- stats::sd(x)
  drift <- final - first
  z <- if (is.finite(full_sd) && full_sd > 0) drift / full_sd else if (drift == 0) 0 else sign(drift) * Inf
  rel <- abs(drift) / max(abs(full_mean), .Machine$double.eps)
  flag <- if (abs(z) < 0.25) "little visible drift" else if (abs(z) < 0.75) "moderate visible drift" else "substantial visible drift"
  cm <- cumsum(x) / seq_along(x)
  time <- seq(0, 1, length.out = n)
  ess <- tryCatch(as.numeric(posterior::ess_bulk(matrix(x, ncol = 1))), error = function(e) NA_real_)
  ac <- stats::acf(x, lag.max = 50, plot = FALSE, na.action = na.pass)$acf
  data.frame(condition = condition, group = group, quantity = quantity,
    mean = full_mean, sd = full_sd, first_third = first, middle_third = middle,
    final_third = final, first_half = mean(x[seq_len(n / 2L)]),
    second_half = mean(x[(n / 2L + 1L):n]), absolute_drift = drift,
    relative_drift = rel, standardized_drift = z,
    linear_trend = unname(stats::coef(stats::lm(x ~ time))[2L]),
    max_cumulative_mean_displacement = max(abs(cm - full_mean)),
    max_cumulative_mean_displacement_sd = if (full_sd > 0) max(abs(cm - full_mean)) / full_sd else 0,
    acf_lag1 = ac[2L], acf_lag10 = ac[11L], acf_lag50 = ac[51L],
    single_chain_bulk_ess = ess, flag = flag, stringsAsFactors = FALSE)
}

study06_gctb_drift_table <- function(off) {
  rows <- list(); add <- function(x, name, group) rows[[length(rows) + 1L]] <<- study06_gctb_drift_one(x, off$id, name, group)
  add(off$histories$hsq, "heritability", "variance")
  add(off$histories$residual_variance, "residual_variance", "variance")
  add(off$histories$total_genetic_variance, "total_genetic_variance", "variance")
  add(off$histories$effect_variance, "effect_variance", "variance")
  add(off$histories$active_count, "active_count", "occupancy")
  for (j in seq_len(ncol(off$histories$pi))) add(off$histories$pi[, j], paste0("pi", j), "mixture")
  for (j in seq_len(ncol(off$histories$component_count))) add(off$histories$component_count[, j], paste0("component_count_", j), "occupancy")
  for (j in seq_len(ncol(off$histories$component_vg_fraction))) add(off$histories$component_vg_fraction[, j], paste0("component_vg_fraction_", j), "architecture")
  for (s in names(off$alpha)) for (j in seq_len(ncol(off$alpha[[s]]))) add(off$alpha[[s]][, j], paste0("alpha_", s, "_", names(off$alpha[[s]])[j]), "alpha")
  for (s in names(off$conditional)) for (j in seq_len(ncol(off$conditional[[s]]))) add(off$conditional[[s]][, j], paste0("conditional_", s, "_", names(off$conditional[[s]])[j]), "annotation_probability")
  for (k in names(off$joint)) for (j in seq_len(ncol(off$joint[[k]]))) add(off$joint[[k]][, j], paste0("joint_", k, "_", names(off$joint[[k]])[j]), "annotation_probability")
  do.call(rbind, rows)
}

study06_gctb_auc <- function(score, truth) {
  truth <- as.logical(truth); n1 <- sum(truth); n0 <- sum(!truth)
  ranks <- rank(score, ties.method = "average")
  auroc <- (sum(ranks[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  ord <- order(score, decreasing = TRUE, method = "radix")
  tp <- cumsum(truth[ord]); fp <- cumsum(!truth[ord])
  recall <- tp / n1; precision <- tp / pmax(tp + fp, 1)
  auprc <- sum(diff(c(0, recall)) * precision)
  c(auprc = auprc, auroc = auroc)
}

study06_gctb_bfdr <- function(pip, q) {
  ord <- order(pip, decreasing = TRUE, method = "radix")
  efdr <- cumsum(1 - pip[ord]) / seq_along(ord)
  k <- max(c(0L, which(efdr <= q)))
  if (k == 0L) integer() else ord[seq_len(k)]
}

study06_gctb_detection_metrics <- function(method, pip, effect, truth) {
  causal <- as.logical(truth$marker_truth$true_nonnull)
  true_effect <- truth$effects
  ord <- order(pip, decreasing = TRUE, method = "radix")
  ranks <- rank(-pip, ties.method = "average")
  auc <- study06_gctb_auc(pip, causal)
  ks <- c(10L, 25L, 50L, 100L, 200L)
  out <- data.frame(method = method, metric = names(auc), value = unname(auc), stringsAsFactors = FALSE)
  for (k in ks) {
    idx <- ord[seq_len(k)]
    out <- rbind(out,
      data.frame(method = method, metric = paste0("precision_at_", k), value = mean(causal[idx])),
      data.frame(method = method, metric = paste0("recall_at_", k), value = sum(causal[idx]) / sum(causal)))
  }
  f5 <- study06_gctb_bfdr(pip, 0.05); f10 <- study06_gctb_bfdr(pip, 0.10)
  extras <- c(
    causal_rank_median = stats::median(ranks[causal]), causal_rank_mean = mean(ranks[causal]),
    mean_causal_pip = mean(pip[causal]), mean_noncausal_pip = mean(pip[!causal]),
    bfdr05_discoveries = length(f5), bfdr05_true = sum(causal[f5]), bfdr05_false = sum(!causal[f5]),
    bfdr10_discoveries = length(f10), bfdr10_true = sum(causal[f10]), bfdr10_false = sum(!causal[f10]),
    effect_truth_correlation = stats::cor(effect, true_effect),
    causal_effect_truth_correlation = stats::cor(effect[causal], true_effect[causal]),
    causal_effect_rmse = sqrt(mean((effect[causal] - true_effect[causal])^2)),
    noncausal_effect_scale = sqrt(mean(effect[!causal]^2))
  )
  rbind(out, data.frame(method = method, metric = names(extras), value = unname(extras), stringsAsFactors = FALSE))
}

study06_gctb_reference <- function(fit_id, cfg) {
  path <- file.path(cfg$reference_dir, "checkpoints", paste0(fit_id, ".rds"))
  study06_gctb_assert(file.exists(path), paste("missing sblr reference", fit_id))
  x <- readRDS(path)
  study06_gctb_assert(identical(x$identity$spec_hash, cfg$specification_hash) && identical(x$identity$truth_hash, cfg$truth_hash), paste("reference identity changed", fit_id))
  r <- x$result
  list(id = fit_id, checkpoint = x, result = r,
    pip = as.numeric(r$estimates$pip), effect = as.numeric(r$estimates$effects),
    chains = r$native_fit$chains)
}

study06_gctb_prediction_metrics <- function(method, effect, truth) {
  pred <- drop(truth$validation_x %*% effect)
  data.frame(method = method, metric = c("validation_genetic_correlation", "phenotype_prediction_correlation", "prediction_slope", "prediction_rmse", "predicted_genetic_variance"),
    value = c(stats::cor(pred, truth$validation_genetic), stats::cor(pred, truth$validation_phenotype),
      unname(stats::coef(stats::lm(truth$validation_genetic ~ pred))[2L]),
      sqrt(mean((pred - truth$validation_genetic)^2)), stats::var(pred)), stringsAsFactors = FALSE)
}

study06_gctb_pair_metrics <- function(comparison, a_name, a_pip, a_effect, b_name, b_pip, b_effect, truth) {
  reg_p <- stats::coef(stats::lm(a_pip ~ b_pip)); reg_e <- stats::coef(stats::lm(a_effect ~ b_effect))
  causal <- as.logical(truth$marker_truth$true_nonnull)
  sign_idx <- abs(a_effect) > 0 | abs(b_effect) > 0
  vals <- c(
    pip_pearson = stats::cor(a_pip, b_pip), pip_spearman = stats::cor(a_pip, b_pip, method = "spearman"),
    pip_regression_intercept = reg_p[1], pip_regression_slope = reg_p[2],
    pip_mad = mean(abs(a_pip - b_pip)), pip_rmse = sqrt(mean((a_pip - b_pip)^2)),
    pip_median_abs_difference = stats::median(abs(a_pip - b_pip)), pip_max_abs_difference = max(abs(a_pip - b_pip)),
    effect_pearson = stats::cor(a_effect, b_effect), effect_spearman = stats::cor(a_effect, b_effect, method = "spearman"),
    effect_regression_intercept = reg_e[1], effect_regression_slope = reg_e[2],
    effect_mae = mean(abs(a_effect - b_effect)), effect_rmse = sqrt(mean((a_effect - b_effect)^2)),
    effect_sign_agreement = mean(sign(a_effect[sign_idx]) == sign(b_effect[sign_idx])),
    causal_effect_correlation = stats::cor(a_effect[causal], b_effect[causal]),
    causal_effect_rmse = sqrt(mean((a_effect[causal] - b_effect[causal])^2)),
    genetic_value_correlation = stats::cor(drop(truth$validation_x %*% a_effect), drop(truth$validation_x %*% b_effect)),
    pip_kendall = stats::cor(a_pip, b_pip, method = "kendall")
  )
  data.frame(comparison = comparison, method_a = a_name, method_b = b_name,
    metric = names(vals), value = unname(vals), stringsAsFactors = FALSE)
}

study06_gctb_top_overlap <- function(comparison, a_name, a_pip, b_name, b_pip, truth) {
  causal <- as.logical(truth$marker_truth$true_nonnull); ks <- c(10L, 25L, 50L, 100L, 200L)
  do.call(rbind, lapply(ks, function(k) {
    a <- order(a_pip, decreasing = TRUE, method = "radix")[seq_len(k)]
    b <- order(b_pip, decreasing = TRUE, method = "radix")[seq_len(k)]
    inter <- intersect(a, b)
    data.frame(comparison = comparison, method_a = a_name, method_b = b_name, k = k,
      overlap = length(inter), jaccard = length(inter) / length(union(a, b)),
      causal_a = sum(causal[a]), causal_b = sum(causal[b]), shared_causal = sum(causal[inter]),
      precision_a = mean(causal[a]), precision_b = mean(causal[b]),
      recall_a = sum(causal[a]) / sum(causal), recall_b = sum(causal[b]) / sum(causal))
  }))
}

study06_gctb_pip_bins <- function(comparison, a_name, a_pip, b_name, b_pip, truth) {
  ord <- order(a_pip, method = "radix")
  bin <- integer(length(a_pip)); bin[ord] <- rep(seq_len(10L), each = length(a_pip) / 10L)
  causal <- as.logical(truth$marker_truth$true_nonnull)
  aggregate(cbind(pip_a = a_pip, pip_b = b_pip, causal = as.numeric(causal)),
    list(comparison = rep(comparison, length(bin)), method_a = rep(a_name, length(bin)),
      method_b = rep(b_name, length(bin)), bin = bin), mean)
}

study06_gctb_official_prior <- function(off, annotation) {
  if (!length(off$alpha)) return(NULL)
  A <- as.matrix(annotation)
  study06_gctb_assert(ncol(A) == 4L && all(A[, 1L] == 1), "annotation matrix does not contain the registered intercept")
  sticks <- lapply(off$alpha, as.matrix)
  k <- length(sticks) + 1L
  acc <- matrix(0, nrow(A), k)
  chunk <- 250L
  for (lo in seq(1L, nrow(sticks[[1L]]), by = chunk)) {
    hi <- min(nrow(sticks[[1L]]), lo + chunk - 1L)
    q1 <- stats::pnorm(A %*% t(sticks[[1L]][lo:hi, , drop = FALSE]))
    acc[, 1L] <- acc[, 1L] + rowSums(1 - q1)
    remaining <- q1
    if (length(sticks) > 1L) for (s in 2:length(sticks)) {
      q <- stats::pnorm(A %*% t(sticks[[s]][lo:hi, , drop = FALSE]))
      acc[, s] <- acc[, s] + rowSums(remaining * (1 - q))
      remaining <- remaining * q
    }
    acc[, k] <- acc[, k] + rowSums(remaining)
  }
  acc / nrow(sticks[[1L]])
}

study06_gctb_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
