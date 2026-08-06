#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
has_flag <- function(flag) flag %in% args

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
helper <- file.path(root, "studies", "06_annotation_models", "gctb-single-trajectory.R")
if (!file.exists(helper)) stop("run from the sblrbench repository root")
source(helper, local = .GlobalEnv)
cfg <- study06_gctb_single_constants(root)

condition <- arg_value("--condition")
validate_only <- identical(tolower(arg_value("--validate-only", "false")), "true")
analyze_only <- identical(tolower(arg_value("--analyze-only", "false")), "true")

preflight <- function() {
  export <- study06_gctb_validate_export(cfg)
  official <- study06_gctb_validate_official(cfg)
  source_head <- system2("git", c("-C", shQuote(cfg$official_source), "rev-parse", "HEAD"), stdout = TRUE)
  source_status <- system2("git", c("-C", shQuote(cfg$official_source), "status", "--short"), stdout = TRUE)
  sibling_head <- system2("git", c("-C", "../sblr", "rev-parse", "HEAD"), stdout = TRUE)
  sibling_status <- system2("git", c("-C", "../sblr", "status", "--short"), stdout = TRUE)
  study06_gctb_assert(identical(source_head, cfg$official_sha) && !length(source_status), "official source checkout changed")
  study06_gctb_assert(identical(sibling_head, cfg$sblr_sha) && !length(sibling_status), "sibling sblr identity changed")
  list(export = export, official = official)
}

fixed_stream_check <- function() {
  long <- readRDS(file.path(cfg$output_dir, "runs", "D0", "D0.rds"))
  smoke <- readRDS(file.path(cfg$smoke_dir, "G0--chain1", "G0--chain1.rds"))
  fields <- c("hsq_hist", "ssq_hist", "pi_hist", "vg_comp_hist", "n_comp_hist")
  checks <- lapply(fields, function(n) {
    a <- long[[n]]; b <- smoke[[n]]
    if (is.null(dim(b))) a <- a[seq_along(b)] else a <- a[seq_len(nrow(b)), , drop = FALSE]
    data.frame(field = n, comparable_values = length(b), exact = identical(a, b),
      max_absolute_difference = max(abs(as.numeric(a) - as.numeric(b))))
  })
  out <- do.call(rbind, checks)
  study06_gctb_write_csv(out, file.path(cfg$output_dir, "analysis", "fixed_stream_reproducibility.csv"))
  study06_gctb_assert(all(out$exact), "official D0 default stream is not reproducible against the previous fresh-process smoke")
  out
}

official_summary <- function(off) {
  h <- off$histories
  data.frame(method = off$id,
    quantity = c("active_count", "heritability", "total_genetic_variance", "residual_variance", "effect_variance",
      paste0("component_count_", seq_len(ncol(h$component_count))), paste0("component_pi_", seq_len(ncol(h$pi))),
      paste0("component_vg_fraction_", seq_len(ncol(h$component_vg_fraction)))),
    mean = c(mean(h$active_count), mean(h$hsq), mean(h$total_genetic_variance), mean(h$residual_variance), mean(h$effect_variance),
      colMeans(h$component_count), colMeans(h$pi), colMeans(h$component_vg_fraction)),
    sd = c(sd(h$active_count), sd(h$hsq), sd(h$total_genetic_variance), sd(h$residual_variance), sd(h$effect_variance),
      apply(h$component_count, 2, sd), apply(h$pi, 2, sd), apply(h$component_vg_fraction, 2, sd)))
}

extract_tuning <- function() {
  path <- file.path(cfg$output_dir, "runs", "D2", "D2.mcmcsamples", "D2_tune.txt")
  study06_gctb_assert(file.exists(path), "D2 tuning table is missing")
  x <- data.table::fread(path, data.table = FALSE)
  candidate <- x[is.finite(x$r) & is.finite(x$rel_r), , drop = FALSE]
  best <- candidate[which.max(candidate$rel_r), , drop = FALSE]
  selected <- if (best$rel_r > 1.25 || candidate$r[candidate$thresh == max(candidate$thresh)][1L] < 0) best$thresh else max(cfg$native_start_pi * 0 + c(0.995, 0.99, 0.95, 0.9, 0)[1:5])
  # The official wrapper's ordinary branch selects max(tuneStep), i.e. 0.995.
  if (!(best$rel_r > 1.25 || candidate$r[candidate$thresh == max(candidate$thresh)][1L] < 0)) selected <- 0.995
  x$selected <- x$thresh == selected
  study06_gctb_write_csv(x, file.path(cfg$output_dir, "analysis", "D2_tuning.csv"))
  list(table = x, selected = selected)
}

plot_comparison <- function(name, a, b, a_label, b_label, bins) {
  pdir <- file.path(cfg$output_dir, "analysis", "plots")
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file.path(pdir, paste0(name, "--pip-scatter.png")), 1200, 1000, res = 150)
  graphics::smoothScatter(a, b, xlab = paste(a_label, "PIP"), ylab = paste(b_label, "PIP"), main = name)
  graphics::abline(0, 1, col = "firebrick", lwd = 2)
  grDevices::dev.off()
  grDevices::png(file.path(pdir, paste0(name, "--pip-bins.png")), 1200, 850, res = 150)
  graphics::matplot(bins$bin, cbind(bins$pip_a, bins$pip_b), type = "b", pch = c(16, 17), lty = 1,
    xlab = paste("PIP quantile bin ordered by", a_label), ylab = "mean PIP", main = paste(name, "binned PIP comparison"))
  graphics::legend("topleft", c(a_label, b_label), pch = c(16, 17), col = 1:2, lty = 1)
  grDevices::dev.off()
}

analyze <- function() {
  pf <- preflight(); truth <- pf$export$truth
  out <- lapply(c("D0", "D1", "D2"), study06_gctb_official_result, cfg = cfg); names(out) <- c("D0", "D1", "D2")
  study06_gctb_assert(all(vapply(out, function(x) identical(as.character(x$snp$SNP), truth$marker_ids), logical(1))), "official SNP alignment failed")
  repro <- fixed_stream_check()
  analysis_dir <- file.path(cfg$output_dir, "analysis"); dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  drift <- do.call(rbind, lapply(out, study06_gctb_drift_table))
  study06_gctb_write_csv(drift, file.path(analysis_dir, "within_trajectory_drift.csv"))
  arch <- do.call(rbind, lapply(out, official_summary))
  study06_gctb_write_csv(arch, file.path(analysis_dir, "official_architecture.csv"))

  refs <- list(
    baseline_block = study06_gctb_reference("baseline--block_eigen", cfg),
    baseline_bed = study06_gctb_reference("baseline--bed", cfg),
    learned_block = study06_gctb_reference("learned_informative--block_eigen", cfg),
    learned_bed = study06_gctb_reference("learned_informative--bed", cfg),
    fixed_block = study06_gctb_reference("fixed_true_alpha--block_eigen", cfg),
    fixed_bed = study06_gctb_reference("fixed_true_alpha--bed", cfg),
    shuffled_block = study06_gctb_reference("learned_shuffled--block_eigen", cfg),
    shuffled_bed = study06_gctb_reference("learned_shuffled--bed", cfg)
  )
  methods <- c(lapply(out, function(x) list(id = x$id, pip = x$snp$PIP, effect = x$snp$BETA)), refs)
  detection <- do.call(rbind, lapply(methods, function(x) study06_gctb_detection_metrics(x$id, x$pip, x$effect, truth)))
  prediction <- do.call(rbind, lapply(methods, function(x) study06_gctb_prediction_metrics(x$id, x$effect, truth)))
  study06_gctb_write_csv(detection, file.path(analysis_dir, "causal_detection_metrics.csv"))
  study06_gctb_write_csv(prediction, file.path(analysis_dir, "prediction_metrics.csv"))

  comparisons <- list(
    C0 = list(out$D0, refs$baseline_block),
    C1 = list(out$D1, refs$learned_block),
    C2 = list(out$D1, refs$learned_bed),
    C3_block = list(out$D1, refs$fixed_block),
    C3_bed = list(out$D1, refs$fixed_bed),
    C4 = list(out$D2, out$D1),
    C5 = list(out$D1, out$D0)
  )
  pair_rows <- overlap_rows <- bin_rows <- list()
  bfdr_rows <- effect_overlap_rows <- list()
  for (nm in names(comparisons)) {
    a <- comparisons[[nm]][[1L]]; b <- comparisons[[nm]][[2L]]
    a_name <- a$id; b_name <- b$id
    a_pip <- if (!is.null(a$snp)) a$snp$PIP else a$pip; b_pip <- if (!is.null(b$snp)) b$snp$PIP else b$pip
    a_effect <- if (!is.null(a$snp)) a$snp$BETA else a$effect; b_effect <- if (!is.null(b$snp)) b$snp$BETA else b$effect
    pair_rows[[nm]] <- study06_gctb_pair_metrics(nm, a_name, a_pip, a_effect, b_name, b_pip, b_effect, truth)
    overlap_rows[[nm]] <- study06_gctb_top_overlap(nm, a_name, a_pip, b_name, b_pip, truth)
    bin_rows[[nm]] <- study06_gctb_pip_bins(nm, a_name, a_pip, b_name, b_pip, truth)
    plot_comparison(nm, a_pip, b_pip, a_name, b_name, bin_rows[[nm]])
    bfdr_rows[[nm]] <- do.call(rbind, lapply(c(0.05, 0.10), function(q) {
      ia <- study06_gctb_bfdr(a_pip, q); ib <- study06_gctb_bfdr(b_pip, q)
      causal <- truth$marker_truth$true_nonnull
      data.frame(comparison = nm, q = q, discoveries_a = length(ia), discoveries_b = length(ib),
        overlap = length(intersect(ia, ib)), true_a = sum(causal[ia]), true_b = sum(causal[ib]),
        false_a = sum(!causal[ia]), false_b = sum(!causal[ib]), shared_true = sum(causal[intersect(ia, ib)]))
    }))
    effect_overlap_rows[[nm]] <- do.call(rbind, lapply(c(10L, 25L, 50L, 100L, 200L), function(k) {
      ia <- order(abs(a_effect), decreasing = TRUE)[seq_len(k)]; ib <- order(abs(b_effect), decreasing = TRUE)[seq_len(k)]
      data.frame(comparison = nm, k = k, overlap = length(intersect(ia, ib)), jaccard = length(intersect(ia, ib)) / length(union(ia, ib)))
    }))
  }
  pair <- do.call(rbind, pair_rows); overlap <- do.call(rbind, overlap_rows); bins <- do.call(rbind, bin_rows)
  study06_gctb_write_csv(pair, file.path(analysis_dir, "cross_method_agreement.csv"))
  study06_gctb_write_csv(overlap, file.path(analysis_dir, "top_marker_overlap.csv"))
  study06_gctb_write_csv(bins, file.path(analysis_dir, "pip_calibration_bins.csv"))
  study06_gctb_write_csv(do.call(rbind, bfdr_rows), file.path(analysis_dir, "bayesian_fdr_overlap.csv"))
  study06_gctb_write_csv(do.call(rbind, effect_overlap_rows), file.path(analysis_dir, "top_effect_overlap.csv"))

  high_effect <- order(abs(truth$effects), decreasing = TRUE)[seq_len(25L)]
  high_effect_ranks <- do.call(rbind, lapply(methods, function(x) data.frame(method = x$id,
    marker_id = truth$marker_ids[high_effect], true_effect = truth$effects[high_effect],
    pip = x$pip[high_effect], pip_rank = rank(-x$pip, ties.method = "average")[high_effect],
    effect = x$effect[high_effect], effect_abs_rank = rank(-abs(x$effect), ties.method = "average")[high_effect])))
  study06_gctb_write_csv(high_effect_ranks, file.path(analysis_dir, "high_effect_causal_marker_ranks.csv"))

  per_chain <- list()
  for (ref_name in c("learned_block", "learned_bed")) {
    ref <- refs[[ref_name]]
    for (i in seq_along(ref$chains)) {
      ch <- ref$chains[[i]]
      tmp <- study06_gctb_pair_metrics(paste0("D1_vs_", ref_name, "_chain", i), "D1", out$D1$snp$PIP, out$D1$snp$BETA,
        paste0(ref$id, "--chain", i), as.numeric(ch$dm), as.numeric(ch$bm), truth)
      top <- study06_gctb_top_overlap(paste0("D1_vs_", ref_name, "_chain", i), "D1", out$D1$snp$PIP,
        paste0(ref$id, "--chain", i), as.numeric(ch$dm), truth)
      alpha <- ch$convergence_trace$alpha
      alpha_mean <- if (ncol(alpha)) colMeans(alpha) else numeric()
      alpha_names <- colnames(alpha) %||% paste0("alpha_", seq_along(alpha_mean))
      has_ncomp <- length(ch$ncomp) == 4L
      meta <- data.frame(reference = ref_name, chain = i,
        active_count = if (has_ncomp) 1500 - ch$ncomp[1L] else sum(ch$dm),
        active_count_source = if (has_ncomp) "posterior mean component occupancy" else "expected active count from chain PIPs; occupancy trace unavailable",
        heritability = mean(ch$vgs / (ch$vgs + ch$ves)), stringsAsFactors = FALSE)
      per_chain[[length(per_chain) + 1L]] <- list(metrics = tmp, top = top, meta = meta,
        alpha = data.frame(reference = ref_name, chain = i, coefficient = alpha_names, mean = unname(alpha_mean)))
    }
  }
  study06_gctb_write_csv(do.call(rbind, lapply(per_chain, `[[`, "metrics")), file.path(analysis_dir, "sblr_chain_agreement.csv"))
  study06_gctb_write_csv(do.call(rbind, lapply(per_chain, `[[`, "top")), file.path(analysis_dir, "sblr_chain_top_overlap.csv"))
  study06_gctb_write_csv(do.call(rbind, lapply(per_chain, `[[`, "meta")), file.path(analysis_dir, "sblr_chain_architecture.csv"))
  study06_gctb_write_csv(do.call(rbind, lapply(per_chain, `[[`, "alpha")), file.path(analysis_dir, "sblr_chain_alpha.csv"))

  alpha_summary <- do.call(rbind, lapply(c("D1", "D2"), function(id) do.call(rbind, lapply(names(out[[id]]$alpha), function(s) {
    x <- out[[id]]$alpha[[s]]
    data.frame(condition = id, stick = s, annotation = names(x), mean = colMeans(x), sd = apply(x, 2, sd), stringsAsFactors = FALSE)
  }))))
  study06_gctb_write_csv(alpha_summary, file.path(analysis_dir, "official_alpha_summary.csv"))
  probability_summary <- do.call(rbind, lapply(c("D1", "D2"), function(id) {
    items <- c(out[[id]]$conditional, out[[id]]$joint)
    do.call(rbind, lapply(names(items), function(n) data.frame(condition = id, probability = n,
      annotation = names(items[[n]]), mean = colMeans(items[[n]]), sd = apply(items[[n]], 2, sd))))
  }))
  study06_gctb_write_csv(probability_summary, file.path(analysis_dir, "official_annotation_probability_summary.csv"))
  prior <- lapply(c("D1", "D2"), function(id) study06_gctb_official_prior(out[[id]], truth$annotations)); names(prior) <- c("D1", "D2")
  saveRDS(prior, file.path(analysis_dir, "official_marker_prior_probabilities.rds"), version = 3)
  prior_summary <- do.call(rbind, lapply(names(prior), function(id) {
    p <- prior[[id]]; causal <- truth$marker_truth$true_nonnull; enriched <- truth$annotations[, "enriched_binary"] == 1
    data.frame(condition = id, quantity = c("expected_active", "nonnull_enriched", "nonnull_unenriched", "nonnull_causal", "nonnull_noncausal"),
      value = c(sum(1 - p[, 1L]), mean(1 - p[enriched, 1L]), mean(1 - p[!enriched, 1L]), mean(1 - p[causal, 1L]), mean(1 - p[!causal, 1L])))
  }))
  study06_gctb_write_csv(prior_summary, file.path(analysis_dir, "official_marker_prior_summary.csv"))

  tuning <- extract_tuning()
  smoke_tuning <- c(0.65984, 0.65773, 0.66519, 0.68057)
  tuning_reproduces_smoke <- isTRUE(all.equal(tuning$table$r, smoke_tuning, tolerance = 1e-5)) && identical(tuning$selected, 0.995)
  runtime <- do.call(rbind, lapply(out, function(x) data.frame(condition = x$id, elapsed_seconds = x$record$elapsed_seconds,
    warnings = paste(x$record$warnings, collapse = " | "))))
  study06_gctb_write_csv(runtime, file.path(analysis_dir, "runtime.csv"))

  all_files <- list.files(cfg$output_dir, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[basename(all_files) != "manifest.json"]
  package_files <- list.files(pf$official$path, recursive = TRUE, full.names = TRUE)
  package_hashes <- setNames(vapply(package_files, study06_gctb_sha256_real, character(1)),
    substring(normalizePath(package_files, winslash = "/"), nchar(normalizePath(pf$official$path, winslash = "/")) + 2L))
  manifest <- list(schema = cfg$schema, profile = cfg$profile, specification_hash = cfg$specification_hash,
    truth_hash = cfg$truth_hash, marker_order_hash = cfg$marker_order_hash, gwas_hash = cfg$gwas_hash,
    annotation_hash = cfg$annotation_hash, ld_block_hash = cfg$ld_block_hash,
    official_version = cfg$official_version, official_sha = cfg$official_sha,
    official_loaded_path = pf$official$path, installed_package_description_sha256 = study06_gctb_sha256_real(file.path(pf$official$path, "DESCRIPTION")),
    installed_package_tree_sha256 = study06_gctb_hash_object(package_hashes),
    fixed_stream_reproducible = all(repro$exact), D2_selected_threshold = tuning$selected,
    D2_tuning_reproduces_smoke = tuning_reproduces_smoke,
    configuration_hashes = setNames(lapply(c("D0", "D1", "D2"), function(id) study06_gctb_hash_object(study06_gctb_condition(id, cfg))), c("D0", "D1", "D2")),
    script_hashes = list(
      runner = study06_gctb_sha256_real(file.path(root, "scripts", "run_study06_gctb_single_trajectory.R")),
      helper = study06_gctb_sha256_real(file.path(root, "studies", "06_annotation_models", "gctb-single-trajectory.R"))),
    files = data.frame(path = substring(normalizePath(all_files, winslash = "/"), nchar(normalizePath(cfg$output_dir, winslash = "/")) + 2L),
      bytes = file.info(all_files)$size, sha256 = vapply(all_files, study06_gctb_sha256_real, character(1)), stringsAsFactors = FALSE))
  jsonlite::write_json(manifest, file.path(cfg$output_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  cat("Analysis complete\n")
  cat("D2 selected threshold:", tuning$selected, "\n")
  cat("Fixed stream reproducible:", all(repro$exact), "\n")
  invisible(list(drift = drift, detection = detection, pair = pair, overlap = overlap, architecture = arch, alpha = alpha_summary, manifest = manifest))
}

pf <- preflight()
if (validate_only) {
  cat("VALID: exact export, official package, source checkout, and sibling identity\n")
  quit(save = "no", status = 0L)
}
if (!is.null(condition)) {
  study06_gctb_run_condition(condition, cfg)
  quit(save = "no", status = 0L)
}
if (analyze_only) {
  analyze()
  quit(save = "no", status = 0L)
}

for (id in c("D0", "D1", "D2")) {
  status <- system2(file.path(R.home("bin"), "Rscript.exe"), c("scripts/run_study06_gctb_single_trajectory.R", "--condition", id),
    env = c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1", "BLAS_NUM_THREADS=1"))
  if (status != 0L) stop(id, " failed; later conditions were not launched")
  if (id == "D0") fixed_stream_check()
}
analyze()
