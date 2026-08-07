# Fast final aggregation for Study 06. Requires estimability-and-contrasts.R.

finalize_study06_estimability <- function(root = ".") {
  out <- file.path(root, "results/local/06_annotation_models/estimability_and_contrasts")
  td <- file.path(out, "tables")
  truth <- readRDS(file.path(root, "results/local/06_annotation_models/gctb_parity/export/truth.rds"))
  A <- truth$annotations
  fit_paths <- c(bed = file.path(root, "results/local/06_annotation_models/v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_bed_bayesrc.rds"),
    block = file.path(root, "results/local/06_annotation_models/v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_block_eigen_sbayesrc.rds"))
  fits <- lapply(fit_paths, function(p) readRDS(p)$result$native_fit)
  source(file.path(root, "studies/06_annotation_models/spec.R"), local = TRUE)
  source(file.path(root, "studies/06_annotation_models/annotation-design.R"), local = TRUE)
  truth_alpha <- construct_annotation_truth(A, spec)$informative_annotations
  truth_pi <- study06_pi(A %*% truth_alpha)
  marker <- readRDS(file.path(out, "marker_derived_chain_summaries.rds"))

  raw <- raw_chain <- expected <- representative <- ranking <- gv <- eta_chain <- probability_crosscheck <- list()
  reps <- study06_representatives(A)
  for (route in names(fits)) {
    alpha <- study06_chain_alpha(fits[[route]])
    rr <- study06_raw_alpha(route, alpha, truth_alpha)
    raw[[route]] <- rr$summary; raw_chain[[route]] <- rr$chains
    q <- alpha$descriptor; active_draws <- vector("list", 4L)
    rep_draws <- lapply(seq_len(4L), function(i) array(NA_real_, c(nrow(alpha$traces[[i]]), length(reps), 5L),
      dimnames = list(NULL, NULL, c("active", paste0("pi", 0:3)))))
    for (ch in seq_len(4L)) {
      tr <- alpha$traces[[ch]]
      idx1 <- which(q$stick_index == 1L)
      active_draws[[ch]] <- rowSums(stats::pnorm(A %*% t(tr[, idx1, drop = FALSE])))
      eta <- lapply(seq_len(3L), function(s) tr[, which(q$stick_index == s), drop = FALSE] %*% t(A[reps, , drop = FALSE]))
      c1 <- stats::pnorm(eta[[1]]); c2 <- stats::pnorm(eta[[2]]); c3 <- stats::pnorm(eta[[3]])
      rep_draws[[ch]][, , 2L] <- 1 - c1
      rep_draws[[ch]][, , 3L] <- c1 * (1 - c2)
      rep_draws[[ch]][, , 4L] <- c1 * c2 * (1 - c3)
      rep_draws[[ch]][, , 5L] <- c1 * c2 * c3
      rep_draws[[ch]][, , 1L] <- c1
    }
    dg <- study06_diag(active_draws); pooled <- unlist(active_draws)
    expected[[route]] <- data.frame(route, posterior_mean = mean(pooled), posterior_sd = stats::sd(pooled),
      q025 = stats::quantile(pooled, .025), q975 = stats::quantile(pooled, .975),
      truth = sum(1 - truth_pi[, 1]), bias = mean(pooled) - sum(1 - truth_pi[, 1]), t(dg))
    for (r in seq_along(reps)) for (target in dimnames(rep_draws[[1]])[[3]]) {
      vals <- lapply(rep_draws, function(z) z[, r, target]); z <- unlist(vals); d <- study06_diag(vals)
      tv <- if (target == "active") 1 - truth_pi[reps[r], 1] else truth_pi[reps[r], match(target, paste0("pi", 0:3))]
      representative[[length(representative) + 1L]] <- data.frame(route, marker_id = rownames(A)[reps[r]], target,
        posterior_mean = mean(z), truth = tv, bias = mean(z) - tv, t(d))
    }
    mm <- marker[[route]]
    for (ch in seq_len(4L)) {
      eta_sd <- matrix(mm[[ch]]$eta_sd, nrow = nrow(A), ncol = 3L)
      for (stick in seq_len(3L)) eta_chain[[length(eta_chain) + 1L]] <- data.frame(route, chain = ch, stick,
        posterior_mean_across_markers = mean(mm[[ch]]$eta_mean[, stick]),
        posterior_mean_marker_sd = stats::sd(mm[[ch]]$eta_mean[, stick]),
        posterior_mean_q025 = stats::quantile(mm[[ch]]$eta_mean[, stick], .025),
        posterior_mean_median = stats::median(mm[[ch]]$eta_mean[, stick]),
        posterior_mean_q975 = stats::quantile(mm[[ch]]$eta_mean[, stick], .975),
        posterior_sd_across_markers = mean(eta_sd[, stick]),
        posterior_sd_q025 = stats::quantile(eta_sd[, stick], .025),
        posterior_sd_median = stats::median(eta_sd[, stick]),
        posterior_sd_q975 = stats::quantile(eta_sd[, stick], .975))
      checked_draws <- unique(c(1L, ceiling(nrow(alpha$traces[[ch]]) / 2), nrow(alpha$traces[[ch]])))
      for (draw in checked_draws) {
        alpha_draw <- matrix(alpha$traces[[ch]][draw, ], nrow = ncol(A), ncol = 3L)
        offline <- study06_pi(A %*% alpha_draw)
        package <- sblr::sbayesrc_marker_pi(A, alpha_draw, gamma = c(0, .01, .1, 1))
        probability_crosscheck[[length(probability_crosscheck) + 1L]] <- data.frame(route, chain = ch,
          draw, maximum_absolute_difference = max(abs(offline - package)),
          comparison = "offline exact transform versus sblr::sbayesrc_marker_pi")
      }
    }
    for (x in 1:3) for (y in (x + 1L):4) {
      quantities <- list(prior_active = 1 - mm[[x]]$pi_mean[, 1],
        largest_component_prior = mm[[x]]$pi_mean[, 4], pip = fits[[route]]$chains[[x]]$dm,
        absolute_posterior_beta = abs(fits[[route]]$chains[[x]]$bm))
      quantities_y <- list(prior_active = 1 - mm[[y]]$pi_mean[, 1],
        largest_component_prior = mm[[y]]$pi_mean[, 4], pip = fits[[route]]$chains[[y]]$dm,
        absolute_posterior_beta = abs(fits[[route]]$chains[[y]]$bm))
      for (nm in names(quantities)) ranking[[length(ranking) + 1L]] <- data.frame(route,
        chain_1 = x, chain_2 = y, quantity = nm,
        t(study06_pair_metrics(quantities[[nm]], quantities_y[[nm]], top = c(50L, 100L))))
      gx <- drop(truth$validation_x %*% fits[[route]]$chains[[x]]$bm)
      gy <- drop(truth$validation_x %*% fits[[route]]$chains[[y]]$bm)
      gv[[length(gv) + 1L]] <- data.frame(route, chain_1 = x, chain_2 = y,
        correlation = stats::cor(gx, gy), rmse = sqrt(mean((gx - gy)^2)))
    }
  }
  raw <- do.call(rbind, raw); raw_chain <- do.call(rbind, raw_chain)
  expected <- do.call(rbind, expected); representative <- do.call(rbind, representative)
  ranking <- do.call(rbind, ranking); gv <- do.call(rbind, gv)
  eta_chain <- do.call(rbind, eta_chain); probability_crosscheck <- do.call(rbind, probability_crosscheck)
  study06_write_csv(raw, file.path(td, "raw_alpha_convergence.csv"))
  study06_write_csv(raw_chain, file.path(td, "raw_alpha_chain.csv"))
  study06_write_csv(expected, file.path(td, "expected_active_count.csv"))
  study06_write_csv(representative, file.path(td, "representative_prior_probability_diagnostics.csv"))
  study06_write_csv(ranking, file.path(td, "ranking_stability.csv"))
  study06_write_csv(gv, file.path(td, "genetic_value_stability.csv"))
  study06_write_csv(eta_chain, file.path(td, "A_alpha_chain_summary.csv"))
  study06_write_csv(probability_crosscheck, file.path(td, "probability_package_crosscheck.csv"))
  study06_write_csv(read.csv(file.path(td, "selected_A_alpha_diagnostics.csv"), check.names = FALSE),
    file.path(td, "selected_alpha_linear_contrasts.csv"))

  eta <- read.csv(file.path(td, "A_alpha_stability.csv"), check.names = FALSE)
  pi_stab <- read.csv(file.path(td, "prior_probability_stability.csv"), check.names = FALSE)
  recovery <- read.csv(file.path(td, "derived_truth_recovery.csv"), check.names = FALSE)
  contrast <- read.csv(file.path(td, "counterfactual_annotation_contrasts.csv"), check.names = FALSE)
  snp <- read.csv(file.path(td, "SNP_chain_stability.csv"), check.names = FALSE)
  perf <- read.csv(file.path(td, "SNP_chain_performance.csv"), check.names = FALSE)
  selected <- read.csv(file.path(td, "selected_A_alpha_diagnostics.csv"), check.names = FALSE)
  hierarchy_row <- function(route, level, d = NULL, corr = NA_real_, truth_error = NA_real_) {
    data.frame(route, level,
      worst_rhat = if (is.null(d) || !"rhat" %in% names(d)) NA else max(d$rhat, na.rm = TRUE),
      median_rhat = if (is.null(d) || !"rhat" %in% names(d)) NA else stats::median(d$rhat, na.rm = TRUE),
      minimum_bulk_ess = if (is.null(d) || !"bulk_ess" %in% names(d)) NA else min(d$bulk_ess, na.rm = TRUE),
      maximum_relative_mcse = if (is.null(d) || !"relative_mcse" %in% names(d)) NA else max(d$relative_mcse, na.rm = TRUE),
      minimum_between_chain_correlation = corr, truth_error = truth_error)
  }
  hierarchy <- list()
  for (route in names(fits)) {
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "raw alpha", raw[raw$route == route, ], truth_error = max(raw$absolute_error[raw$route == route]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "selected alpha linear contrasts", selected[selected$route == route, ], truth_error = max(abs(selected$bias[selected$route == route])))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "A alpha", selected[selected$route == route, ], min(eta$pearson[eta$route == route]), max(recovery$rmse[recovery$route == route & recovery$level == "eta"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "prior active probability", expected[expected$route == route, ], min(pi_stab$pearson[pi_stab$route == route & pi_stab$quantity == "p_active"]), abs(expected$bias[expected$route == route]) / nrow(A))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "component probabilities", representative[representative$route == route, ], min(pi_stab$pearson[pi_stab$route == route & grepl("^pi", pi_stab$quantity)]), max(recovery$rmse[recovery$route == route & recovery$level == "pi"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "counterfactual annotation contrasts", contrast[contrast$route == route, ], truth_error = max(abs(contrast$bias[contrast$route == route])))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "SNP PIP", corr = min(snp$pearson[snp$route == route & snp$quantity == "pip"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "posterior beta", corr = min(snp$pearson[snp$route == route & snp$quantity == "beta_signed"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "genetic value", corr = min(gv$correlation[gv$route == route]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "prediction", truth_error = 1 - min(perf$phenotype_prediction_correlation[perf$route == route]))
  }
  hierarchy <- do.call(rbind, hierarchy)
  study06_write_csv(hierarchy, file.path(td, "hierarchy_of_stability.csv"))

  official_raw <- read.csv(file.path(td, "official_alpha_summary.csv"), check.names = FALSE)
  official_contrast <- read.csv(file.path(td, "official_contrast_summary.csv"), check.names = FALSE)
  contrast_truth <- function(annotation, target) {
    pair <- if (annotation == "enriched_binary") c(0, 1) else c(-1, 1)
    aa <- A; aa[, annotation] <- pair[1]; p0 <- study06_pi(aa %*% truth_alpha)
    aa[, annotation] <- pair[2]; p1 <- study06_pi(aa %*% truth_alpha)
    values <- c(active = mean(1 - p1[, 1]) - mean(1 - p0[, 1]), colMeans(p1 - p0))
    unname(values[target])
  }
  official_contrast$truth <- mapply(contrast_truth, official_contrast$annotation, official_contrast$target)
  official_contrast$bias <- official_contrast$mean - official_contrast$truth
  official_contrast$coverage <- official_contrast$truth >= official_contrast$q025 & official_contrast$truth <= official_contrast$q975
  official_contrast$p_lt_0 <- 1 - official_contrast$p_gt_0
  official_contrast$relative_mcse <- 1 / sqrt(official_contrast$ess)
  study06_write_csv(official_contrast, file.path(td, "official_contrast_summary.csv"))
  ofiles <- file.path(root, "results/local/06_annotation_models/gctb_single_trajectory/runs/D1/D1.mcmcsamples",
    paste0("D1.mcmcsamples.AnnoEffects_p", 1:3))
  otr <- lapply(ofiles, function(f) as.matrix(read.table(f, header = TRUE, check.names = FALSE)))
  n <- min(vapply(otr, nrow, integer(1))); half <- floor(n/2); official_level <- list()
  for (part in list(first = seq_len(half), second = (n-half+1L):n)) {
    am <- sapply(otr, function(z) colMeans(z[part, , drop = FALSE])); rownames(am) <- colnames(otr[[1]])
    ep <- A %*% am; pp <- study06_pi(ep)
    official_level[[length(official_level)+1L]] <- list(eta = ep, pi = pp)
  }
  official_compare <- data.frame(level = c("raw alpha", "A alpha", "component probabilities", "annotation contrasts"),
    official_evidence = c(sprintf("single-chain ESS %.1f-%.1f; max half shift %.2f SD", min(official_raw$ess), max(official_raw$ess), max(abs(official_raw$half_mean_shift_sd))),
      sprintf("half-vector minimum Pearson %.3f", min(vapply(1:3, function(k) cor(official_level[[1]]$eta[,k], official_level[[2]]$eta[,k]), numeric(1)))),
      sprintf("half-vector minimum Pearson %.3f", min(vapply(1:4, function(k) suppressWarnings(cor(official_level[[1]]$pi[,k], official_level[[2]]$pi[,k])), numeric(1)), na.rm = TRUE)),
      sprintf("active contrast max half shift %.2f SD", max(abs(official_contrast$half_mean_shift_sd[official_contrast$target == "active"])))),
    sblr_evidence = c(sprintf("max R-hat %.3f", max(raw$rhat)),
      sprintf("selected eta max R-hat %.3f; min vector Pearson %.3f", max(selected$rhat), min(eta$pearson)),
      sprintf("representative max R-hat %.3f; min vector Pearson %.3f", max(representative$rhat), min(pi_stab$pearson)),
      sprintf("max R-hat %.3f", max(contrast$rhat))),
    limitation = c("official trajectories cannot be independently seeded", rep("single official trajectory only", 3)))
  study06_write_csv(official_compare, file.path(td, "official_sbayesrc_comparison.csv"))

  p1 <- ggplot2::ggplot(raw_chain, ggplot2::aes(chain, mean, colour = factor(chain))) + ggplot2::geom_point() +
    ggplot2::facet_grid(annotation ~ route + stick, scales = "free_y") + ggplot2::theme_bw() + ggplot2::guides(colour = "none")
  ggplot2::ggsave(file.path(out, "figures/alpha_chain_disagreement.png"), p1, width = 12, height = 8, dpi = 160)
  hp <- rbind(data.frame(route = hierarchy$route, level = hierarchy$level, metric = "R-hat", value = hierarchy$worst_rhat),
    data.frame(route = hierarchy$route, level = hierarchy$level, metric = "between-chain correlation", value = hierarchy$minimum_between_chain_correlation))
  p5 <- ggplot2::ggplot(hp[is.finite(hp$value), ], ggplot2::aes(level, value, fill = route)) + ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~metric, scales = "free_y") + ggplot2::coord_flip() + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(out, "figures/hierarchy_of_stability.png"), p5, width = 10, height = 8, dpi = 160)
  invisible(list(hierarchy = hierarchy, official = official_compare))
}
