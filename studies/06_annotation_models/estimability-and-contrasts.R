# Study 06 final estimability and annotation-contrast analysis.
# This script consumes frozen retained chains; it never runs a sampler.

study06_pi <- function(eta) {
  stopifnot(is.matrix(eta), ncol(eta) == 3L, all(is.finite(eta)))
  continuation <- stats::pnorm(eta)
  remaining <- rep(1, nrow(eta))
  out <- matrix(0, nrow(eta), 4L)
  for (stick in seq_len(3L)) {
    out[, stick] <- remaining * (1 - continuation[, stick])
    remaining <- remaining * continuation[, stick]
  }
  out[, 4L] <- remaining
  colnames(out) <- paste0("pi", 0:3)
  out
}

study06_diag <- function(chains) {
  n <- min(vapply(chains, length, integer(1)))
  mat <- do.call(cbind, lapply(chains, function(x) as.numeric(x)[seq_len(n)]))
  draws <- posterior::as_draws_array(array(mat, dim = c(n, ncol(mat), 1L)))
  pooled <- as.numeric(mat)
  sd_value <- stats::sd(pooled)
  mcse <- posterior::mcse_mean(draws)
  c(rhat = posterior::rhat(draws), bulk_ess = posterior::ess_bulk(draws),
    tail_ess = posterior::ess_tail(draws), mcse = mcse,
    relative_mcse = if (is.finite(sd_value) && sd_value > 0) mcse / sd_value else NA_real_)
}

study06_pair_metrics <- function(x, y, top = c(25L, 50L, 100L)) {
  out <- c(pearson = stats::cor(x, y), spearman = stats::cor(x, y,
    method = "spearman"), rmse = sqrt(mean((x - y)^2)),
    mae = mean(abs(x - y)), max_abs = max(abs(x - y)))
  for (k in top) {
    ix <- head(order(x, decreasing = TRUE), min(k, length(x)))
    iy <- head(order(y, decreasing = TRUE), min(k, length(y)))
    out[paste0("top", k, "_overlap")] <- length(intersect(ix, iy)) / min(k, length(x))
  }
  out
}

study06_auc <- function(score, truth) {
  truth <- as.logical(truth)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - sum(seq_len(sum(truth)))) / (sum(truth) * sum(!truth))
}

study06_auprc <- function(score, truth) {
  truth <- as.logical(truth)
  ord <- order(score, decreasing = TRUE)
  hit <- truth[ord]
  recall <- cumsum(hit) / sum(hit)
  precision <- cumsum(hit) / seq_along(hit)
  sum((recall - c(0, head(recall, -1L))) * precision)
}

study06_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

study06_chain_alpha <- function(fit) {
  quantities <- fit$convergence_traces$quantities
  descriptor <- quantities[quantities$parameter_name == "alpha" &
    quantities$annotation_index > 0L & quantities$stick_index > 0L, , drop = FALSE]
  descriptor <- descriptor[order(descriptor$quantity_index), , drop = FALSE]
  traces <- lapply(fit$chains, function(chain) {
    x <- chain$convergence_trace$alpha
    stopifnot(ncol(x) == nrow(descriptor), nrow(x) == chain$retained_draw_count)
    colnames(x) <- paste(descriptor$annotation_name, descriptor$stick_name, sep = "::")
    x
  })
  list(traces = traces, descriptor = descriptor)
}

study06_raw_alpha <- function(route, alpha, truth) {
  desc <- alpha$descriptor
  rows <- vector("list", nrow(desc))
  chain_rows <- list()
  for (j in seq_len(nrow(desc))) {
    values <- lapply(alpha$traces, function(x) x[, j])
    pooled <- unlist(values, use.names = FALSE)
    target <- truth[desc$annotation_index[j], desc$stick_index[j]]
    dg <- study06_diag(values)
    rows[[j]] <- data.frame(route, stick = desc$stick_name[j],
      annotation = desc$annotation_name[j], posterior_mean = mean(pooled),
      posterior_sd = stats::sd(pooled), median = stats::median(pooled),
      q025 = stats::quantile(pooled, .025), q975 = stats::quantile(pooled, .975),
      truth = target, bias = mean(pooled) - target,
      absolute_error = abs(mean(pooled) - target), t(dg), check.names = FALSE)
    chain_rows[[j]] <- do.call(rbind, lapply(seq_along(values), function(ch) {
      z <- values[[ch]]
      data.frame(route, chain = ch, stick = desc$stick_name[j],
        annotation = desc$annotation_name[j], mean = mean(z), sd = stats::sd(z),
        q025 = stats::quantile(z, .025), q975 = stats::quantile(z, .975))
    }))
  }
  list(summary = do.call(rbind, rows), chains = do.call(rbind, chain_rows))
}

study06_design_dependence <- function(route, A, alpha) {
  design <- data.frame(route, statistic = c("rank", "condition_number",
      paste0("singular_value_", seq_len(ncol(A)))),
    value = c(qr(A)$rank, kappa(A), svd(A, nu = 0, nv = 0)$d))
  correlation <- as.data.frame(as.table(stats::cor(A)))
  names(correlation) <- c("annotation_1", "annotation_2", "correlation")
  correlation$route <- route
  covariance <- correlation_rows <- eigen_rows <- list()
  q <- alpha$descriptor
  for (stick in seq_len(3L)) {
    idx <- which(q$stick_index == stick)
    mats <- lapply(alpha$traces, function(x) x[, idx, drop = FALSE])
    names_here <- q$annotation_name[idx]
    for (ch in seq_along(mats)) {
      cv <- stats::cov(mats[[ch]]); cr <- stats::cor(mats[[ch]])
      covariance[[length(covariance) + 1L]] <- transform(as.data.frame(as.table(cv)),
        route = route, chain = ch, stick = stick, matrix = "covariance")
      correlation_rows[[length(correlation_rows) + 1L]] <- transform(as.data.frame(as.table(cr)),
        route = route, chain = ch, stick = stick, matrix = "correlation")
    }
    pooled <- do.call(rbind, mats); cv <- stats::cov(pooled); cr <- stats::cor(pooled)
    covariance[[length(covariance) + 1L]] <- transform(as.data.frame(as.table(cv)),
      route = route, chain = 0L, stick = stick, matrix = "covariance")
    correlation_rows[[length(correlation_rows) + 1L]] <- transform(as.data.frame(as.table(cr)),
      route = route, chain = 0L, stick = stick, matrix = "correlation")
    eg <- eigen(cv, symmetric = TRUE)
    for (direction in seq_along(eg$values)) {
      for (a in seq_along(names_here)) eigen_rows[[length(eigen_rows) + 1L]] <-
        data.frame(route, stick, direction, eigenvalue = eg$values[direction],
          annotation = names_here[a], loading = eg$vectors[a, direction])
    }
  }
  list(design = design, design_correlation = correlation,
    posterior_matrices = rbind(do.call(rbind, covariance), do.call(rbind, correlation_rows)),
    eigen = do.call(rbind, eigen_rows))
}

study06_representatives <- function(A) {
  target <- expand.grid(enriched_binary = 0:1,
    continuous_signal = c(-1, 0, 1), null_annotation = 0)
  unique(vapply(seq_len(nrow(target)), function(i) {
    distance <- rowSums((A[, c("enriched_binary", "continuous_signal", "null_annotation")] -
      as.numeric(target[i, ]))^2)
    which.min(distance)
  }, integer(1)))
}

study06_derive_route <- function(route, fit, alpha, A, truth_alpha, truth_pi,
                                  chunk_size = 500L) {
  m <- nrow(A); q <- alpha$descriptor; representative <- study06_representatives(A)
  mods <- list(enriched_binary = c(0, 1), continuous_signal = c(-1, 1),
    null_annotation = c(-1, 1))
  marker <- vector("list", length(alpha$traces)); scalar <- vector("list", length(alpha$traces))
  contrast <- vector("list", length(alpha$traces)); observed <- vector("list", length(alpha$traces))
  normalization <- list()
  for (ch in seq_along(alpha$traces)) {
    message(sprintf("%s: deriving chain %d/%d", route, ch, length(alpha$traces)))
    tr <- alpha$traces[[ch]]; nd <- nrow(tr)
    eta_sum <- matrix(0, m, 3L); eta_sq <- matrix(0, m, 3L)
    pi_sum <- matrix(0, m, 4L); pi_sq <- matrix(0, m, 4L)
    eta_selected <- array(NA_real_, c(nd, length(representative), 3L))
    contrast_draw <- array(NA_real_, c(nd, length(mods), 5L),
      dimnames = list(NULL, names(mods), c("active", paste0("pi", 0:3))))
    observed_draw <- matrix(NA_real_, nd, 5L,
      dimnames = list(NULL, c("active", paste0("pi", 0:3))))
    max_row_error <- 0; min_probability <- 1; max_probability <- 0
    for (lo in seq.int(1L, nd, by = chunk_size)) {
      ix <- lo:min(nd, lo + chunk_size - 1L); nc <- length(ix)
      eta <- array(NA_real_, c(m, nc, 3L)); coef <- vector("list", 3L)
      for (stick in seq_len(3L)) {
        idx <- which(q$stick_index == stick)
        coef[[stick]] <- tr[ix, idx, drop = FALSE]
        eta[, , stick] <- A %*% t(coef[[stick]])
        eta_sum[, stick] <- eta_sum[, stick] + rowSums(eta[, , stick])
        eta_sq[, stick] <- eta_sq[, stick] + rowSums(eta[, , stick]^2)
        eta_selected[ix, , stick] <- t(eta[representative, , stick, drop = FALSE][, , 1L])
      }
      cont <- lapply(seq_len(3L), function(stick) stats::pnorm(eta[, , stick]))
      pi <- array(0, c(m, nc, 4L)); rem <- matrix(1, m, nc)
      for (stick in seq_len(3L)) { pi[, , stick] <- rem * (1 - cont[[stick]]); rem <- rem * cont[[stick]] }
      pi[, , 4L] <- rem
      for (k in seq_len(4L)) { pi_sum[, k] <- pi_sum[, k] + rowSums(pi[, , k]); pi_sq[, k] <- pi_sq[, k] + rowSums(pi[, , k]^2) }
      row_error <- abs(apply(pi, c(1, 2), sum) - 1)
      max_row_error <- max(max_row_error, row_error); min_probability <- min(min_probability, pi)
      max_probability <- max(max_probability, pi)
      enriched <- A[, "enriched_binary"] == 1
      for (k in seq_len(4L)) observed_draw[ix, k + 1L] <- colMeans(pi[enriched, , k]) - colMeans(pi[!enriched, , k])
      observed_draw[ix, 1L] <- -observed_draw[ix, 2L]
      for (a in seq_along(mods)) {
        annotation <- names(mods)[a]; pair <- mods[[a]]; pi_cf <- vector("list", 2L)
        annotation_index <- match(annotation, colnames(A))
        for (side in 1:2) {
          cfs <- lapply(seq_len(3L), function(stick) stats::pnorm(eta[, , stick] +
            outer(pair[side] - A[, annotation], coef[[stick]][, annotation_index])))
          z <- array(0, c(m, nc, 4L)); rr <- matrix(1, m, nc)
          for (stick in seq_len(3L)) { z[, , stick] <- rr * (1 - cfs[[stick]]); rr <- rr * cfs[[stick]] }
          z[, , 4L] <- rr; pi_cf[[side]] <- z
        }
        for (k in seq_len(4L)) contrast_draw[ix, a, k + 1L] <-
          colMeans(pi_cf[[2L]][, , k]) - colMeans(pi_cf[[1L]][, , k])
        contrast_draw[ix, a, 1L] <- -contrast_draw[ix, a, 2L]
      }
    }
    colnames(eta_sum) <- paste0("stick_", 1:3); colnames(pi_sum) <- paste0("pi", 0:3)
    marker[[ch]] <- list(eta_mean = eta_sum / nd,
      eta_sd = sqrt(pmax(0, (eta_sq - eta_sum^2 / nd) / (nd - 1))),
      pi_mean = pi_sum / nd, pi_sd = sqrt(pmax(0, (pi_sq - pi_sum^2 / nd) / (nd - 1))))
    scalar[[ch]] <- eta_selected; contrast[[ch]] <- contrast_draw; observed[[ch]] <- observed_draw
    normalization[[ch]] <- data.frame(route, chain = ch, min_probability,
      max_probability, max_row_sum_error = max_row_error, all_finite = TRUE)
    message(sprintf("%s: completed chain %d/%d", route, ch, length(alpha$traces)))
  }
  eta_stability <- pi_stability <- list()
  for (x in seq_len(length(marker) - 1L)) for (y in (x + 1L):length(marker)) {
    for (stick in seq_len(3L)) eta_stability[[length(eta_stability) + 1L]] <-
      data.frame(route, chain_1 = x, chain_2 = y, quantity = paste0("eta_stick_", stick),
        t(study06_pair_metrics(marker[[x]]$eta_mean[, stick], marker[[y]]$eta_mean[, stick])))
    for (k in seq_len(4L)) pi_stability[[length(pi_stability) + 1L]] <-
      data.frame(route, chain_1 = x, chain_2 = y, quantity = paste0("pi", k - 1L),
        t(study06_pair_metrics(marker[[x]]$pi_mean[, k], marker[[y]]$pi_mean[, k])))
    pi_stability[[length(pi_stability) + 1L]] <- data.frame(route, chain_1 = x,
      chain_2 = y, quantity = "p_active", t(study06_pair_metrics(
        1 - marker[[x]]$pi_mean[, 1L], 1 - marker[[y]]$pi_mean[, 1L])))
  }
  eta_scalar <- list()
  for (r in seq_along(representative)) for (stick in seq_len(3L)) {
    vals <- lapply(scalar, function(z) z[, r, stick]); dg <- study06_diag(vals)
    target <- sum(A[representative[r], ] * truth_alpha[, stick])
    eta_scalar[[length(eta_scalar) + 1L]] <- data.frame(route,
      marker_id = rownames(A)[representative[r]], stick, truth = target,
      posterior_mean = mean(unlist(vals)), bias = mean(unlist(vals)) - target, t(dg))
  }
  contrast_table <- list(); observed_table <- list()
  truth_contrast <- function(annotation, target) {
    pair <- mods[[annotation]]; aa <- A; aa[, annotation] <- pair[1]; p0 <- study06_pi(aa %*% truth_alpha)
    aa[, annotation] <- pair[2]; p1 <- study06_pi(aa %*% truth_alpha)
    values <- c(active = mean(1 - p1[, 1]) - mean(1 - p0[, 1]), colMeans(p1 - p0))
    values[target]
  }
  for (a in seq_along(mods)) for (target in dimnames(contrast[[1]])[[3]]) {
    vals <- lapply(contrast, function(z) z[, a, target]); pooled <- unlist(vals)
    dg <- study06_diag(vals); tv <- truth_contrast(names(mods)[a], target)
    contrast_table[[length(contrast_table) + 1L]] <- data.frame(route,
      annotation = names(mods)[a], target, comparison = if (a == 1L) "1_vs_0" else "+1SD_vs_-1SD",
      posterior_mean = mean(pooled), posterior_sd = stats::sd(pooled), median = stats::median(pooled),
      q025 = stats::quantile(pooled, .025), q975 = stats::quantile(pooled, .975),
      p_gt_0 = mean(pooled > 0), p_lt_0 = mean(pooled < 0), t(dg), truth = tv,
      bias = mean(pooled) - tv, coverage = tv >= stats::quantile(pooled, .025) && tv <= stats::quantile(pooled, .975))
    ovals <- lapply(observed, function(z) z[, target]); op <- unlist(ovals); od <- study06_diag(ovals)
    observed_table[[length(observed_table) + 1L]] <- data.frame(route,
      annotation = "enriched_binary", target, posterior_mean = mean(op),
      q025 = stats::quantile(op, .025), q975 = stats::quantile(op, .975), t(od))
  }
  pooled_eta <- Reduce(`+`, lapply(marker, `[[`, "eta_mean")) / length(marker)
  pooled_pi <- Reduce(`+`, lapply(marker, `[[`, "pi_mean")) / length(marker)
  eta_truth <- A %*% truth_alpha
  recovery <- rbind(data.frame(route, level = "eta", quantity = paste0("stick_", 1:3),
    rmse = sqrt(colMeans((pooled_eta - eta_truth)^2)), mae = colMeans(abs(pooled_eta - eta_truth))),
    data.frame(route, level = "pi", quantity = paste0("pi", 0:3),
      rmse = sqrt(colMeans((pooled_pi - truth_pi)^2)), mae = colMeans(abs(pooled_pi - truth_pi))))
  list(marker = marker, eta_stability = do.call(rbind, eta_stability),
    pi_stability = do.call(rbind, pi_stability), eta_scalar = do.call(rbind, eta_scalar),
    contrasts = do.call(rbind, contrast_table), observed = unique(do.call(rbind, observed_table)),
    normalization = do.call(rbind, normalization), recovery = recovery,
    representatives = data.frame(route, marker_index = representative, marker_id = rownames(A)[representative]))
}

study06_snp_stability <- function(route, fit, truth) {
  rows <- list(); chain_metrics <- list(); rank_rows <- list()
  for (ch in seq_along(fit$chains)) {
    z <- fit$chains[[ch]]; pip <- z$dm; beta <- z$bm
    gv <- drop(truth$validation_x %*% beta)
    chain_metrics[[ch]] <- data.frame(route, chain = ch,
      auprc = study06_auprc(pip, truth$marker_truth$true_nonnull),
      auroc = study06_auc(pip, truth$marker_truth$true_nonnull),
      validation_genetic_correlation = stats::cor(gv, truth$validation_genetic),
      phenotype_prediction_correlation = stats::cor(gv, truth$validation_phenotype))
  }
  for (x in seq_len(length(fit$chains) - 1L)) for (y in (x + 1L):length(fit$chains)) {
    for (quantity in c("pip", "abs_beta")) {
      a <- if (quantity == "pip") fit$chains[[x]]$dm else abs(fit$chains[[x]]$bm)
      b <- if (quantity == "pip") fit$chains[[y]]$dm else abs(fit$chains[[y]]$bm)
      rows[[length(rows) + 1L]] <- data.frame(route, chain_1 = x, chain_2 = y, quantity,
        t(study06_pair_metrics(a, b, top = c(50L, 100L))))
    }
    bx <- fit$chains[[x]]$bm; by <- fit$chains[[y]]$bm
    rows[[length(rows) + 1L]] <- data.frame(route, chain_1 = x, chain_2 = y,
      quantity = "beta_signed", pearson = stats::cor(bx, by),
      spearman = stats::cor(bx, by, method = "spearman"), rmse = sqrt(mean((bx - by)^2)),
      mae = mean(abs(bx - by)), max_abs = max(abs(bx - by)), top50_overlap = NA, top100_overlap = NA)
  }
  list(stability = do.call(rbind, rows), performance = do.call(rbind, chain_metrics))
}

study06_official <- function(A, truth_alpha, base) {
  files <- file.path(base, paste0("D1.mcmcsamples.AnnoEffects_p", 1:3))
  if (!all(file.exists(files))) return(list(status = "alpha traces unavailable"))
  traces <- lapply(files, function(f) as.matrix(utils::read.table(f, header = TRUE,
    check.names = FALSE)))
  n <- min(vapply(traces, nrow, integer(1))); traces <- lapply(traces, function(x) tail(x, n))
  raw <- list()
  for (stick in seq_len(3L)) for (a in seq_len(4L)) {
    z <- traces[[stick]][, a]; first <- head(z, floor(n / 2)); last <- tail(z, floor(n / 2))
    raw[[length(raw) + 1L]] <- data.frame(stick, annotation = colnames(traces[[stick]])[a],
      mean = mean(z), sd = stats::sd(z), q025 = stats::quantile(z, .025), q975 = stats::quantile(z, .975),
      truth = truth_alpha[a, stick], bias = mean(z) - truth_alpha[a, stick],
      ess = coda::effectiveSize(z), relative_mcse = stats::sd(z) / sqrt(coda::effectiveSize(z)) / stats::sd(z),
      half_mean_shift_sd = (mean(last) - mean(first)) / stats::sd(z))
  }
  contrast <- list(); chunk <- 500L
  mods <- list(enriched_binary = c(0, 1), continuous_signal = c(-1, 1), null_annotation = c(-1, 1))
  draws <- array(NA_real_, c(n, length(mods), 5L), dimnames = list(NULL, names(mods), c("active", paste0("pi", 0:3))))
  for (lo in seq.int(1L, n, by = chunk)) {
    ix <- lo:min(n, lo + chunk - 1L)
    for (a in seq_along(mods)) {
      ps <- vector("list", 2L)
      for (side in 1:2) {
        aa <- A; aa[, names(mods)[a]] <- mods[[a]][side]
        eta <- array(NA_real_, c(nrow(A), length(ix), 3L))
        for (stick in 1:3) eta[, , stick] <- aa %*% t(traces[[stick]][ix, , drop = FALSE])
        cc <- lapply(1:3, function(s) stats::pnorm(eta[, , s])); p <- array(0, c(nrow(A), length(ix), 4L)); rem <- matrix(1, nrow(A), length(ix))
        for (s in 1:3) { p[, , s] <- rem * (1 - cc[[s]]); rem <- rem * cc[[s]] }; p[, , 4] <- rem; ps[[side]] <- p
      }
      for (k in 1:4) draws[ix, a, k + 1L] <- colMeans(ps[[2]][, , k]) - colMeans(ps[[1]][, , k])
      draws[ix, a, 1L] <- -draws[ix, a, 2L]
    }
  }
  for (a in seq_along(mods)) for (target in dimnames(draws)[[3]]) {
    z <- draws[, a, target]; first <- head(z, floor(n / 2)); last <- tail(z, floor(n / 2))
    contrast[[length(contrast) + 1L]] <- data.frame(annotation = names(mods)[a], target,
      mean = mean(z), sd = stats::sd(z), q025 = stats::quantile(z, .025), q975 = stats::quantile(z, .975),
      p_gt_0 = mean(z > 0), ess = coda::effectiveSize(z),
      half_mean_shift_sd = (mean(last) - mean(first)) / stats::sd(z))
  }
  list(status = "single deterministic native trajectory; multitrajectory inference unavailable",
    raw = do.call(rbind, raw), contrasts = do.call(rbind, contrast))
}

run_study06_estimability <- function(root = ".") {
  stopifnot(basename(normalizePath(root)) == "sblrbench",
    as.character(utils::packageVersion("sblr")) == "0.2.0")
  suppressPackageStartupMessages({ library(posterior); library(ggplot2) })
  source(file.path(root, "studies/06_annotation_models/spec.R"), local = TRUE)
  source(file.path(root, "studies/06_annotation_models/annotation-design.R"), local = TRUE)
  out <- file.path(root, "results/local/06_annotation_models/estimability_and_contrasts")
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
  truth_path <- file.path(root, "results/local/06_annotation_models/gctb_parity/export/truth.rds")
  truth <- readRDS(truth_path); A <- truth$annotations
  stopifnot(identical(rownames(A), truth$marker_ids), identical(colnames(A), spec$annotation_design$columns))
  truth_object <- construct_annotation_truth(A, spec); truth_alpha <- truth_object$informative_annotations
  truth_pi <- study06_pi(A %*% truth_alpha)
  stopifnot(max(abs(truth_pi - as.matrix(truth$marker_truth[paste0("true_prior_component_", 0:3)]))) < 1e-10)
  fit_paths <- c(bed = file.path(root, "results/local/06_annotation_models/v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_bed_bayesrc.rds"),
    block = file.path(root, "results/local/06_annotation_models/v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_block_eigen_sbayesrc.rds"))
  fits <- lapply(fit_paths, function(p) readRDS(p)$result$native_fit)
  all_raw <- all_raw_chain <- all_design <- all_design_cor <- all_postmat <- all_eigen <- list()
  all_eta <- all_pi <- all_eta_scalar <- all_contrast <- all_observed <- all_norm <- all_recovery <- all_snp <- all_perf <- list(); derived <- list()
  for (route in names(fits)) {
    message(sprintf("Starting retained-draw analysis for %s", route))
    alpha <- study06_chain_alpha(fits[[route]])
    raw <- study06_raw_alpha(route, alpha, truth_alpha); dep <- study06_design_dependence(route, A, alpha)
    drv <- study06_derive_route(route, fits[[route]], alpha, A, truth_alpha, truth_pi); derived[[route]] <- drv
    snp <- study06_snp_stability(route, fits[[route]], truth)
    all_raw[[route]] <- raw$summary; all_raw_chain[[route]] <- raw$chains
    all_design[[route]] <- dep$design; all_design_cor[[route]] <- dep$design_correlation
    all_postmat[[route]] <- dep$posterior_matrices; all_eigen[[route]] <- dep$eigen
    all_eta[[route]] <- drv$eta_stability; all_pi[[route]] <- drv$pi_stability; all_eta_scalar[[route]] <- drv$eta_scalar
    all_contrast[[route]] <- drv$contrasts; all_observed[[route]] <- drv$observed; all_norm[[route]] <- drv$normalization
    all_recovery[[route]] <- drv$recovery; all_snp[[route]] <- snp$stability; all_perf[[route]] <- snp$performance
    message(sprintf("Completed retained-draw analysis for %s", route))
  }
  tables <- list(raw_alpha_convergence = do.call(rbind, all_raw), raw_alpha_chain = do.call(rbind, all_raw_chain),
    annotation_design = do.call(rbind, all_design), annotation_design_correlations = do.call(rbind, all_design_cor),
    alpha_posterior_matrices = do.call(rbind, all_postmat), alpha_eigendirections = do.call(rbind, all_eigen),
    A_alpha_stability = do.call(rbind, all_eta), selected_A_alpha_diagnostics = do.call(rbind, all_eta_scalar),
    prior_probability_stability = do.call(rbind, all_pi), probability_normalization = do.call(rbind, all_norm),
    derived_truth_recovery = do.call(rbind, all_recovery), counterfactual_annotation_contrasts = do.call(rbind, all_contrast),
    observed_group_contrasts = do.call(rbind, all_observed), SNP_chain_stability = do.call(rbind, all_snp),
    SNP_chain_performance = do.call(rbind, all_perf))
  official_base <- file.path(root, "results/local/06_annotation_models/gctb_single_trajectory/runs/D1/D1.mcmcsamples")
  official <- study06_official(A, truth_alpha, official_base)
  message("Completed official single-trajectory derived summaries")
  if (!is.null(official$raw)) tables$official_alpha_summary <- official$raw
  if (!is.null(official$contrasts)) tables$official_contrast_summary <- official$contrasts
  for (nm in names(tables)) study06_write_csv(tables[[nm]], file.path(out, "tables", paste0(nm, ".csv")))
  saveRDS(lapply(derived, `[[`, "marker"), file.path(out, "marker_derived_chain_summaries.rds"))
  manifest <- data.frame(role = c("truth", names(fit_paths), "official_D1_alpha_p1", "official_D1_alpha_p2", "official_D1_alpha_p3"),
    path = c(truth_path, fit_paths, file.path(official_base, paste0("D1.mcmcsamples.AnnoEffects_p", 1:3))), stringsAsFactors = FALSE)
  manifest$bytes <- file.info(manifest$path)$size
  manifest$sha256 <- vapply(manifest$path, digest::digest, character(1), algo = "sha256", file = TRUE)
  manifest$specification_hash <- "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56"
  manifest$truth_hash <- "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb"
  study06_write_csv(manifest, file.path(out, "analysis_input_manifest.csv"))
  p1 <- ggplot(tables$raw_alpha_chain, aes(chain, mean, colour = factor(chain))) + geom_point() +
    facet_grid(annotation ~ route + stick, scales = "free_y") + theme_bw() + guides(colour = "none")
  ggsave(file.path(out, "figures/alpha_chain_disagreement.png"), p1, width = 12, height = 8, dpi = 160)
  p2 <- ggplot(tables$A_alpha_stability, aes(pearson, rmse, colour = route)) + geom_point() + facet_wrap(~quantity, scales = "free") + theme_bw()
  ggsave(file.path(out, "figures/A_alpha_chain_agreement.png"), p2, width = 9, height = 5, dpi = 160)
  p3 <- ggplot(tables$prior_probability_stability, aes(pearson, rmse, colour = route)) + geom_point() + facet_wrap(~quantity, scales = "free") + theme_bw()
  ggsave(file.path(out, "figures/marker_prior_probability_agreement.png"), p3, width = 9, height = 5, dpi = 160)
  p4 <- ggplot(tables$counterfactual_annotation_contrasts[tables$counterfactual_annotation_contrasts$target == "active", ],
    aes(annotation, posterior_mean, ymin = q025, ymax = q975, colour = route)) + geom_pointrange(position = position_dodge(.4)) + geom_hline(yintercept = 0, linetype = 2) + theme_bw()
  ggsave(file.path(out, "figures/annotation_contrast_posteriors.png"), p4, width = 8, height = 5, dpi = 160)
  invisible(list(tables = tables, derived = derived, official = official, manifest = manifest))
}

if (sys.nframe() == 0L) run_study06_estimability()
