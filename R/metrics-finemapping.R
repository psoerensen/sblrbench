# Fine-mapping summaries preserve Study 01's native marginal-PIP credible-set
# semantics. Native sets are constructed by sblr; these helpers only normalize
# and evaluate their tidy representation.

finemapping_marker_table <- function(result, simulation, method) {
  effects <- extract_marker_effects(result)
  probabilities <- extract_marker_probabilities(result)$posterior_inclusion
  if (is.null(probabilities))
    stop("Posterior inclusion probabilities are required for fine-mapping.",
      call. = FALSE)
  probabilities <- align_traits(align_markers(probabilities, rownames(effects)),
    colnames(effects))
  causal <- simulation$truth$causal$all
  out <- data.frame(study = simulation$scenario$study,
    scenario = simulation$scenario$architecture,
    replicate = simulation$scenario$replicate, method = method,
    marker = rownames(effects), trait = colnames(effects)[[1L]],
    true_effect = simulation$truth$effects[rownames(effects), 1L],
    posterior_mean_effect = effects[, 1L],
    posterior_inclusion_probability = probabilities[, 1L],
    causal = rownames(effects) %in% causal, stringsAsFactors = FALSE)
  out$causal_rank <- rank(-out$posterior_inclusion_probability,
    ties.method = "min")
  out
}

finemapping_credible_set_members <- function(native) {
  if (is.null(native$sets)) return(data.frame())
  rows <- list()
  for (locus in names(native$sets)) for (set_id in names(native$sets[[locus]])) {
    value <- native$sets[[locus]][[set_id]]
    members <- if (is.character(value)) value else if (!is.null(value$markers))
      value$markers else if (!is.null(value$marker)) value$marker else rownames(value)
    rows[[length(rows) + 1L]] <- data.frame(locus = locus, set = set_id,
      marker = as.character(members), member_order = seq_along(members),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

evaluate_finemapping_credible_sets <- function(native, simulation, positions,
    ld, method, min_r2 = 0.5, max_locus_distance = 1e6) {
  members <- finemapping_credible_set_members(native)
  if (!nrow(members)) return(data.frame())
  causals <- simulation$truth$causal$all
  causal_positions <- positions[match(causals, names(positions))]
  keys <- unique(members[c("locus", "set")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    z <- members[members$locus == keys$locus[[i]] &
      members$set == keys$set[[i]], , drop = FALSE]
    locus_row <- native$loci[match(keys$locus[[i]], native$loci$locus),,
      drop = FALSE]
    center <- if (nrow(locus_row)) mean(c(locus_row$start, locus_row$end)) else
      mean(positions[z$marker], na.rm = TRUE)
    eligible <- which(abs(causal_positions - center) <= max_locus_distance)
    causal_i <- if (length(eligible)) eligible[order(
      abs(causal_positions[eligible] - center), eligible)][[1L]] else
      which.min(abs(causal_positions - center))
    causal <- causals[[causal_i]]
    exact <- causal %in% z$marker
    usable <- intersect(z$marker, colnames(ld))
    proxy <- exact || (causal %in% rownames(ld) && length(usable) &&
      any(ld[causal, usable]^2 >= min_r2))
    summary_row <- native$summary[native$summary$locus == keys$locus[[i]] &
      native$summary$set_id == keys$set[[i]], , drop = FALSE]
    cumulative <- if (nrow(summary_row) &&
      "cumulative_pip" %in% names(summary_row))
      summary_row$cumulative_pip[[1L]] else NA_real_
    data.frame(study = simulation$scenario$study,
      scenario = simulation$scenario$architecture,
      replicate = simulation$scenario$replicate, method = method,
      trait = simulation$data$trait_names[[1L]], locus = causal,
      set = paste(keys$locus[[i]], keys$set[[i]], sep = ":"),
      set_size = nrow(z), cumulative_pip = cumulative,
      exact_covered = exact, ld_proxy_covered = proxy,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

summarise_finemapping_metrics <- function(marker_results, credible_sets) {
  marker <- evaluate_metrics_from_tables(marker_results)
  if (is.null(credible_sets) || !nrow(credible_sets)) return(marker)
  groups <- split(credible_sets, interaction(credible_sets$scenario,
    credible_sets$replicate, credible_sets$method, drop = TRUE))
  cs <- do.call(rbind, lapply(groups, function(x) {
    values <- c(number_of_credible_sets = nrow(x),
      mean_set_size = mean(x$set_size),
      median_set_size = stats::median(x$set_size),
      exact_coverage_fraction = mean(x$exact_covered),
      ld_proxy_coverage_fraction = mean(x$ld_proxy_covered),
      causal_loci_detected = length(unique(x$locus[x$exact_covered])),
      sets_without_exact_causal = sum(!x$exact_covered),
      sets_without_causal_or_proxy = sum(!x$ld_proxy_covered))
    data.frame(study = x$study[[1L]], scenario = x$scenario[[1L]],
      replicate = x$replicate[[1L]], method = x$method[[1L]],
      metric = names(values), value = unname(values), status = "ok",
      reason = "", stringsAsFactors = FALSE)
  }))
  rbind(marker, cs)
}

evaluate_metrics_from_tables <- function(marker_results) {
  groups <- split(marker_results, interaction(marker_results$scenario,
    marker_results$replicate, marker_results$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(x) {
    truth <- as.numeric(x$causal)
    pip <- x$posterior_inclusion_probability
    effect <- x$posterior_mean_effect
    ordering <- order(-pip, seq_along(pip))
    precision <- cumsum(truth[ordering]) / seq_along(ordering)
    causal_ranks <- sort(x$causal_rank[x$causal])
    values <- c(pip_brier = mean((pip - truth)^2),
      effect_rmse = sqrt(mean((effect - x$true_effect)^2)),
      average_precision = mean(precision[truth[ordering] == 1]),
      causal_top_10_recall = mean(causal_ranks <= 10))
    data.frame(study=x$study[[1L]], scenario=x$scenario[[1L]],
      replicate=x$replicate[[1L]], method=x$method[[1L]],
      metric=names(values),value=unname(values),status="ok",reason="",
      stringsAsFactors=FALSE)
  }))
}
