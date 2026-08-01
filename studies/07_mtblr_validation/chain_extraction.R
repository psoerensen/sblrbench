.study07_extract_draws <- function(run, architecture, replicate) {
  if (!identical(run$status, "ok"))
    stop("Cannot extract a failed Study 07 fit.", call. = FALSE)
  fit <- run$fit; bundle <- fit$convergence_traces
  values <- bundle$values; q <- bundle$quantities
  if (length(dim(values)) != 3L || dim(values)[2L] != 4L ||
      nrow(q) != dim(values)[3L])
    stop("Study 07 convergence trace dimensions are invalid.", call. = FALSE)
  trait_names <- colnames(fit$bm)
  name_one <- function(i) {
    group <- q$group[[i]]
    t1 <- q$trait_index[[i]]
    t2 <- if ("trait2_index" %in% names(q)) q$trait2_index[[i]] else -1L
    pattern <- if ("pattern_name" %in% names(q)) q$pattern_name[[i]] else NA
    if (group %in% c("vbs", "vgs", "ves", "vle", "vld"))
      return(paste(group, trait_names[[t1]], sep = "__"))
    if (group %in% c("cov_b", "cov_g", "cov_e"))
      return(paste(group, trait_names[[t1]], trait_names[[t2]], sep = "__"))
    if (group %in% c("pattern_pi", "pi_active"))
      return(paste("state_probability", pattern, sep = "__"))
    paste0(group, "__", i)
  }
  rows <- lapply(seq_len(nrow(q)), function(i) {
    x <- values[, , i, drop = TRUE]
    data.frame(architecture = architecture, replicate = replicate,
      implementation = run$implementation$id, estimand = name_one(i),
      iteration = rep(seq_len(nrow(x)), times = ncol(x)),
      chain = rep(seq_len(ncol(x)), each = nrow(x)),
      value = as.vector(x), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  get_matrix <- function(name) {
    ids <- which(q$group == name)
    setNames(lapply(ids, function(i) values[, , i, drop = TRUE]),
      vapply(ids, name_one, ""))
  }
  vgs <- get_matrix("vgs"); ves <- get_matrix("ves")
  derived <- list()
  for (trait in trait_names) {
    g <- vgs[[paste("vgs", trait, sep = "__")]]
    e <- ves[[paste("ves", trait, sep = "__")]]
    if (!is.null(g) && !is.null(e))
      derived[[paste("heritability", trait, sep = "__")]] <- g / (g + e)
  }
  cg <- get_matrix("cov_g")
  if (length(cg) && length(vgs) == 2L) {
    derived$genetic_covariance <- cg[[1L]]
    derived$genetic_correlation <- cg[[1L]] /
      sqrt(vgs[[1L]] * vgs[[2L]])
  }
  if (length(derived)) out <- rbind(out, do.call(rbind,
    lapply(names(derived), function(name) {
      x <- derived[[name]]
      data.frame(architecture = architecture, replicate = replicate,
        implementation = run$implementation$id, estimand = name,
        iteration = rep(seq_len(nrow(x)), times = ncol(x)),
        chain = rep(seq_len(ncol(x)), each = nrow(x)),
        value = as.vector(x), stringsAsFactors = FALSE)
    })))
  if (any(!is.finite(out$value)) ||
      !identical(sort(unique(out$chain)), 1:4))
    stop("Study 07 retained draws are non-finite or chains are ambiguous.",
      call. = FALSE)
  out
}

.study07_truth_for_estimand <- function(estimand, simulation) {
  traits <- colnames(simulation$effects)
  if (estimand == paste("vgs", traits[[1L]], sep = "__"))
    return(simulation$truth$cov_g[1L, 1L])
  if (estimand == paste("vgs", traits[[2L]], sep = "__"))
    return(simulation$truth$cov_g[2L, 2L])
  if (estimand == paste("ves", traits[[1L]], sep = "__"))
    return(simulation$truth$cov_e[1L, 1L])
  if (estimand == paste("ves", traits[[2L]], sep = "__"))
    return(simulation$truth$cov_e[2L, 2L])
  if (estimand == paste("heritability", traits[[1L]], sep = "__"))
    return(simulation$truth$heritability[[1L]])
  if (estimand == paste("heritability", traits[[2L]], sep = "__"))
    return(simulation$truth$heritability[[2L]])
  if (estimand == "genetic_covariance") return(simulation$truth$cov_g[1L, 2L])
  if (estimand == "genetic_correlation")
    return(simulation$truth$genetic_correlation)
  if (startsWith(estimand, "state_probability__")) {
    labels <- c("neither", "trait1_only", "trait2_only", "both")
    label <- sub("^state_probability__", "", estimand)
    return(as.numeric(simulation$truth$state_probabilities[match(label,
      labels)]))
  }
  NA_real_
}

.study07_parameter_estimates <- function(draws, simulation) {
  groups <- split(draws, draws$estimand)
  do.call(rbind, lapply(groups, function(x) {
    truth <- .study07_truth_for_estimand(x$estimand[[1L]], simulation)
    lower <- unname(stats::quantile(x$value, .025))
    upper <- unname(stats::quantile(x$value, .975))
    data.frame(architecture = x$architecture[[1L]],
      replicate = x$replicate[[1L]],
      implementation = x$implementation[[1L]],
      estimand = x$estimand[[1L]], available = TRUE,
      unavailable_reason = "", posterior_mean = mean(x$value),
      posterior_sd = stats::sd(x$value),
      posterior_median = stats::median(x$value), lower_95 = lower,
      upper_95 = upper, truth = truth,
      bias = if (is.finite(truth)) mean(x$value) - truth else NA_real_,
      squared_error = if (is.finite(truth))
        (mean(x$value) - truth)^2 else NA_real_,
      interval_coverage = if (is.finite(truth))
        truth >= lower && truth <= upper else NA,
      chain_count = length(unique(x$chain)),
      draws_per_chain = unique(as.integer(table(x$chain)))[[1L]],
      stringsAsFactors = FALSE)
  }))
}
