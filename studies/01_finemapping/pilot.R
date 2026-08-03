`%||%` <- function(x, y) if (is.null(x)) y else x

.study01_assert_ids <- function(x, label) {
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop(label, " must be unique, non-missing character identifiers.", call. = FALSE)
  x
}

select_separated_causal_markers <- function(Glist, chr, marker_ids, n_causal,
                                             min_distance_bp, seed,
                                             min_maf = NULL, max_maf = NULL) {
  marker_ids <- .study01_assert_ids(marker_ids, "marker_ids")
  chr <- as.integer(chr); n_causal <- as.integer(n_causal); seed <- as.integer(seed)
  if (length(chr) != 1L || is.na(chr) || chr < 1L || length(n_causal) != 1L ||
      is.na(n_causal) || n_causal < 1L || !is.finite(min_distance_bp) || min_distance_bp < 0 ||
      length(seed) != 1L || is.na(seed)) stop("Invalid causal-selection scalar setting.", call. = FALSE)
  chromosome_ids <- .study01_assert_ids(Glist$rsids[[chr]], "chromosome marker IDs")
  absent <- marker_ids[is.na(match(marker_ids, chromosome_ids))]
  if (length(absent)) stop("Candidate markers are absent from the selected chromosome: ", paste(absent, collapse = ", "), call. = FALSE)
  idx <- match(marker_ids, chromosome_ids)
  pos <- Glist$pos[[chr]][idx]; maf <- Glist$maf[[chr]][idx]
  if (length(pos) != length(marker_ids) || any(!is.finite(pos))) stop("Aligned marker positions must be present and finite.", call. = FALSE)
  if (length(maf) != length(marker_ids) || any(!is.finite(maf)) || any(maf < 0 | maf > 0.5)) stop("Aligned MAF values must be finite and lie in [0, 0.5].", call. = FALSE)
  keep <- rep(TRUE, length(marker_ids))
  if (!is.null(min_maf)) keep <- keep & maf >= min_maf
  if (!is.null(max_maf)) keep <- keep & maf <= max_maf
  candidates <- which(keep)
  if (length(candidates) < n_causal) stop("Insufficient eligible candidates for the requested causal count.", call. = FALSE)
  set.seed(seed); starts <- sample(candidates)
  selected <- integer()
  for (i in starts) {
    if (!length(selected) || all(abs(pos[i] - pos[selected]) >= min_distance_bp)) selected <- c(selected, i)
    if (length(selected) == n_causal) break
  }
  if (length(selected) != n_causal) stop("Insufficient separated candidates for the requested causal count and distance.", call. = FALSE)
  selected <- sort(selected)
  minimum <- if (length(selected) < 2L) Inf else min(abs(outer(pos[selected], pos[selected], "-")[lower.tri(matrix(0, length(selected), length(selected)))]))
  list(marker_ids = marker_ids[selected], marker_indices = idx[selected], chromosome = chr,
       positions = pos[selected], maf = maf[selected], pairwise_min_distance_bp = minimum,
       selection_seed = seed, settings = list(n_causal = n_causal,
       min_distance_bp = min_distance_bp, min_maf = min_maf, max_maf = max_maf))
}

.study01_replicate_specs <- function(config, override = Sys.getenv("SBLR_BENCH_REPLICATES", "")) {
  n <- if (nzchar(override)) suppressWarnings(as.integer(override)) else config$replicate_counts[["development"]]
  if (length(n) != 1L || is.na(n) || !n %in% unname(config$replicate_counts)) stop("SBLR_BENCH_REPLICATES must be 1, 5, or 10.", call. = FALSE)
  lapply(seq_len(n), function(i) list(replicate = i, causal_seed = config$simulation$base_seed + i,
    simulation_seed = config$simulation$base_seed + 1000L + i))
}

.study01_method_specs <- function(config) {
  map <- list(st_bed_bayesc = c("stblr_bed", "bayesc"), st_bed_bayesr = c("stblr_bed", "bayesr"),
              st_csr_sbayesc = c("stblr_csr", "sbayesc"), st_csr_sbayesr = c("stblr_csr", "sbayesr"))
  lapply(seq_along(config$methods), function(i) {
    id <- config$methods[[i]]; x <- map[[id]]
    list(id = id, interface = x[[1]], native_method = x[[2]], method_index = i,
         capabilities = c("posterior_effects", "pip", "scalar_trait", if (x[[1]] == "stblr_bed") "individual_level" else "summary_statistics"))
  })
}

.study01_simulate <- function(spec, selection, Z, config) {
  raw <- sblr::mtsim(W = Z, standardize_W = FALSE, nt = 1L,
    n_shared = length(selection$marker_ids), n_specific = 0L,
    causal_rsids = selection$marker_ids, h2 = config$simulation$h2,
    seed = spec$simulation_seed)
  if (!identical(sort(raw$causal$all), sort(selection$marker_ids))) stop("mtsim() did not use the exact selected causal set.", call. = FALSE)
  sim <- sblrbench::as_sblrbench_simulation(raw, study = config$study,
    architecture = config$architecture, replicate = as.integer(spec$replicate), seed = spec$simulation_seed)
  # mtsim() 0.1.2 returns a vector for a one-trait phenotype. Materialize the
  # contract's sample-by-trait matrix without changing values or scale.
  if (is.null(dim(sim$truth$phenotypes))) sim$truth$phenotypes <- matrix(sim$truth$phenotypes,
    ncol = 1L, dimnames = list(sim$data$sample_ids, sim$data$trait_names))
  sim$extras$causal_selection <- selection
  sblrbench::validate_sblrbench_simulation(sim)
  sim
}

.study01_summary_stats <- function(simulation, Glist, config) {
  y <- simulation$truth$phenotypes
  if (!identical(rownames(y), Glist$ids)) stop("Phenotype sample names do not match Glist$ids.", call. = FALSE)
  stats <- sblr::make_summary_stats(Glist = Glist, y = y, chr = config$chr, scale = TRUE, nthreads = 1L)
  ids <- stats$marker_names %||% stats$rsids %||% rownames(stats$bhat) %||% rownames(stats$B)
  if (!identical(as.character(ids), simulation$data$marker_ids)) stop("Summary-statistic marker order does not match canonical marker order.", call. = FALSE)
  stats
}

.study01_fit <- function(method, simulation, stats, Glist, config) {
  seed <- as.integer(config$mcmc$seed_offset + simulation$scenario$replicate * 100L + method$method_index)
  controls <- config$mcmc[c("nit", "nburn", "nthin", "nchains", "ncores", "convergence_control")]
  controls$seed <- seed; controls$convergence <- "core"; controls$verbose <- FALSE
  spec <- sblrbench::new_sblr_native_method(method$id, method$id, method$interface,
    method$native_method, method$capabilities, metadata = list(development_settings = TRUE))
  fit_inputs <- if (method$interface == "stblr_bed") list(y = simulation$truth$phenotypes, Glist = Glist) else list(stats = stats, Glist = Glist)
  tryCatch({
    result <- sblrbench::run_sblrbench_method(spec, fit_inputs = fit_inputs, controls = controls)
    sblrbench::validate_sblrbench_result(result, simulation)
    consistency <- sblr::check_stblr_consistency(result$native_fit, verbose = FALSE)
    list(status = "ok", reason = "", method = method, mcmc_seed = seed, controls = controls,
         result = result, consistency = consistency)
  }, error = function(e) list(status = "failed", reason = conditionMessage(e), method = method,
      mcmc_seed = seed, controls = controls, result = NULL, consistency = NULL))
}

.study01_marker_metrics <- function(fit_record, simulation) {
  if (fit_record$status != "ok") return(data.frame(study = simulation$scenario$study, scenario = simulation$scenario$architecture,
    replicate = simulation$scenario$replicate, method_id = fit_record$method$id, trait = simulation$data$trait_names,
    metric = "method_fit", value = NA_real_, status = "failed", reason = fit_record$reason, stringsAsFactors = FALSE))
  sblrbench::evaluate_metrics(simulation, fit_record$result,
    metrics = c("pip_brier", "effect_rmse", "average_precision", "causal_ranks"))
}

.study01_credible_sets <- function(fit_record, Glist, config) {
  if (fit_record$status != "ok") return(list(status = "failed", reason = fit_record$reason, native = NULL))
  cs <- config$credible_sets
  tryCatch(list(status = "ok", reason = "", native = sblr::make_credible_sets(
    fit = fit_record$result$native_fit, Glist = Glist, trait = 1L, coverage = cs$coverage,
    min_r2 = cs$min_r2, pip_cutoff = cs$pip_cutoff, locus_pip_cutoff = cs$locus_pip_cutoff,
    max_locus_distance = cs$max_locus_distance)),
    error = function(e) list(status = "failed", reason = conditionMessage(e), native = NULL))
}

.study01_cs_members <- function(cs) {
  if (is.null(cs$sets)) return(list())
  out <- list()
  for (locus in names(cs$sets)) for (set_id in names(cs$sets[[locus]])) {
    z <- cs$sets[[locus]][[set_id]]
    members <- if (is.character(z)) z else z$markers %||% z$marker %||% rownames(z) %||% character()
    out[[length(out) + 1L]] <- list(locus = locus, set_id = set_id, members = as.character(members))
  }
  out
}

evaluate_credible_sets <- function(credible_sets, simulation, positions, LD, method_id,
                                   min_r2 = 0.5, max_locus_distance = 1e6) {
  causals <- simulation$truth$causal$all
  traits <- simulation$data$trait_names
  empty <- data.frame(study = simulation$scenario$study, scenario = simulation$scenario$architecture,
    replicate = simulation$scenario$replicate, method_id = method_id, trait = traits[[1]], locus_id = NA_character_,
    set_id = NA_character_, coverage_type = c("exact", "ld_proxy"), covered = NA, set_size = NA_integer_,
    cumulative_pip = NA_real_, status = credible_sets$status, reason = credible_sets$reason, stringsAsFactors = FALSE)
  if (credible_sets$status != "ok") return(empty)
  members <- .study01_cs_members(credible_sets$native)
  if (!length(members)) { empty$status <- "ok"; empty$reason <- "no credible sets"; return(empty) }
  pos <- positions[match(causals, names(positions))]
  if (any(!is.finite(pos))) stop("Causal positions are unavailable.", call. = FALSE)
  summary <- credible_sets$native$summary
  rows <- list()
  for (x in members) {
    locus_row <- credible_sets$native$loci[match(x$locus, credible_sets$native$loci$locus), , drop = FALSE]
    center <- if (nrow(locus_row)) mean(c(locus_row$start, locus_row$end)) else mean(positions[x$members], na.rm = TRUE)
    eligible <- which(abs(pos - center) <= max_locus_distance)
    causal_i <- if (length(eligible)) eligible[order(abs(pos[eligible] - center), eligible)][1L] else which.min(abs(pos - center))
    causal <- causals[[causal_i]]
    exact <- causal %in% x$members
    proxy <- exact
    if (!proxy && causal %in% rownames(LD)) {
      usable <- intersect(x$members, colnames(LD)); if (length(usable)) proxy <- any(LD[causal, usable]^2 >= min_r2)
    }
    sr <- if (nrow(summary)) summary[summary$locus == x$locus & summary$set_id == x$set_id, , drop = FALSE] else summary
    cp <- if (nrow(sr) && "cumulative_pip" %in% names(sr)) sr$cumulative_pip[[1]] else NA_real_
    for (type in c("exact", "ld_proxy")) rows[[length(rows)+1L]] <- data.frame(study = simulation$scenario$study,
      scenario = simulation$scenario$architecture, replicate = simulation$scenario$replicate, method_id = method_id,
      trait = traits[[1]], locus_id = causal, set_id = paste(x$locus, x$set_id, sep = ":"), coverage_type = type,
      covered = if (type == "exact") exact else proxy, set_size = length(x$members), cumulative_pip = cp,
      status = "ok", reason = "", stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

.study01_compact_fit <- function(fit_record, simulation, config) {
  conv <- if (fit_record$status == "ok") fit_record$result$diagnostics$convergence else NULL
  data.frame(architecture = config$architecture, replicate = simulation$scenario$replicate,
    method = fit_record$method$id, trait = simulation$data$trait_names[[1]], simulation_seed = simulation$provenance$seed,
    causal_selection_seed = simulation$extras$causal_selection$selection_seed, mcmc_seed = fit_record$mcmc_seed,
    sblr_version = as.character(utils::packageVersion("sblr")), qgg_version = as.character(utils::packageVersion("qgg")),
    sblrbench_commit = sblrbench::sblrbench_git_commit("."), method_controls = jsonlite::toJSON(fit_record$controls, auto_unbox = TRUE),
    runtime = if (fit_record$status == "ok") fit_record$result$computation$elapsed_seconds else NA_real_,
    convergence_status = if (fit_record$status == "ok") if (is.null(conv)) "unavailable" else "reported" else "failed",
    status = fit_record$status, reason = fit_record$reason, stringsAsFactors = FALSE)
}

.study01_write_csv <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); utils::write.csv(x, path, row.names = FALSE); path }

.study01_metric_summary <- function(marker_metrics, credible_set_summary) {
  metric_rows <- rbind(
    marker_metrics[, c("method", "metric", "value")],
    credible_set_summary[, c("method", "metric", "value")]
  )
  metric_rows <- metric_rows[is.finite(metric_rows$value), , drop = FALSE]
  split_rows <- split(metric_rows, interaction(metric_rows$method,
    metric_rows$metric, drop = TRUE, lex.order = TRUE))
  out <- do.call(rbind, lapply(split_rows, function(x) data.frame(
    method = x$method[[1]], metric = x$metric[[1]], n = nrow(x),
    mean = mean(x$value), sd = if (nrow(x) > 1L) stats::sd(x$value) else NA_real_,
    median = stats::median(x$value), min = min(x$value), max = max(x$value),
    stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out[order(out$method, out$metric), , drop = FALSE]
}

.study01_evaluate_cs <- function(cs, fit_record, simulation, Z, Glist, config) {
  causals <- simulation$truth$causal$all
  members <- unique(unlist(lapply(.study01_cs_members(cs$native), `[[`, "members"), use.names = FALSE))
  ids <- unique(c(causals, members)); ids <- ids[!is.na(match(ids, colnames(Z)))]
  LD <- if (length(ids)) stats::cor(Z[, ids, drop = FALSE]) else matrix(numeric(), 0L, 0L)
  idx <- match(simulation$data$marker_ids, Glist$rsids[[config$chr]])
  positions <- stats::setNames(Glist$pos[[config$chr]][idx], simulation$data$marker_ids)
  evaluate_credible_sets(cs, simulation, positions, LD, fit_record$method$id,
    config$credible_sets$min_r2, config$credible_sets$max_locus_distance)
}

.study01_cs_summary <- function(x) {
  keys <- unique(x[, c("replicate", "method_id", "trait")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    z <- merge(keys[i, , drop = FALSE], x, by = c("replicate", "method_id", "trait"))
    ok <- z$status == "ok" & !is.na(z$set_id)
    sets <- unique(z$set_id[ok]); exact <- z[ok & z$coverage_type == "exact", ]; proxy <- z[ok & z$coverage_type == "ld_proxy", ]
    vals <- c(number_of_credible_sets = length(sets), mean_set_size = if (nrow(exact)) mean(exact$set_size) else NA,
      median_set_size = if (nrow(exact)) stats::median(exact$set_size) else NA,
      exact_coverage_fraction = if (nrow(exact)) mean(exact$covered) else NA,
      ld_proxy_coverage_fraction = if (nrow(proxy)) mean(proxy$covered) else NA,
      causal_loci_detected = length(unique(exact$locus_id[exact$covered])),
      sets_without_exact_causal = if (nrow(exact)) sum(!exact$covered) else NA,
      sets_without_causal_or_proxy = if (nrow(proxy)) sum(!proxy$covered) else NA)
    data.frame(replicate = keys$replicate[[i]], method_id = keys$method_id[[i]],
      trait = keys$trait[[i]], metric = names(vals), value = unname(vals),
      stringsAsFactors = FALSE)
  }))
}
