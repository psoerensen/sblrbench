.study05_paths <- function() list(
  glist_path = Sys.getenv("SBLR_BENCH_GLIST", ""),
  data_dir = Sys.getenv("SBLR_BENCH_DATA_DIR",
    file.path("results", "local", "02_prediction", "data")),
  genotype_output_dir = Sys.getenv("SBLR_BENCH_LD_DIR",
    file.path("results", "local", "02_prediction", "genotype_setup")),
  local_dir = Sys.getenv("SBLR_BENCH_STUDY05_LOCAL",
    file.path("results", "local", "study05_annotation_models")))

.study05_load_scaled_genotypes <- function(base_glist, chr, sample_ids,
                                           marker_ids, split) {
  prior_store <- file.path("results", "local", "five_replicate_overnight",
    "_targets_02")
  if (dir.exists(prior_store)) {
    cached <- try(targets::tar_read(prediction_scaled_genotypes,
      store = prior_store), silent = TRUE)
    if (!inherits(cached, "try-error") &&
        identical(rownames(cached$all), sample_ids) &&
        identical(colnames(cached$all), marker_ids) &&
        identical(rownames(cached$train), split$train_ids) &&
        identical(rownames(cached$test), split$test_ids))
      return(cached)
  }
  raw <- sblrbench:::benchmark_extract_raw_genotypes(base_glist, chr,
    sample_ids, marker_ids)
  sblrbench::training_scaled_genotypes(raw, split$train_rows)
}

.study05_summary_stats <- function(simulation, Glist, split, config) {
  started <- proc.time()[["elapsed"]]
  y <- simulation$phenotype[split$train_ids, , drop = FALSE]
  stats <- sblr::make_summary_stats(Glist = Glist, y = y, chr = config$chr,
    rows = split$train_rows, scale = TRUE, nthreads = 1L)
  elapsed <- proc.time()[["elapsed"]] - started
  if (!identical(stats$marker_names, simulation$marker_truth$marker_id) ||
      !identical(stats$trait_names, config$trait) ||
      !identical(as.integer(stats$n), length(split$train_ids)) ||
      !isTRUE(all.equal(unname(stats$af[[1L]]),
        unname(Glist$sparseLD$af[[1L]]), tolerance = 0)))
    stop("Study 05 training summary-statistic contract failed.", call. = FALSE)
  list(stats = stats, elapsed_seconds = elapsed)
}

.study05_specs <- function(config) {
  out <- list()
  for (scenario in config$scenarios)
    for (replicate in seq_len(config$simulation$replicate_count))
      out[[length(out) + 1L]] <- list(scenario = scenario, replicate = replicate)
  out
}

.study05_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- x[do.call(order, unname(x)), , drop = FALSE]
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

.study05_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null")
  path
}

.study05_fit_checkpoint <- function(run, phase) {
  root <- Sys.getenv("SBLR_BENCH_STUDY05_LOCAL",
    file.path("results", "local", "study05_annotation_models"))
  dir <- file.path(root, "fit_checkpoints", phase)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- paste0(file.path(dir, paste(run$scenario, run$replicate,
    run$method$id, sep = "__")), ".csv")
  x <- data.frame(phase = phase, scenario = run$scenario,
    replicate = run$replicate, method = run$method$id,
    status = run$status, started_at = format(run$started_at, tz = "UTC",
      usetz = TRUE), finished_at = format(run$finished_at, tz = "UTC",
      usetz = TRUE), elapsed_seconds = run$runtime,
    fit_seed = run$seeds[["fit_seed"]],
    chain_seeds = paste(run$seeds[paste0("chain_", 1:4)], collapse = ";"),
    effective_chain_seeds = if (identical(run$status, "ok"))
      paste(vapply(run$fit$chains, `[[`, integer(1), "seed"),
        collapse = ";") else "",
    chain_context = "native four-chain fit",
    error_message = run$reason, stringsAsFactors = FALSE)
  tmp <- tempfile(".fit-", dir, ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  backup <- paste0(path, ".previous")
  if (file.exists(backup)) unlink(backup)
  if (file.exists(path) && !file.rename(path, backup))
    stop("Could not preserve prior fit checkpoint: ", path, call. = FALSE)
  if (!file.rename(tmp, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    unlink(tmp)
    stop("Atomic fit checkpoint replacement failed: ", path, call. = FALSE)
  }
  if (file.exists(backup)) unlink(backup)
  if (!identical(run$status, "ok"))
    stop("Study 05 fit failed; structured checkpoint written to ",
      path, ": ", run$reason, call. = FALSE)
  path
}

.study05_seed_registry <- function(config) {
  methods <- config$methods
  rows <- list()
  for (scenario in config$scenarios) for (replicate in seq_len(5L)) {
    sim <- .study05_simulation_seeds(scenario, replicate, config)
    for (method in methods) {
      fit <- .study05_fit_seeds(scenario, replicate, method, config)
      for (chain in seq_len(4L))
        rows[[length(rows) + 1L]] <- data.frame(
          scenario = scenario, replicate = replicate, method = method,
          data_selection_seed = config$seeds$data_selection,
          annotation_seed = config$seeds$annotation,
          component_allocation_seed = sim[["component_allocation"]],
          effect_generation_seed = sim[["effect_generation"]],
          phenotype_residual_seed = sim[["phenotype_residuals"]],
          fit_seed = fit[["fit_seed"]], chain = chain,
          chain_seed = fit[[paste0("chain_", chain)]],
          stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

.study05_effective_seed_registry <- function(runs, config) {
  registry <- .study05_seed_registry(config)
  registry$effective_chain_seed <- NA_integer_
  for (run in runs) {
    if (!identical(run$status, "ok") || length(run$fit$chains) != 4L) next
    effective <- vapply(run$fit$chains, `[[`, integer(1), "seed")
    for (chain in seq_len(4L)) {
      i <- registry$scenario == run$scenario &
        registry$replicate == run$replicate &
        registry$method == run$method$id & registry$chain == chain
      registry$effective_chain_seed[i] <- effective[chain]
    }
  }
  registry
}
