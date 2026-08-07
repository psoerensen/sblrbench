#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
fit_arg <- sub("^--fit=", "", args[grepl("^--fit=", args)])
validate_only <- "--validate-only" %in% args
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
library_path <- file.path(root, "results", "local", "06_annotation_models",
  "gctb_block_contract_validation", "rlib")
if (!dir.exists(library_path)) stop("The pinned isolated sblr library is absent.")
.libPaths(c(library_path, .libPaths()))
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "studies", "06_annotation_models", "power-isolation.R"))
source(file.path(root, "studies", "06_annotation_models",
  "gctb-single-trajectory.R"))
source(file.path(root, "studies", "06_annotation_models",
  "gctb-block-contract-validation.R"))

cfg <- study06_gctb_block_constants(root)
registry <- study06_gctb_block_registry(cfg)
dir.create(cfg$output, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name) utils::write.csv(x,
  file.path(cfg$output, name), row.names = FALSE, na = "")
write_json <- function(x, name) jsonlite::write_json(x,
  file.path(cfg$output, name), auto_unbox = TRUE, pretty = TRUE,
  null = "null", digits = 16)

source_sha <- system2("git", c("-C", "../sblr", "rev-parse", "HEAD"),
  stdout = TRUE)
source_status <- system2("git", c("-C", "../sblr", "status", "--short"),
  stdout = TRUE)
if (!identical(source_sha, cfg$sblr_sha) || length(source_status))
  stop("Pinned sibling source identity or cleanliness failed.")
provenance <- benchmark_package_provenance("sblr", lib.loc = library_path)
if (!identical(provenance$version, cfg$sblr_version) ||
    !identical(provenance$installed_tree_sha256,
      cfg$installed_tree_sha256)) stop("Installed sblr identity failed.")

official_cfg <- study06_gctb_single_constants(root)
export <- study06_gctb_validate_export(official_cfg)
for (id in c("D0", "D1")) {
  record <- jsonlite::read_json(file.path(cfg$official_output, "runs", id,
    "run-record.json"), simplifyVector = TRUE)
  if (!identical(record$status, "complete") ||
      !identical(record$official_sha, cfg$official_sha))
    stop("Pinned official output is missing or has changed: ", id)
}

glist_cache <- file.path(root, "results", "local", "06_annotation_models",
  "checkpoints", "data", "human_glist.rds")
if (!file.exists(glist_cache)) stop("Small Study 06 genotype cache is absent.")
Sys.setenv(SBLR_BENCH_GLIST = glist_cache, OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLAS_NUM_THREADS = "1")
spec <- read_benchmark_spec(file.path(root, "studies",
  "06_annotation_models", "spec.R"))
data <- prepare_prediction_data(spec, cfg$output)
logic <- .annotation_logic(spec)
annotations <- logic$construct_annotation_design(data$markers$marker_ids, spec)
annotation_truth <- logic$construct_annotation_truth(annotations, spec)
seed_grid <- benchmark_annotation_seeds(spec, "benchmark", "qualification")
seed_row <- seed_grid[seed_grid$scenario == "informative_annotations" &
  seed_grid$replicate == 1L & seed_grid$method == "st_bed_bayesrc", ,
  drop = FALSE]
coordinate <- as.list(seed_row[1L, c("scenario", "replicate",
  "component_seed", "effect_seed", "residual_seed")])
simulation <- logic$simulate_annotation_architecture(coordinate,
  data$scaled$all, data$split$train_rows, annotations, annotation_truth, spec)
test_simulation <- subset_sblrbench_simulation_samples(simulation,
  data$split$test_ids)
stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
  spec$data)
bundle <- list(simulation = simulation, test_simulation = test_simulation,
  stats = stats, annotations = annotations, annotation_truth = annotation_truth,
  marker_truth = simulation$extras$marker_truth)
truth_identity <- study06_power_truth_identity(data, bundle)
truth_hash <- benchmark_hash_object(truth_identity)
reference_identity <- readRDS(file.path(root, "results", "local",
  "06_annotation_models", "v2_paired_power_isolation", "checkpoints",
  "baseline--block_eigen.rds"))$identity
identity_checks <- c(
  specification = identical(benchmark_annotation_spec_hash(spec),
    cfg$specification_hash),
  truth = identical(reference_identity$truth_hash, cfg$truth_hash) &&
    identical(export$manifest$truth_hash, cfg$truth_hash),
  marker_order = identical(data$markers$marker_ids,
    export$truth$marker_ids), markers = nrow(annotations) == 1500L,
  blocks = length(data$block_start) == 15L,
  effects = isTRUE(all.equal(as.numeric(simulation$truth$effects),
    export$truth$effects, tolerance = 0)),
  phenotype = isTRUE(all.equal(as.numeric(simulation$truth$phenotypes[
    data$split$train_ids, 1L]), export$truth$training_y, tolerance = 0)))
if (!all(identity_checks)) stop("Small comparison immutable identity failed: ",
  paste(names(identity_checks)[!identity_checks], collapse = ", "), ".")

input_audit <- list(specification_hash = cfg$specification_hash,
  truth_hash = cfg$truth_hash, marker_count = nrow(annotations),
  block_count = length(data$block_start), block_sizes = as.list(table(
    export$truth$block)), all_positive_modes_retained = all(
      export$ld_audit$rank == 100L & export$ld_audit$reader_k == 100L),
  gamma = as.list(cfg$matched_gamma %||% c(0, .01, .1, 1)),
  package = provenance, sibling_sha = source_sha,
  official_sha = cfg$official_sha, gctb_sha = cfg$gctb_sha)
write_json(input_audit, "input_audit.json")
write_csv(study06_residual_semantic_crosswalk(), "semantic_crosswalk.csv")
write_csv(registry, "registry.csv")
if (validate_only) {
  message("Validated pinned small contract registry; no fit was run.")
  quit(save = "no", status = 0L)
}

methods <- resolve_benchmark_methods(spec)
names(methods) <- vapply(methods, `[[`, character(1), "id")
ids <- if (length(fit_arg)) fit_arg else registry$fit_id
if (any(!ids %in% registry$fit_id)) stop("Unknown --fit value.")
dir.create(file.path(cfg$output, "checkpoints"), recursive = TRUE,
  showWarnings = FALSE)

for (id in ids) {
  row <- registry[registry$fit_id == id, , drop = FALSE]
  controls <- study06_gctb_block_controls(spec, row, annotation_truth)
  identity <- list(schema = cfg$schema, fit = as.list(row),
    specification_hash = cfg$specification_hash, truth_hash = cfg$truth_hash,
    sblr_sha = cfg$sblr_sha,
    installed_tree_sha256 = cfg$installed_tree_sha256,
    residual_policy = cfg$residual_policy, block_ve_mode = cfg$block_ve_mode,
    resam_thresh = cfg$resam_thresh,
    minimum_ve_ratio = cfg$minimum_ve_ratio,
    niter = cfg$niter, package_nit = cfg$retained,
    burn = cfg$burn, retained = cfg$retained)
  semantic_hash <- benchmark_hash_object(identity)
  checkpoint <- file.path(cfg$output, "checkpoints", paste0(id, ".rds"))
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    if (!identical(saved$semantic_hash, semantic_hash))
      stop("Checkpoint identity mismatch: ", id)
    message("Reusing identity-matched ", id)
    next
  }
  started <- Sys.time()
  result <- study06_gctb_block_fit(methods[[row$method_id]], controls,
    simulation, stats, data$ld_glist, data$split, annotations,
    annotation_truth, data$block_start)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  native <- result$native_fit
  required_policy <- identical(native$input$residual_policy,
      cfg$residual_policy) && identical(native$input$block_ve_mode,
      cfg$block_ve_mode) && identical(native$input$resam_thresh,
      cfg$resam_thresh) && identical(native$input$minimum_ve_ratio,
      cfg$minimum_ve_ratio)
  if (!required_policy || is.null(native$block_ve) ||
      is.null(native$heritability_summary))
    stop("Resolved block residual contract or output is incomplete: ", id)
  benchmark_atomic_save_rds(list(identity = identity,
    semantic_hash = semantic_hash, elapsed_seconds = elapsed, result = result),
    checkpoint, compress = FALSE, temporary_prefix = ".gctb-block-")
  message("Completed ", id, " in ", round(elapsed, 2), " seconds")
  rm(result); invisible(gc())
}

if (!all(file.exists(file.path(cfg$output, "checkpoints",
    paste0(registry$fit_id, ".rds"))))) quit(save = "no", status = 0L)

metric_rows <- architecture_rows <- alpha_rows <- sigma_rows <- list()
component_rows <- stick_rows <- prior_rows <- list()
runtime_rows <- list()
for (i in seq_len(nrow(registry))) {
  row <- registry[i, , drop = FALSE]
  saved <- readRDS(file.path(cfg$output, "checkpoints",
    paste0(row$fit_id, ".rds")))
  s <- study06_gctb_block_extract(saved$result, export$truth$marker_ids,
    export$truth)
  o <- study06_gctb_block_official(row$official_id, cfg, export$truth)
  metric_rows[[i]] <- study06_gctb_block_compare(s, o, export$truth,
    row$fit_id)
  architecture_rows[[i]] <- data.frame(fit_id = row$fit_id,
    source = rep(c("sblr", "official"), each = 1L),
    active_count = c(s$active_count, o$active_count),
    mean_block_ve = c(s$mean_block_ve, o$mean_block_ve),
    median_block_ve = c(s$median_block_ve, o$median_block_ve),
    block_ve_min = c(s$block_ve_min, o$block_ve_min),
    block_ve_max = c(s$block_ve_max, o$block_ve_max),
    heritability = c(s$heritability, o$heritability),
    resamples = c(s$resamples, NA_real_), resets = c(s$resets, NA_real_))
  component_rows[[i]] <- rbind(
    data.frame(fit_id = row$fit_id, source = "sblr",
      component = seq_along(s$component_count) - 1L,
      mean_count = as.numeric(s$component_count)),
    data.frame(fit_id = row$fit_id, source = "official",
      component = seq_along(o$component_count) - 1L,
      mean_count = as.numeric(o$component_count)))
  tr <- saved$result$native_fit$chains[[1L]]$convergence_trace
  stick_rows[[i]] <- data.frame(fit_id = row$fit_id,
    stick = seq_len(ncol(tr$stick_eligible_count)),
    eligible = colMeans(tr$stick_eligible_count),
    continuation = colMeans(tr$stick_continue_count),
    stopping = colMeans(tr$stick_stop_count))
  runtime_rows[[i]] <- data.frame(fit_id = row$fit_id,
    elapsed_seconds = saved$elapsed_seconds)
  if (row$annotation_aware) {
    trace <- saved$result$native_fit$chains[[1L]]$convergence_trace
    if (!identical(dim(trace$alpha), c(cfg$retained, 12L)) ||
        !identical(dim(trace$sigmaSqAlpha), c(cfg$retained, 3L)))
      stop("S1-new retained alpha/sigmaSqAlpha shapes are invalid.")
    annotations_order <- colnames(annotations)
    alpha_map <- expand.grid(annotation = annotations_order,
      stick = paste0("stick_", 1:3), KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE)
    alpha_rows[[i]] <- do.call(rbind, lapply(seq_len(12L), function(j) {
      x <- trace$alpha[, j]
      data.frame(parameter = "alpha", annotation = alpha_map$annotation[j],
        stick = alpha_map$stick[j], mean = mean(x), sd = stats::sd(x),
        lower = unname(stats::quantile(x, .025)), median = stats::median(x),
        upper = unname(stats::quantile(x, .975)))
    }))
    sigma_rows[[i]] <- do.call(rbind, lapply(seq_len(3L), function(j) {
      x <- trace$sigmaSqAlpha[, j]
      data.frame(quantity = paste0("sigmaSqAlpha_stick_", j),
        mean = mean(x), sd = stats::sd(x),
        lower = unname(stats::quantile(x, .025)),
        median = stats::median(x),
        upper = unname(stats::quantile(x, .975)))
    }))
    alpha_mean <- saved$result$native_fit$annotation_effects[[1L]]
    q <- stats::pnorm(annotations %*% alpha_mean)
    pi <- cbind(1 - q[, 1L], q[, 1L] * (1 - q[, 2L]),
      q[, 1L] * q[, 2L] * (1 - q[, 3L]),
      q[, 1L] * q[, 2L] * q[, 3L])
    active_prior <- 1 - pi[, 1L]
    enriched <- annotations[, "enriched_binary"] == 1
    prior_rows[[i]] <- rbind(
      data.frame(fit_id = row$fit_id, summary = paste0("mean_pi_", 0:3),
        value = colMeans(pi)),
      data.frame(fit_id = row$fit_id,
        summary = c("expected_active", "active_prior_enriched",
          "active_prior_unenriched"),
        value = c(sum(active_prior), mean(active_prior[enriched]),
          mean(active_prior[!enriched]))))
  }
}
metrics <- do.call(rbind, metric_rows)
gate <- study06_gctb_block_gate(metrics, cfg)
write_csv(metrics, "comparison_metrics.csv")
write_csv(do.call(rbind, architecture_rows), "architecture.csv")
write_csv(do.call(rbind, component_rows), "component_occupancy.csv")
write_csv(do.call(rbind, stick_rows), "stick_summary.csv")
write_csv(do.call(rbind, runtime_rows), "runtime.csv")
if (length(alpha_rows)) write_csv(do.call(rbind, alpha_rows),
  "alpha_summary.csv")
if (length(sigma_rows)) write_csv(do.call(rbind, sigma_rows),
  "sigmaSqAlpha_summary.csv")
if (length(prior_rows)) write_csv(do.call(rbind, prior_rows),
  "annotation_probability_summary.csv")
decision <- list(schema = cfg$schema, decision = gate$decision,
  pass_by_fit = as.list(gate$pass), proceed_to_large = all(gate$pass),
  specification_hash = cfg$specification_hash, truth_hash = cfg$truth_hash,
  sblr_sha = cfg$sblr_sha,
  installed_tree_sha256 = cfg$installed_tree_sha256,
  official_sha = cfg$official_sha, gctb_sha = cfg$gctb_sha,
  residual_policy = cfg$residual_policy, block_ve_mode = cfg$block_ve_mode,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
write_json(decision, "decision.json")
message("Stage A decision: ", gate$decision)
