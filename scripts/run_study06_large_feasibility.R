#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) default else args[[hit + 1L]]
}
stage <- value("--stage", "all")
fit_id <- value("--fit", "")
resume <- tolower(value("--resume", "true")) == "true"
phase <- value("--phase", "initial")
if (!phase %in% c("initial", "continuation"))
  stop("--phase must be initial or continuation.")
output <- value("--output-dir",
  "results/local/06_annotation_models/large_feasibility")
isolated <- normalizePath("results/local/current_benchmark_refresh/rlib",
  winslash = "/", mustWork = TRUE)
.libPaths(c(isolated, .libPaths()))

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required.")
pkgload::load_all(".", quiet = TRUE)
source("studies/06_annotation_models/large-feasibility.R")
spec <- study06_large_spec()
registry <- study06_large_registry(spec)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE,
    null = "null", na = "null", digits = 16)
}
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
sha_file <- function(path) digest::digest(file = path, algo = "sha256")
bundle_path <- file.path(output, "prepared_bundle.rds")

prepare <- function() {
  message("Preparing immutable large-feasibility truth; no MCMC is running.")
  glist_path <- spec$source$glist
  glist <- readRDS(glist_path)
  if (length(glist$ids) != spec$source$sample_count || anyDuplicated(glist$ids))
    stop("Canonical sample audit failed.")
  markers <- benchmark_filter_markers(glist, spec$source$chromosome,
    spec$source$qc)
  if (markers$marker_count_after < 37000L || markers$marker_count_after > 39000L)
    stop("Canonical QC marker count is outside the registered range.")
  idx <- match(markers$marker_ids, glist$rsids[[spec$source$chromosome]])
  if (anyNA(idx)) stop("QC marker alignment failed.")
  analysis_order <- order(markers$positions, idx, method = "radix")
  markers$marker_ids <- markers$marker_ids[analysis_order]
  markers$positions <- markers$positions[analysis_order]
  markers$af <- markers$af[analysis_order]
  markers$maf <- markers$maf[analysis_order]
  idx <- idx[analysis_order]
  markers$marker_order_id <- paste(spec$source$chromosome,
    length(markers$marker_ids), markers$marker_ids[[1L]],
    markers$marker_ids[[length(markers$marker_ids)]], sep = ":")
  sample_ids <- as.character(glist$ids)
  a1 <- glist$a1[[1L]][idx]; a2 <- glist$a2[[1L]][idx]
  nmiss <- glist$nmiss[[1L]][idx]
  if (anyDuplicated(markers$marker_ids) || anyNA(c(a1, a2)) ||
      any(!is.finite(markers$af)) || any(markers$maf <= 0))
    stop("Canonical marker/allele audit failed.")
  glist <- benchmark_set_glist_marker_order(glist, spec$source$chromosome,
    markers$marker_ids)
  blocks <- study06_large_blocks(markers$marker_ids, markers$positions, spec)
  annotation <- study06_large_annotations(markers$marker_ids,
    blocks$panel$block_id, spec)
  calibrated <- study06_large_calibrate_alpha(annotation$matrix, spec)
  truth <- study06_large_truth(annotation$matrix, calibrated, glist,
    markers$marker_ids, sample_ids, blocks$panel$block_id, spec)
  gates <- study06_large_truth_gates(truth, annotation$matrix, calibrated, spec)
  if (!gates$pass) stop("LARGE-F6 truth gate: ", paste(gates$reasons,
    collapse = "; "))
  trace_panel <- study06_large_trace_panel(truth$marker_truth,
    annotation$matrix, markers$positions, spec)
  message("Truth gates passed; auditing 76 full-rank block operators.")
  ld_audit <- study06_large_ld_audit(glist, sample_ids, blocks$panel, spec)
  gwas <- study06_large_gwas(glist, truth$phenotype, markers$marker_ids, spec)
  initial_marker_count <- length(glist$rsids[[spec$source$chromosome]])
  marker_audit <- data.frame(initial_individuals = length(sample_ids),
    final_individuals = length(sample_ids), initial_markers = initial_marker_count,
    final_markers = length(markers$marker_ids), missing_genotypes = sum(nmiss),
    missing_rate = sum(nmiss) / (length(sample_ids) * length(markers$marker_ids)),
    maf_min = min(markers$maf), maf_max = max(markers$maf),
    af_min = min(markers$af), af_max = max(markers$af),
    chromosome_count = length(unique(blocks$panel$chromosome)),
    spacing_min = min(diff(markers$positions)),
    spacing_median = median(diff(markers$positions)),
    spacing_max = max(diff(markers$positions)),
    genotype_variance_min = min(2 * markers$af * (1 - markers$af)),
    genotype_variance_max = max(2 * markers$af * (1 - markers$af)))
  identities <- list(
    specification_hash = study06_large_hash(spec),
    truth_hash = study06_large_hash(list(marker_truth = truth$marker_truth,
      genetic = truth$genetic_value, residual = truth$residual,
      phenotype = truth$phenotype)),
    marker_order_hash = study06_large_hash(markers$marker_ids),
    sample_order_hash = study06_large_hash(sample_ids),
    allele_hash = study06_large_hash(data.frame(marker_id = markers$marker_ids,
      a1 = a1, a2 = a2)), annotation_hash = annotation$audit$annotation_hash,
    block_hash = blocks$block_hash, gwas_hash = gwas$hash,
    alpha_truth_hash = calibrated$alpha_hash,
    trace_panel_hash = trace_panel$hash,
    genotype_source_hash = sha_file(glist$bedfiles[[1L]]))
  source_files <- unique(c(glist$bedfiles, glist$bimfiles, glist$famfiles))
  source_hashes <- data.frame(path = normalizePath(source_files, winslash = "/"),
    size = file.info(source_files)$size,
    sha256 = vapply(source_files, sha_file, character(1)))
  bundle <- list(spec = spec, registry = registry, glist = glist,
    sample_ids = sample_ids, markers = markers,
    alleles = data.frame(marker_id = markers$marker_ids, a1 = a1, a2 = a2),
    blocks = blocks, annotations = annotation$matrix,
    annotation_audit = annotation$audit, calibrated = calibrated, truth = truth,
    gates = gates, trace_panel = trace_panel, ld_audit = ld_audit,
    gwas = gwas, identities = identities, marker_audit = marker_audit,
    source_hashes = source_hashes)
  benchmark_atomic_save_rds(bundle, bundle_path, compress = FALSE,
    temporary_prefix = ".large-feasibility-")
  write_csv(marker_audit, file.path(output, "marker_sample_audit.csv"))
  write_csv(source_hashes, file.path(output, "genotype_source_hashes.csv"))
  write_csv(blocks$panel, file.path(output, "block_panel.csv"))
  write_csv(ld_audit$blocks, file.path(output, "block_ld_audit.csv"))
  write_json(ld_audit$summary, file.path(output, "block_ld_summary.json"))
  write_csv(as.data.frame(annotation$matrix), file.path(output,
    "annotations.csv"))
  write_json(annotation$audit, file.path(output, "annotation_audit.json"))
  write_csv(as.data.frame(calibrated$alpha), file.path(output, "alpha_truth.csv"))
  write_csv(truth$marker_truth, file.path(output, "marker_truth.csv"))
  write_csv(gates$stick, file.path(output, "truth_stick_audit.csv"))
  write_json(list(pass = gates$pass, gates = as.list(gates$gates),
    component_count = as.list(gates$component_count),
    expected_count = as.list(gates$expected_count),
    realized_h2 = truth$realized_h2, genetic_variance = truth$genetic_variance,
    residual_variance = truth$residual_variance),
    file.path(output, "truth_gate_summary.json"))
  write_csv(trace_panel$panel, file.path(output, "selected_trace_panel.csv"))
  write_csv(registry, file.path(output, "model_registry.csv"))
  write_json(identities, file.path(output, "prepared_identities.json"))
  message("Prepared bundle: ", bundle_path)
  bundle
}

load_bundle <- function() {
  if (!file.exists(bundle_path)) stop("Prepared bundle is absent; run --stage prepare.")
  x <- readRDS(bundle_path)
  if (!identical(x$identities$specification_hash, study06_large_hash(spec)))
    stop("Prepared specification identity mismatch.")
  x
}

run_one <- function(id, smoke = FALSE) {
  bundle <- load_bundle()
  row <- registry[registry$fit_id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown fit id: ", id)
  controls <- study06_large_controls(row, bundle$calibrated$alpha, spec, smoke)
  suffix <- if (smoke) "smoke" else "fit"
  checkpoint_name <- if (phase == "initial") paste0(id, "_", suffix, ".rds") else
    paste0("continuation_", id, "_", suffix, ".rds")
  path <- file.path(output, "checkpoints", checkpoint_name)
  identity <- list(schema = spec$schema, fit = as.list(row),
    specification_hash = bundle$identities$specification_hash,
    truth_hash = bundle$identities$truth_hash,
    marker_hash = bundle$identities$marker_order_hash,
    annotation_hash = bundle$identities$annotation_hash,
    block_hash = bundle$identities$block_hash, gwas_hash = bundle$identities$gwas_hash,
    trace_hash = bundle$identities$trace_panel_hash, controls = controls,
    smoke = smoke, phase = phase,
    sblr_sha = system2("git", c("-C", "../sblr", "rev-parse", "HEAD"),
      stdout = TRUE),
    sblr_diff_sha256 = digest::digest(system2("git",
      c("-C", "../sblr", "diff", "--binary"), stdout = TRUE),
      algo = "sha256"),
    installed_package = benchmark_package_provenance("sblr"))
  semantic_hash <- study06_large_hash(identity)
  if (resume && file.exists(path)) {
    saved <- readRDS(path)
    if (!identical(saved$semantic_hash, semantic_hash))
      stop("Checkpoint identity mismatch: ", id)
    message("Reused identity-matched ", suffix, " checkpoint for ", id)
    return(saved)
  }
  started <- Sys.time()
  fit <- study06_large_fit_call(row, controls, bundle, spec)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (any(!is.finite(fit$bm)) || any(!is.finite(fit$dm)))
    stop("Non-finite posterior output for ", id)
  saved <- list(identity = identity, semantic_hash = semantic_hash,
    elapsed_seconds = elapsed, fit = fit)
  benchmark_atomic_save_rds(saved, path, compress = FALSE,
    temporary_prefix = paste0(".large-", id, "-"))
  message("Completed ", id, " ", suffix, " in ", round(elapsed, 2), " seconds")
  saved
}

write_blocked_manifest <- function() {
  bundle <- load_bundle()
  smoke_paths <- file.path(output, "checkpoints",
    paste0(c("E0", "E2", "E1"), "_smoke.rds"))
  installed <- find.package("sblr")
  installed_files <- sort(list.files(installed, recursive = TRUE,
    full.names = TRUE))
  installed_hashes <- vapply(installed_files, sha_file, character(1))
  installed_tree_hash <- digest::digest(paste(
    sub(normalizePath(installed, winslash = "/"), "",
      normalizePath(installed_files, winslash = "/")), installed_hashes,
    sep = "="), algo = "sha256")
  manifest <- list(schema = "sblrbench-study06-large-feasibility-manifest-v1",
    decision = "LARGE-F6", identities = bundle$identities,
    sblrbench_sha = system2("git", "rev-parse HEAD", stdout = TRUE),
    sblr_sha = system2("git", c("-C", "../sblr", "rev-parse", "HEAD"),
      stdout = TRUE), sblr_version = as.character(packageVersion("sblr")),
    sblr_path = normalizePath(installed, winslash = "/"),
    installed_sblr_tree_sha256 = installed_tree_hash,
    source_hashes = bundle$source_hashes,
    model_configuration_hashes = stats::setNames(as.list(registry$config_hash),
      registry$fit_id),
    script_hashes = list(
      implementation = sha_file("studies/06_annotation_models/large-feasibility.R"),
      runner = sha_file("scripts/run_study06_large_feasibility.R")),
    raw_output_hashes = stats::setNames(as.list(vapply(smoke_paths, sha_file,
      character(1))), basename(smoke_paths)),
    smoke_runtime_seconds = list(E0 = readRDS(smoke_paths[1L])$elapsed_seconds,
      E2 = readRDS(smoke_paths[2L])$elapsed_seconds,
      E1 = readRDS(smoke_paths[3L])$elapsed_seconds),
    failed_smoke = list(fit_id = "B0",
      error = paste("BayesR operator residual scale is invalid.",
        "trait=0, chain=0, iter=0")),
    trace_storage = list(free_bytes_before_smoke = 54945726464,
      all_marker_component_r_array_bytes_per_fit =
        length(bundle$markers$marker_ids) * spec$mcmc$nit *
          spec$mcmc$nchains * 8,
      six_fit_bytes = length(bundle$markers$marker_ids) * spec$mcmc$nit *
        spec$mcmc$nchains * 8 * 6,
      compact_aggregate_publicly_available = FALSE),
    scientific_fit_count = 0L, generated_at = format(Sys.time(), tz = "UTC",
      usetz = TRUE))
  write_json(manifest, file.path(output, "manifest.json"))
  message("Wrote blocked-run manifest; no scientific fit was launched.")
  invisible(manifest)
}

if (stage == "prepare") { prepare(); quit(save = "no") }
if (stage == "smoke") {
  if (nzchar(fit_id)) run_one(fit_id, TRUE) else
    invisible(lapply(registry$fit_id, run_one, smoke = TRUE))
  quit(save = "no")
}
if (stage == "fit") {
  if (!nzchar(fit_id)) stop("--fit is required for --stage fit.")
  run_one(fit_id, FALSE); quit(save = "no")
}
if (stage == "manifest") {
  write_blocked_manifest(); quit(save = "no")
}
if (stage == "all") {
  if (!file.exists(bundle_path)) prepare()
  invisible(lapply(registry$fit_id, run_one, smoke = TRUE))
  invisible(lapply(registry$fit_id, run_one, smoke = FALSE))
  quit(save = "no")
}
stop("Unknown stage: ", stage)
