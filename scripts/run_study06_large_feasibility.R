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
if (!phase %in% c("initial", "continuation", "gctb_block"))
  stop("--phase must be initial, continuation, or gctb_block.")
default_output <- if (identical(phase, "gctb_block"))
  "results/local/06_annotation_models/large_feasibility/gctb_block" else
    "results/local/06_annotation_models/large_feasibility"
output <- value("--output-dir", default_output)
isolated <- normalizePath(value("--library",
  "results/local/06_annotation_models/gctb_block_contract_validation/rlib"),
  winslash = "/", mustWork = TRUE)
.libPaths(c(isolated, .libPaths()))

if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required.")
pkgload::load_all(".", quiet = TRUE)
source("studies/06_annotation_models/large-feasibility.R")
spec <- study06_large_spec()
registry <- study06_large_registry(spec)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

required_sblr_sha <- "0c89234273389e14112ba0e08ef9d11d3e1819dc"
required_sblr_tree <- "e723528e7d5d570a31b5b1d1c90551896ac48f86ab05261c181c8109af971fd0"
source_sha <- trimws(system2("git", c("-C", "../sblr", "rev-parse", "HEAD"),
  stdout = TRUE))
source_status <- system2("git", c("-C", "../sblr", "status", "--short"),
  stdout = TRUE)
installed_provenance <- benchmark_package_provenance("sblr", lib.loc = isolated)
if (!identical(source_sha, required_sblr_sha) || length(source_status) ||
    !identical(as.character(utils::packageVersion("sblr", lib.loc = isolated)),
      "0.2.0") ||
    !identical(installed_provenance$installed_tree_sha256,
      required_sblr_tree)) {
  stop("Exact clean sblr source or isolated installed-tree identity mismatch.")
}

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
default_bundle <- if (identical(phase, "gctb_block"))
  "results/local/06_annotation_models/large_feasibility/prepared_bundle.rds" else
    file.path(output, "prepared_bundle.rds")
bundle_path <- normalizePath(value("--bundle", default_bundle),
  winslash = "/", mustWork = stage != "prepare")

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
  expected <- c(
    specification_hash = "b001bc36a5531e5e6b342286a253fc1fd34dad4265359d89d2feaa026d4533df",
    truth_hash = "e94a511540f600e61ef47b52947836f19a15388f5e8ce795c179929956817507",
    marker_order_hash = "e394f4324f89ca7ad88691284a6216a24e26c691d7cb59d822d5d7f006b096f2",
    sample_order_hash = "f3c14e98845fc565424a973a5a5850768c587f9d5a5306cc58ff36c5329a256b",
    allele_hash = "996fb56147c14801f798a7bb21692a962aba18829c8a3a7abfa4dc768b06b082",
    annotation_hash = "5cc2a5dc64140b5cb2c3e77044d8a0b7bd698e2b33f1f8cb0f5b39ac278f7bd3",
    block_hash = "3b119c38fcababc5b70892f8fc29dcd773fb73e4e2a2b1cd6daf94556e16294d",
    gwas_hash = "f00c01326d20d32c1b389f239ebe520327f4b8c859acf7320b10d531243cacfd",
    alpha_truth_hash = "4766d00b77653825e9130a32ebcde1b16754ee99103f2ab4c4d3f1d715fbbf82",
    trace_panel_hash = "0ae8cc37d0418d54cf52e4cf5271c5859d01506759a6d798f6c512a66b08438f")
  got <- unlist(x$identities[names(expected)])
  if (!identical(unname(got), unname(expected)) ||
      length(x$sample_ids) != 5000L || length(x$markers$marker_ids) != 37991L ||
      length(x$blocks$block_start) != 76L ||
      !identical(as.integer(x$truth$component_count),
        c(36791L, 618L, 392L, 190L)) ||
      sum(x$truth$marker_truth$active) != 1200L ||
      !isTRUE(all.equal(x$truth$realized_h2, .5, tolerance = 0)))
    stop("Frozen large bundle identity or truth-count audit failed.")
  x
}

run_one <- function(id, smoke = FALSE) {
  bundle <- load_bundle()
  row <- registry[registry$fit_id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Unknown fit id: ", id)
  controls <- study06_large_controls(row, bundle$calibrated$alpha, spec, smoke)
  if (smoke && identical(phase, "gctb_block")) {
    controls$seed <- spec$seeds$fit
    controls$chain_seeds <- spec$seeds$chain
  }
  suffix <- if (smoke) "smoke" else "fit"
  checkpoint_name <- if (phase == "initial") paste0(id, "_", suffix, ".rds") else
    paste0(phase, "_", id, "_", suffix, ".rds")
  path <- file.path(output, "checkpoints", checkpoint_name)
  identity <- list(schema = spec$schema, fit = as.list(row),
    specification_hash = bundle$identities$specification_hash,
    truth_hash = bundle$identities$truth_hash,
    marker_hash = bundle$identities$marker_order_hash,
    annotation_hash = bundle$identities$annotation_hash,
    block_hash = bundle$identities$block_hash, gwas_hash = bundle$identities$gwas_hash,
    trace_hash = bundle$identities$trace_panel_hash, controls = controls,
    block_contract = if (row$route == "block_eigen")
      study06_large_block_contract() else NULL,
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
  if (smoke) study06_large_validate_smoke(fit, row,
    bundle$calibrated$alpha, nrow(bundle$trace_panel$panel), 12L)
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
if (stage == "analyze") {
  source("studies/06_annotation_models/large-feasibility-analysis.R")
  bundle <- load_bundle()
  checkpoint_dir <- file.path(output, "checkpoints")
  analysis_dir <- file.path(output, "analysis")
  result <- study06_large_analyze_completed(checkpoint_dir, analysis_dir,
    bundle, spec)
  message("Analyzed six completed identity-matched fits; no MCMC was run.")
  quit(save = "no")
}
if (stage == "all") {
  if (!file.exists(bundle_path)) prepare()
  invisible(lapply(registry$fit_id, run_one, smoke = TRUE))
  invisible(lapply(registry$fit_id, run_one, smoke = FALSE))
  quit(save = "no")
}
stop("Unknown stage: ", stage)
