#!/usr/bin/env Rscript

# External-reference diagnostic for official zhilizheng/SBayesRC. This script
# cannot update the formal Study 06 qualification or launch final mode.

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION")))
  stop("Run this script from the sblrbench repository root.", call. = FALSE)

value_after <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) return(default)
  args[index + 1L]
}
validate_only <- "--validate-only" %in% args
export_only <- "--export-only" %in% args
smoke <- "--smoke" %in% args
child <- "--child" %in% args
condition_arg <- value_after("--condition")
chain_arg <- as.integer(value_after("--chain", NA_character_))

output <- file.path(root, "results", "local", "06_annotation_models",
  "gctb_parity")
official_library <- file.path(output, "rlib")
benchmark_library <- file.path(root, "results", "local",
  "current_benchmark_refresh", "rlib")
if (!dir.exists(official_library) || !dir.exists(benchmark_library))
  stop("The official or benchmark isolated R library is unavailable.",
    call. = FALSE)
.libPaths(c(official_library, benchmark_library, .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1")

pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "studies/06_annotation_models/power-isolation.R"),
  local = FALSE)
source(file.path(root, "studies/06_annotation_models/gctb-parity.R"),
  local = FALSE)
profile <- study06_gctb_profile()
spec <- read_benchmark_spec(file.path(root,
  "studies/06_annotation_models/spec.R"))
if (!identical(benchmark_annotation_spec_hash(spec), profile$spec_hash))
  stop("The Study 06 specification hash differs from the parity contract.",
    call. = FALSE)
if (!identical(as.character(packageVersion("SBayesRC")), "0.2.6") ||
    !startsWith(normalizePath(find.package("SBayesRC"), winslash = "/"),
      normalizePath(official_library, winslash = "/")))
  stop("Official SBayesRC 0.2.6 is not loaded from the isolated library.",
    call. = FALSE)
if (!identical(spec$packages$sblr$sha,
    benchmark_package_provenance("sblr")$sha))
  stop("The loaded sblr package does not match the Study 06 pin.",
    call. = FALSE)

glist_cache <- file.path(root, "results/local/06_annotation_models",
  "checkpoints/data/human_glist.rds")
if (!nzchar(Sys.getenv("SBLR_BENCH_GLIST", "")) && file.exists(glist_cache))
  Sys.setenv(SBLR_BENCH_GLIST = glist_cache)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

if (child) {
  registry_path <- file.path(output, "registry.csv")
  manifest_path <- file.path(output, "export", "manifest.json")
  if (!file.exists(registry_path) || !file.exists(manifest_path))
    stop("Child mode requires the validated parent export and registry.",
      call. = FALSE)
  registry <- data.table::fread(registry_path, data.table = FALSE)
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  truth_hash <- manifest$truth_hash
  if (!identical(truth_hash, profile$truth_hash) ||
      !identical(manifest$specification_hash, profile$spec_hash) ||
      !identical(manifest$official_package_sha, profile$official_sha))
    stop("Child export identity differs from the parity contract.",
      call. = FALSE)
} else {
  data <- prepare_prediction_data(spec, output)
  logic <- .annotation_logic(spec)
  annotations <- logic$construct_annotation_design(data$markers$marker_ids,
    spec)
  annotation_truth <- logic$construct_annotation_truth(annotations, spec)
  seed_grid <- benchmark_annotation_seeds(spec, "benchmark", "qualification")
  seed_row <- seed_grid[seed_grid$scenario == "informative_annotations" &
    seed_grid$replicate == 1L & seed_grid$method == "st_bed_bayesrc", ,
    drop = FALSE]
  coordinate <- as.list(seed_row[1L, c("scenario", "replicate",
    "component_seed", "effect_seed", "residual_seed")])
  simulation <- logic$simulate_annotation_architecture(coordinate,
    data$scaled$all, data$split$train_rows, annotations, annotation_truth, spec)
  stats <- benchmark_summary_stats(simulation, data$ld_glist, data$split,
    spec$data)
  bundle <- list(spec = spec, simulation = simulation,
    stats = stats, annotations = annotations,
    annotation_truth = annotation_truth,
    marker_truth = simulation$extras$marker_truth)
  truth_identity <- study06_power_truth_identity(data, bundle)
  truth_hash <- benchmark_hash_object(truth_identity)
  if (!identical(truth_hash, profile$truth_hash))
    stop("The reconstructed informative truth hash differs from the contract: ",
      truth_hash, call. = FALSE)
  chain_seeds <- as.integer(seed_row$chain_seeds[[1L]])
  registry <- study06_gctb_registry(chain_seeds, profile)
  data.table::fwrite(registry, file.path(output, "registry.csv"))
  exported <- study06_gctb_export(data, simulation, stats, annotations, spec,
    output, official_library)
}

if (validate_only || export_only) {
  message("Validated Study 06 official export and 12-coordinate registry; ",
    "no official fit run.")
  quit(save = "no", status = 0L)
}

run_child <- function(condition, chain, is_smoke) {
  row <- registry[registry$condition == condition & registry$chain == chain, ]
  if (nrow(row) != 1L) stop("Unknown official diagnostic coordinate.",
    call. = FALSE)
  run_dir <- file.path(output, if (is_smoke) "smoke" else "runs",
    row$fit_id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(run_dir, row$fit_id)
  for (path in list.files(run_dir, full.names = TRUE, all.files = TRUE,
      no.. = TRUE)) unlink(path, recursive = TRUE, force = TRUE)
  set.seed(row$requested_seed)
  use_annotation <- condition != "G0"
  niter <- if (is_smoke) profile$smoke_niter else profile$niter
  burn <- if (is_smoke) profile$smoke_burn else profile$burn
  started <- Sys.time()
  warning_messages <- character()
  error_message <- NULL
  tryCatch(withCallingHandlers(
    SBayesRC::sbayesrc(
      mafile = file.path(output, "export", "study06_informative.ma"),
      LDdir = file.path(output, "export", "ld"), outPrefix = prefix,
      annot = if (use_annotation) file.path(output, "export",
        "study06_informative.annot") else "",
      log2file = TRUE, bTune = condition == "G2",
      tuneIter = profile$tune_iter, tuneBurn = profile$tune_burn,
      thresh = if (condition == "G2") .995 else 1,
      tuneStep = profile$tune_step, bTunePrior = FALSE,
      niter = niter, burn = burn, starth2 = .5,
      startPi = if (condition == "G2") profile$start_pi_native else
        profile$start_pi_matched,
      gamma = if (condition == "G2") profile$gamma_native else
        profile$gamma_matched,
      sSamVe = "allMixVe", twopq = "nbsq", bOutDetail = TRUE,
      seed = row$requested_seed, outFreq = 1L, annoSigmaScale = 1,
      bOutBeta = TRUE), warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) error_message <<- conditionMessage(e))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  rds_path <- paste0(prefix, ".rds")
  fit <- if (is.null(error_message) && file.exists(rds_path))
    readRDS(rds_path) else NULL
  signature <- if (is.null(fit)) NA_character_ else benchmark_hash_object(list(
    fit$betaMean, fit$pip, fit$n_comp_hist, fit$pi_hist, fit$hsq_hist,
    fit$vg_hist, fit$ve_hist, fit$alpha))
  tune_path <- file.path(paste0(prefix, ".mcmcsamples"),
    paste0(row$fit_id, "_tune.txt"))
  tuning <- if (file.exists(tune_path)) data.table::fread(tune_path) else NULL
  selected_threshold <- if (is.null(tuning)) if (condition == "G2")
    NA_real_ else 1 else {
      candidate <- tuning[is.finite(r) & is.finite(rel_r)]
      if (!nrow(candidate)) NA_real_ else {
        best <- candidate[which.max(rel_r)]
        if (best$rel_r > 1.25 || candidate[thresh == max(thresh)]$r < 0)
          best$thresh else max(profile$tune_step)
      }
    }
  record <- list(schema = "sblrbench-study06-gctb-smoke-v1",
    condition = condition, chain = chain,
    seed_source = row$seed_source, requested_seed = row$requested_seed,
    warning = warning_messages, error = error_message,
    status = if (is.null(error_message) && !is.null(fit)) "complete" else
      "failed", wall_seconds = elapsed, result_signature = signature,
    selected_threshold = selected_threshold,
    returned_names = if (is.null(fit)) character() else names(fit),
    official_version = as.character(packageVersion("SBayesRC")),
    official_sha = profile$official_sha, truth_hash = truth_hash,
    specification_hash = profile$spec_hash)
  saveRDS(record, file.path(run_dir, "smoke-record.rds"))
  jsonlite::write_json(record, file.path(run_dir, "smoke-record.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!is.null(error_message)) stop(error_message, call. = FALSE)
}

if (child) {
  if (!condition_arg %in% c("G0", "G1", "G2") || is.na(chain_arg))
    stop("Child mode requires --condition G0/G1/G2 and --chain 1..4.",
      call. = FALSE)
  run_child(condition_arg, chain_arg, smoke)
  quit(save = "no", status = 0L)
}

if (!smoke)
  stop("Full official registry is not started implicitly; run --smoke first.",
    call. = FALSE)

# Two G0 processes with distinct preregistered seeds audit whether the official
# native RNG honors the wrapper's public seed contract. One G1 and one G2 smoke
# then validate the annotation and native-workflow paths.
coordinates <- data.frame(condition = c("G0", "G0", "G1", "G2"),
  chain = c(1L, 2L, 1L, 1L))
rscript <- file.path(R.home("bin"), "Rscript.exe")
statuses <- integer(nrow(coordinates))
for (i in seq_len(nrow(coordinates))) {
  log <- file.path(output, "smoke", paste0(coordinates$condition[i],
    "--chain", coordinates$chain[i], ".process.log"))
  dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  statuses[i] <- system2(rscript, c(file.path(root,
    "scripts/run_study06_gctb_parity.R"), "--child", "--smoke",
    "--condition", coordinates$condition[i], "--chain",
    coordinates$chain[i]), stdout = log, stderr = log)
}
records <- lapply(seq_len(nrow(coordinates)), function(i) readRDS(file.path(
  output, "smoke", paste0(coordinates$condition[i], "--chain",
    coordinates$chain[i]), "smoke-record.rds")))
g0_identical <- identical(records[[1L]]$result_signature,
  records[[2L]]$result_signature)
decision <- list(schema = "sblrbench-study06-gctb-parity-decision-v1",
  decision = if (all(statuses == 0L) && g0_identical) "GCTB-P5" else
    if (any(statuses != 0L)) "GCTB-P5" else "smoke-passed",
  classification = if (g0_identical)
    "blocked: official public seed is not connected to native RNG" else
    if (any(statuses != 0L)) "blocked: official smoke failure" else
      "smoke passed; full registry separately executable",
  formal_qualification_unchanged = TRUE, final_benchmark_authorized = FALSE,
  specification_hash = profile$spec_hash, truth_hash = profile$truth_hash,
  official_sha = profile$official_sha, process_status = statuses,
  requested_seeds = vapply(records, `[[`, integer(1), "requested_seed"),
  result_signatures = vapply(records, `[[`, character(1),
    "result_signature"), g0_distinct_requested_seeds =
      records[[1L]]$requested_seed != records[[2L]]$requested_seed,
  g0_fresh_process_results_identical = g0_identical,
  full_registry_run = FALSE,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
jsonlite::write_json(decision, file.path(output, "diagnostic_decision.json"),
  pretty = TRUE, auto_unbox = TRUE)
if (g0_identical)
  message("Official smoke completed, but independent chains are blocked: ",
    "distinct wrapper seeds produced identical native results.") else
  message("Official smoke completed without detecting an RNG identity blocker.")
