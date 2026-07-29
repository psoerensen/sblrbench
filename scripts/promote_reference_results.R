local_dir <- file.path("results", "local", "01_finemapping", "separated")
snapshot_dir <- file.path("results", "reference", "01_finemapping",
                          "separated-development-v1")
expected <- c("marker_metrics.csv", "credible_set_metrics.csv",
              "credible_set_summary.csv", "computational_summary.csv",
              "replicate_status.csv", "pilot_manifest.json")

fail <- function(...) stop(..., call. = FALSE)
if (!all(file.exists(file.path(local_dir, expected))))
  fail("Missing local result files: ", paste(expected[!file.exists(file.path(local_dir, expected))], collapse = ", "))
if (dir.exists(snapshot_dir) && !identical(tolower(Sys.getenv("SBLRBENCH_OVERWRITE_REFERENCE")), "true"))
  fail("Reference snapshot already exists. Set SBLRBENCH_OVERWRITE_REFERENCE=true to replace it.")

manifest <- jsonlite::read_json(file.path(local_dir, "pilot_manifest.json"), simplifyVector = TRUE)
marker <- read.csv(file.path(local_dir, "marker_metrics.csv"), check.names = FALSE)
cs_metrics <- read.csv(file.path(local_dir, "credible_set_metrics.csv"), check.names = FALSE)
cs_summary <- read.csv(file.path(local_dir, "credible_set_summary.csv"), check.names = FALSE)
computational <- read.csv(file.path(local_dir, "computational_summary.csv"), check.names = FALSE)
status <- read.csv(file.path(local_dir, "replicate_status.csv"), check.names = FALSE)

required_metrics <- c("pip_brier", "effect_rmse", "average_precision",
  "causal_rank_mean", "causal_rank_median", "causal_rank_best",
  "causal_rank_worst", "causal_top_10_recall", "causal_top_20_recall",
  "causal_top_50_recall")
if (!identical(manifest$study, "01_finemapping") || !identical(manifest$architecture, "separated") ||
    !isTRUE(manifest$development_settings) || manifest$replicate_count != 10L || length(manifest$methods) != 4L)
  fail("Manifest does not identify the reviewed 10-replicate separated development pilot.")
if (nrow(status) != 40L || any(status$status != "ok") || nrow(computational) != 40L ||
    any(computational$status != "ok")) fail("The 40 required method fits are not all successful.")
if (nrow(marker) != 400L || length(setdiff(required_metrics, unique(marker$metric))) ||
    any(!is.finite(marker$value))) fail("Marker metrics are incomplete or non-finite.")
if (nrow(unique(cs_summary[c("replicate", "method")])) != 40L || !nrow(cs_metrics))
  fail("Credible-set summaries do not cover all method-replicate combinations.")

scan_paths <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  grepl("[A-Za-z]:[/\\\\]Users[/\\\\]|/Users/|/home/", txt, perl = TRUE)
}
if (any(vapply(file.path(local_dir, expected), scan_paths, logical(1))))
  fail("An absolute user path was found in compact local results.")

if (dir.exists(snapshot_dir)) unlink(snapshot_dir, recursive = TRUE)
dir.create(snapshot_dir, recursive = TRUE)
for (f in setdiff(expected, "pilot_manifest.json"))
  file.copy(file.path(local_dir, f), file.path(snapshot_dir, f), overwrite = FALSE)

repo_commit <- trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
fit_commit <- unique(computational$sblrbench_commit)
if (length(fit_commit) != 1L) fail("Fit records contain inconsistent sblrbench commits.")
sblr_commit <- tryCatch(trimws(system2("git", c("-C", "../sblr", "rev-parse", "HEAD"), stdout = TRUE)), error = function(e) NA_character_)
manifest$benchmark_version <- "separated-development-v1"
manifest$source_snapshot_commit <- repo_commit
manifest$fit_provenance_commit <- fit_commit
manifest$sblr_git_commit <- sblr_commit
manifest$canonical_marker_count <- 37991L
manifest$validation <- list(successful_fits = 40L, expected_fits = 40L,
  marker_metric_rows = 400L, credible_set_method_replicates = 40L,
  exact_causals_verified = TRUE, oracle_checks_passed = 10L,
  marker_alignment_validated = TRUE, sample_alignment_validated = TRUE,
  minimum_observed_causal_distance_bp = 1118217L,
  distinct_causal_sets = 10L)
jsonlite::write_json(manifest, file.path(snapshot_dir, "benchmark_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)

method_labels <- c(st_bed_bayesc = "ST-BED BayesC", st_bed_bayesr = "ST-BED BayesR",
  st_csr_sbayesc = "ST-CSR SBayesC", st_csr_sbayesr = "ST-CSR SBayesR")
metric_stat <- function(method, metric, fun) fun(marker$value[marker$method == method & marker$metric == metric])
cs_stat <- function(method, metric, fun) fun(cs_summary$value[cs_summary$method == method & cs_summary$metric == metric])
summary_rows <- lapply(manifest$methods, function(method) {
  comp <- computational[computational$method == method, ]
  data.frame(method = method, method_label = unname(method_labels[method]),
    n_replicates = length(unique(comp$replicate)), n_successful_fits = sum(comp$status == "ok"),
    n_failed_fits = sum(comp$status != "ok"),
    average_precision_mean = metric_stat(method, "average_precision", mean), average_precision_sd = metric_stat(method, "average_precision", sd),
    pip_brier_mean = metric_stat(method, "pip_brier", mean), pip_brier_sd = metric_stat(method, "pip_brier", sd),
    effect_rmse_mean = metric_stat(method, "effect_rmse", mean), effect_rmse_sd = metric_stat(method, "effect_rmse", sd),
    causal_top_10_recall_mean = metric_stat(method, "causal_top_10_recall", mean), causal_top_10_recall_sd = metric_stat(method, "causal_top_10_recall", sd),
    causal_top_20_recall_mean = metric_stat(method, "causal_top_20_recall", mean), causal_top_20_recall_sd = metric_stat(method, "causal_top_20_recall", sd),
    causal_top_50_recall_mean = metric_stat(method, "causal_top_50_recall", mean), causal_top_50_recall_sd = metric_stat(method, "causal_top_50_recall", sd),
    causal_rank_median_mean = metric_stat(method, "causal_rank_median", mean), causal_rank_worst_mean = metric_stat(method, "causal_rank_worst", mean),
    credible_sets_mean = cs_stat(method, "number_of_credible_sets", mean), credible_sets_sd = cs_stat(method, "number_of_credible_sets", sd),
    causal_loci_detected_mean = cs_stat(method, "causal_loci_detected", mean), causal_loci_detected_sd = cs_stat(method, "causal_loci_detected", sd),
    causal_locus_detection_fraction = cs_stat(method, "causal_loci_detected", mean) / manifest$simulation$n_causal,
    conditional_exact_coverage_mean = cs_stat(method, "exact_coverage_fraction", mean),
    conditional_ld_proxy_coverage_mean = cs_stat(method, "ld_proxy_coverage_fraction", mean),
    mean_set_size = cs_stat(method, "mean_set_size", mean), runtime_mean_seconds = mean(comp$runtime),
    runtime_sd_seconds = sd(comp$runtime), runtime_min_seconds = min(comp$runtime), runtime_max_seconds = max(comp$runtime),
    stringsAsFactors = FALSE)
})
write.csv(do.call(rbind, summary_rows), file.path(snapshot_dir, "benchmark_summary.csv"), row.names = FALSE)

config <- manifest[c("study", "architecture", "development_settings")]
config$replicate_count <- as.integer(manifest$replicate_count); config$methods <- manifest$methods
config$simulation <- manifest$simulation; config$causal_selection <- manifest$causal_selection
config$mcmc <- manifest$mcmc; config$credible_sets <- manifest$credible_sets
config_lines <- c("# This is a frozen configuration supporting the published reference benchmark.",
  "# Do not edit this file in place. Create a new benchmark snapshot for changed settings.",
  capture.output(dput(config)))
writeLines(config_lines, file.path(snapshot_dir, "config.R"))

run_lines <- c("# Exact workflow entry point for the separated development benchmark.",
  "# Genotype input is not redistributed. This reproduces the workflow only when", "# the expected input is available; a compatible Glist may be substituted for testing.",
  "# The frozen reference files reproduce the report without genotype data.", "Sys.setenv(",
  "  SBLR_BENCH_STUDY = \"01_finemapping\",", "  SBLR_BENCH_REPLICATES = \"10\"", ")", "", "targets::tar_make()", "",
  "targets::tar_manifest()", "targets::tar_read(marker_metrics)",
  "targets::tar_read(credible_set_summary)", "targets::tar_read(computational_summary)", "targets::tar_read(replicate_status)")
writeLines(run_lines, file.path(snapshot_dir, "run_benchmark.R"))

minimal <- c("# Minimal, fast contract demonstration; no large Glist or sampler fit is required.",
  "library(sblr)", "library(sblrbench)", "", "set.seed(42)",
  "W <- matrix(sample(0:2, 240, replace = TRUE), nrow = 30, ncol = 8,",
  "  dimnames = list(paste0(\"sample\", 1:30), paste0(\"marker\", 1:8)))",
  "sim_raw <- sblr::mtsim(W = W, nt = 1L, n_shared = 2L, n_specific = 0L,",
  "  h2 = 0.2, seed = 2001L)",
  "simulation <- as_sblrbench_simulation(sim_raw, study = \"minimal_example\",",
  "  architecture = \"synthetic\", replicate = 1L, seed = 2001L)",
  "check_oracle_genetic_values(simulation)", "", "# Lightweight standard-result demonstration using known truth.",
  "causal <- matrix(as.numeric(rownames(simulation$truth$effects) %in% sim_raw$causal$all),",
  "  ncol = 1L, dimnames = dimnames(simulation$truth$effects))",
  "result <- new_sblrbench_result(\"truth_demo\", effects = simulation$truth$effects, pip = causal)",
  "metric_effect_rmse(simulation, result)", "metric_pip_brier(simulation, result)",
  "# A full public sblr sampler requires its corresponding BED/Glist or CSR/LD input;")
writeLines(minimal, file.path(snapshot_dir, "minimal_example.R"))

source_paths <- c("_targets.R", "studies/01_finemapping/config.R", "studies/01_finemapping/targets.R",
  "studies/01_finemapping/pilot.R", "R/metrics.R", "R/alignment.R", "R/simulation.R", "R/adapter-sblr.R")
roles <- c("dispatcher", "study configuration", "targets pipeline", "study helpers", "metric implementation",
  "alignment contract", "simulation contract", "native adapter")
blob <- vapply(source_paths, function(p) trimws(system2("git", c("rev-parse", paste0("HEAD:", p)), stdout = TRUE)), character(1))
source_inventory <- data.frame(path = source_paths, git_blob_sha = blob,
  file_md5 = unname(tools::md5sum(source_paths)), role = roles, stringsAsFactors = FALSE)
write.csv(source_inventory, file.path(snapshot_dir, "source_files.csv"), row.names = FALSE)

quarto <- tryCatch(trimws(system2("quarto", "--version", stdout = TRUE, stderr = FALSE)), error = function(e) "unavailable")
session <- c(capture.output(sessionInfo()), "", paste("sblr:", packageVersion("sblr")),
  paste("sblrbench:", packageVersion("sblrbench")), paste("qgg:", packageVersion("qgg")),
  paste("targets:", packageVersion("targets")), paste("Quarto:", quarto),
  paste("sblrbench source snapshot commit:", repo_commit), paste("fit provenance commit:", fit_commit),
  paste("sblr Git commit:", sblr_commit))
writeLines(session, file.path(snapshot_dir, "session_info.txt"))

readme <- c("# Separated-locus fine-mapping development pilot", "",
  "This frozen reference capsule records a 10-replicate structural benchmark of four single-trait sblr methods on 37,991 canonical markers with 10 separated causal markers per replicate.", "",
  "The 500-iteration, 250-burn-in, one-chain settings are for development validation only. The example markers have limited LD; these results do not support scientific method rankings.", "",
  "## Workflow reproduction", "", "The code and configuration can be used to rerun the same analytical workflow with compatible genotype input.", "",
  "## Exact numerical reproduction", "", "Exact numerical reproduction requires the same original genotype data, software versions, seeds and computing environment. The genotype data are not included in the public reference snapshot.", "",
  "Run `source(\"minimal_example.R\")` for the small synthetic contract demonstration. Run `source(\"run_benchmark.R\")` only in a suitable clone with the required non-redistributed genotype/Glist input.", "",
  "The capsule includes compact metrics, summaries, statuses, the frozen manifest/configuration, reproduction scripts, source inventory, session information and checksums. It excludes genotype data, sparse LD, fits, posterior samples and `_targets/`.", "",
  paste0("Source snapshot: `", repo_commit, "`; fit provenance: `", fit_commit, "`; sblr `", manifest$package_versions$sblr, "` (`", sblr_commit, "`); qgg `", manifest$package_versions$qgg, "`."), "",
  "Verify a file with `tools::md5sum()` and compare it with `checksums.csv`. The checksum table intentionally does not checksum itself.")
writeLines(readme, file.path(snapshot_dir, "README.md"))

checksum_files <- setdiff(list.files(snapshot_dir, recursive = TRUE), "checksums.csv")
info <- file.info(file.path(snapshot_dir, checksum_files))
checksums <- data.frame(file = checksum_files, size_bytes = info$size,
  md5 = unname(tools::md5sum(file.path(snapshot_dir, checksum_files))), row.names = NULL)
write.csv(checksums, file.path(snapshot_dir, "checksums.csv"), row.names = FALSE)
message("Created validated reference snapshot: ", snapshot_dir)
