# Promote the reviewed Study 01 pilot into the improved v1.1 capsule.
# Numerical outputs are copied from immutable v1. Benchmark-specific data
# metadata are optional for future studies; this benchmark uses public qgdata.

source_dir <- file.path("results", "reference", "01_finemapping",
                        "separated-development-v1")
snapshot_dir <- file.path("results", "reference", "01_finemapping",
                          "separated-development-v1.1")
overwrite <- identical(tolower(Sys.getenv("SBLRBENCH_OVERWRITE_REFERENCE")), "true")
fail <- function(...) stop(..., call. = FALSE)

if (!dir.exists(source_dir)) fail("Immutable v1 reference capsule is missing.")
if (dir.exists(snapshot_dir) && !overwrite)
  fail("Reference snapshot already exists. Set SBLRBENCH_OVERWRITE_REFERENCE=true to replace it.")

v1_checksums <- read.csv(file.path(source_dir, "checksums.csv"), stringsAsFactors = FALSE)
actual_md5 <- unname(tools::md5sum(file.path(source_dir, v1_checksums$file)))
if (!identical(actual_md5, v1_checksums$md5)) fail("Immutable v1 checksum validation failed.")

result_files <- c("benchmark_summary.csv", "marker_metrics.csv",
  "credible_set_metrics.csv", "credible_set_summary.csv",
  "computational_summary.csv", "replicate_status.csv")
if (!all(file.exists(file.path(source_dir, result_files)))) fail("v1 result files are incomplete.")

marker <- read.csv(file.path(source_dir, "marker_metrics.csv"), check.names = FALSE)
cs_summary <- read.csv(file.path(source_dir, "credible_set_summary.csv"), check.names = FALSE)
computational <- read.csv(file.path(source_dir, "computational_summary.csv"), check.names = FALSE)
status <- read.csv(file.path(source_dir, "replicate_status.csv"), check.names = FALSE)
required_metrics <- c("pip_brier", "effect_rmse", "average_precision",
  "causal_rank_mean", "causal_rank_median", "causal_rank_best",
  "causal_rank_worst", "causal_top_10_recall", "causal_top_20_recall",
  "causal_top_50_recall")
if (nrow(status) != 40L || any(status$status != "ok") || nrow(computational) != 40L ||
    any(computational$status != "ok") || nrow(marker) != 400L ||
    any(!is.finite(marker$value)) || length(setdiff(required_metrics, marker$metric)) ||
    nrow(unique(cs_summary[c("replicate", "method")])) != 40L)
  fail("Reviewed benchmark results failed promotion validation.")

data_source <- list(
  type = "public",
  downloadable = TRUE,
  repository = "psoerensen/qgdata",
  commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
  subdirectory = "simulated_human_data",
  licence = NA_character_
)
data_files <- data.frame(
  filename = c("human.bed", "human.bim", "human.fam", "human.pheno", "human.covar"),
  size_bytes = c(62500003, 1882359, 117786, 92786, 641513),
  md5 = c("e89bea9a6cedd9eeef3fd0a5c807db81", "0105119b04c67b7ac7f66cc5e6680963",
    "3c5db3d9eb7f3fc893c75f6f2b89836d", "6a9e7cb1162e43999c170a363863176d",
    "d06002aa2b1b79bdc4c0e92c21f27ae5"),
  role = c("PLINK genotype", "PLINK marker metadata", "PLINK sample metadata",
    "example phenotype", "example covariates"), stringsAsFactors = FALSE
)
data_files$repository <- data_source$repository
data_files$commit <- data_source$commit
data_files$path <- file.path(data_source$subdirectory, data_files$filename)
data_files$download_url <- paste0("https://raw.githubusercontent.com/",
  data_source$repository, "/", data_source$commit, "/", data_files$path)
data_manifest <- data_files[c("repository", "commit", "path", "filename",
  "size_bytes", "md5", "download_url", "role")]

if (dir.exists(snapshot_dir)) unlink(snapshot_dir, recursive = TRUE)
dir.create(snapshot_dir, recursive = TRUE)
for (file in result_files)
  if (!file.copy(file.path(source_dir, file), file.path(snapshot_dir, file)))
    fail("Could not copy reviewed result file: ", file)
write.csv(data_manifest, file.path(snapshot_dir, "example_data_manifest.csv"), row.names = FALSE)

manifest <- jsonlite::read_json(file.path(source_dir, "benchmark_manifest.json"), simplifyVector = TRUE)
manifest$benchmark_version <- "separated-development-v1.1"
manifest$results_identical_to <- "separated-development-v1"
manifest$reproducibility_update <- "Same numerical benchmark results; improved public-data provenance and worked examples."
manifest$public_data_source <- data_source
jsonlite::write_json(manifest, file.path(snapshot_dir, "benchmark_manifest.json"),
  pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)

config <- source(file.path(source_dir, "config.R"), local = TRUE)$value
config$example_data <- list(repository = data_source$repository,
  commit = data_source$commit, subdirectory = data_source$subdirectory,
  files = data_files$filename, size_bytes = stats::setNames(data_files$size_bytes, data_files$filename),
  md5 = stats::setNames(data_files$md5, data_files$filename))
writeLines(c("# Frozen configuration for reference benchmark v1.1.",
  "# Do not edit in place; create a new snapshot for changed settings.",
  capture.output(dput(config))), file.path(snapshot_dir, "config.R"))

invisible(file.copy("scripts/worked_finemapping_example.R",
          file.path(snapshot_dir, "worked_finemapping_example.R")))
invisible(file.copy("scripts/contract_smoke_test.R",
          file.path(snapshot_dir, "contract_smoke_test.R")))
invisible(file.copy("scripts/run_benchmark.R", file.path(snapshot_dir, "run_benchmark.R")))

source_paths <- c("_targets.R", "studies/01_finemapping/config.R",
  "studies/01_finemapping/targets.R", "studies/01_finemapping/setup_example_data.R",
  "studies/01_finemapping/pilot.R", "R/example-data.R", "R/metrics.R",
  "R/alignment.R", "R/simulation.R", "R/adapter-sblr.R",
  "scripts/worked_finemapping_example.R", "scripts/contract_smoke_test.R",
  "scripts/run_benchmark.R")
roles <- c("dispatcher", "study configuration", "targets pipeline", "public data setup",
  "study helpers", "pinned public-data downloader", "metric implementation",
  "alignment contract", "simulation contract", "native adapter", "worked user example",
  "developer contract smoke test", "complete benchmark entry point")
git_blob <- vapply(source_paths, function(path)
  trimws(system2("git", c("hash-object", path), stdout = TRUE)), character(1))
source_inventory <- data.frame(path = source_paths, git_blob_sha = git_blob,
  file_md5 = unname(tools::md5sum(source_paths)), role = roles, stringsAsFactors = FALSE)
write.csv(source_inventory, file.path(snapshot_dir, "source_files.csv"), row.names = FALSE)

quarto <- tryCatch(trimws(system2("quarto", "--version", stdout = TRUE, stderr = FALSE)),
                   error = function(e) "unavailable")
session_lines <- c(capture.output(sessionInfo()), "",
  paste("sblr:", packageVersion("sblr")), paste("sblrbench:", packageVersion("sblrbench")),
  paste("qgg:", packageVersion("qgg")), paste("targets:", packageVersion("targets")),
  paste("Quarto:", quarto), paste("benchmark source snapshot commit:", manifest$source_snapshot_commit),
  paste("fit provenance commit:", manifest$fit_provenance_commit),
  paste("sblr Git commit:", manifest$sblr_git_commit),
  paste("qgdata commit:", data_source$commit))
writeLines(session_lines, file.path(snapshot_dir, "session_info.txt"))

readme <- c("# Separated-locus fine-mapping development pilot v1.1", "",
  "This capsule contains the same numerical benchmark results as `separated-development-v1`. Version 1.1 improves public-data provenance, examples, and reproduction guidance; it does not rerun or alter the benchmark.", "",
  "The benchmark remains a structural development run with limited LD, 500 iterations, 250 burn-in iterations, one chain, and no method-ranking claims.", "",
  "## Worked user example", "",
  "`worked_finemapping_example.R` downloads the public simulated PLINK files, creates and QC-filters a Glist, simulates one phenotype, fits a real ST-BED BayesC model, and evaluates causal recovery. Its short settings demonstrate the workflow and are not convergence recommendations.", "",
  "## Developer contract smoke test", "",
  "`contract_smoke_test.R` quickly checks simulation, oracle, result-object, and metric contracts without fitting a sampler.", "",
  "## Complete benchmark reproduction", "",
  "`run_benchmark.R` runs the configuration-driven 10-replicate targets workflow. Downloaded qgdata files are cached and checksum-validated; `_targets/` should normally be retained so only missing or outdated work runs.", "",
  "## Public simulated data", "",
  paste0("The five inputs are publicly accessible from `", data_source$repository, "` at commit `", data_source$commit, "`. `example_data_manifest.csv` records pinned URLs, sizes, MD5 checksums, and roles. The cached benchmark files match that revision exactly."), "",
  "No explicit data licence was found in the qgdata repository root at the pinned revision. Public accessibility is verified, but reuse terms should be clarified with the repository owner.", "",
  "Exact numerical reproduction uses the pinned files, checksums, package versions, seeds, platform, compiler, and numerical libraries. Small numerical differences may occur across platforms.", "",
  "The capsule excludes genotype data, sparse LD, fit objects, posterior samples, and `_targets/`. Verify files against `checksums.csv`; that table intentionally excludes itself.")
writeLines(readme, file.path(snapshot_dir, "README.md"))

checksum_files <- setdiff(list.files(snapshot_dir), "checksums.csv")
info <- file.info(file.path(snapshot_dir, checksum_files))
checksums <- data.frame(file = checksum_files, size_bytes = info$size,
  md5 = unname(tools::md5sum(file.path(snapshot_dir, checksum_files))), row.names = NULL)
write.csv(checksums, file.path(snapshot_dir, "checksums.csv"), row.names = FALSE)
message("Created validated reference snapshot: ", snapshot_dir)
