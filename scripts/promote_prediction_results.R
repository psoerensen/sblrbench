# Promote the complete one-replicate Study 02 development benchmark.
source(file.path("studies", "02_prediction", "promotion.R"), local = TRUE)

git_status <- system2("git", c("status", "--porcelain"), stdout = TRUE)
if (length(git_status)) stop("Prediction promotion requires a clean Git working tree.", call. = FALSE)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
config <- source(file.path("studies", "02_prediction", "config.R"), local = TRUE)$value
profile <- config$reference_profiles$one_replicate_development
if (profile$replicate_count != 1L) stop("The reference profile must contain one replicate.", call. = FALSE)

source_files <- c(
  "studies/02_prediction/config.R", "studies/02_prediction/pilot.R",
  "studies/02_prediction/targets.R", "studies/02_prediction/promotion.R",
  "R/prediction.R", "scripts/worked_prediction_example.R",
  "scripts/prediction_contract_smoke_test.R", "scripts/run_prediction_benchmark.R",
  "scripts/promote_prediction_results.R"
)
roles <- c("study configuration", "study helpers", "targets pipeline",
  "promotion validation", "prediction contracts", "worked example",
  "developer smoke test", "benchmark entry point", "promotion entry point")
published <- vapply(source_files, function(path)
  system2("git", c("cat-file", "-e", paste0("HEAD:", path))) == 0L, logical(1))
if (!all(published)) stop("Required source files are absent from HEAD: ",
  paste(source_files[!published], collapse = ", "), call. = FALSE)

source_dir <- file.path("results", "local", "02_prediction")
destination <- file.path("results", "reference", "02_prediction", profile$capsule_id)
result_files <- c("prediction_metrics.csv", "paired_method_differences.csv",
  "computational_summary.csv", "replicate_status.csv", "simulation_summary.csv",
  "benchmark_summary.csv", "paired_comparison_summary.csv", "prediction_manifest.json")
missing <- result_files[!file.exists(file.path(source_dir, result_files))]
if (length(missing)) stop("Missing local prediction outputs: ", paste(missing, collapse = ", "), call. = FALSE)

metrics <- read.csv(file.path(source_dir, "prediction_metrics.csv"))
paired <- read.csv(file.path(source_dir, "paired_method_differences.csv"))
computation <- read.csv(file.path(source_dir, "computational_summary.csv"))
status <- read.csv(file.path(source_dir, "replicate_status.csv"))
simulations <- read.csv(file.path(source_dir, "simulation_summary.csv"))
benchmark_summary <- read.csv(file.path(source_dir, "benchmark_summary.csv"))
paired_summary <- read.csv(file.path(source_dir, "paired_comparison_summary.csv"))
manifest <- jsonlite::read_json(file.path(source_dir, "prediction_manifest.json"), simplifyVector = TRUE)
if (!identical(manifest$sblrbench_source_commit, commit) || !isTRUE(manifest$source_tree_clean))
  stop("Manifest source provenance must equal clean HEAD.", call. = FALSE)
.study02_validate_promotion_tables(config, metrics, paired, computation, status,
  simulations, manifest, benchmark_summary, paired_summary)

if (dir.exists(destination)) stop("Reference destination already exists: ", destination, call. = FALSE)
dir.create(destination, recursive = TRUE)
ok <- file.copy(file.path(source_dir, result_files), destination)
if (!all(ok)) stop("Failed to copy compact prediction outputs.", call. = FALSE)
file.rename(file.path(destination, "prediction_manifest.json"),
  file.path(destination, "benchmark_manifest.json"))

file.copy("studies/02_prediction/config.R", file.path(destination, "config.R"))
file.copy("studies/02_prediction/pilot.R", file.path(destination, "pilot.R"))
file.copy("studies/02_prediction/targets.R", file.path(destination, "targets.R"))
file.copy("scripts/run_prediction_benchmark.R", file.path(destination, "run_prediction_benchmark.R"))
file.copy("scripts/worked_prediction_example.R", file.path(destination, "worked_prediction_example.R"))
file.copy("scripts/prediction_contract_smoke_test.R", file.path(destination, "prediction_contract_smoke_test.R"))

data_files <- config$example_data$files
data_manifest <- data.frame(repository = config$example_data$repository,
  commit = config$example_data$commit,
  path = file.path(config$example_data$subdirectory, data_files), filename = data_files,
  size_bytes = unname(config$example_data$size_bytes[data_files]),
  md5 = unname(config$example_data$md5[data_files]),
  download_url = paste0("https://raw.githubusercontent.com/", config$example_data$repository,
    "/", config$example_data$commit, "/", config$example_data$subdirectory, "/", data_files),
  role = c("PLINK genotype", "PLINK marker metadata", "PLINK sample metadata",
    "example phenotype", "example covariates"))
write.csv(data_manifest, file.path(destination, "example_data_manifest.csv"), row.names = FALSE)

blobs <- vapply(source_files, function(path)
  system2("git", c("rev-parse", paste0(commit, ":", path)), stdout = TRUE), character(1))
source_inventory <- data.frame(path = source_files, role = roles,
  git_blob_sha = blobs, md5 = unname(tools::md5sum(source_files)),
  source_commit = commit,
  source_url = paste0("https://github.com/psoerensen/sblrbench/blob/", commit, "/", source_files))
write.csv(source_inventory, file.path(destination, "source_files.csv"), row.names = FALSE)

session <- c(capture.output(sessionInfo()), "",
  paste("sblr:", packageVersion("sblr")), paste("sblrbench:", packageVersion("sblrbench")),
  paste("qgg:", packageVersion("qgg")), paste("targets:", packageVersion("targets")),
  paste("source commit:", commit))
writeLines(session, file.path(destination, "session_info.txt"))

readme <- c("# Single-trait prediction one-replicate development benchmark", "",
  "## Benchmark status", "", "This is a complete full-size one-replicate development benchmark.", "",
  "## Scientific question", "", "BayesC versus BayesR and BED versus CSR for held-out single-trait prediction.", "",
  "## What is complete", "", "Two architectures, one replicate per architecture, four methods, eight successful fits, complete metrics and paired comparisons, and passing leakage and oracle checks.", "",
  "## What is not claimed", "", "This capsule provides no replicate-to-replicate uncertainty, convergence evidence, final method ranking, general runtime superiority, or comparison beyond single-trait models.", "",
  "## Reproduction", "", paste("The public simulated qgdata files are pinned to commit", config$example_data$commit, "and checksum validated."),
  paste("The clean analysis source commit is", commit, "."), "Run `run_prediction_benchmark.R` from a compatible clone. Valid targets are reused; development settings are short and single-chain. Small numerical differences may occur across platforms and numerical libraries.", "",
  "No explicit data licence was identified for qgdata; clarify reuse terms before redistribution.")
writeLines(readme, file.path(destination, "README.md"))

text_files <- list.files(destination, full.names = TRUE)
contents <- unlist(lapply(text_files, function(path) readLines(path, warn = FALSE)), use.names = FALSE)
if (any(grepl("[A-Za-z]:[/\\\\]Users[/\\\\]", contents))) stop("Absolute user path detected.", call. = FALSE)
files <- sort(setdiff(list.files(destination), "checksums.csv"))
info <- file.info(file.path(destination, files))
checksums <- data.frame(file = files, size_bytes = info$size,
  md5 = unname(tools::md5sum(file.path(destination, files))))
write.csv(checksums, file.path(destination, "checksums.csv"), row.names = FALSE)
.study02_validate_capsule(destination)
message("Created clean reference capsule at ", destination)
