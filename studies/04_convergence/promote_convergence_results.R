source(file.path("studies", "04_convergence", "promotion.R"))
source_dir <- Sys.getenv("SBLR_BENCH_OUTPUT_DIR",
  file.path("results", "local", "04_convergence"))
destination <- Sys.getenv("SBLR_BENCH_CAPSULE_DESTINATION",
  file.path("results", "reference", "04_convergence", "current-selection"))
if (dir.exists(destination)) stop("Reference capsule already exists.", call. = FALSE)

result_files <- c("benchmark_manifest.json", "diagnostic_registry.csv",
  "diagnostic_thresholds.csv", "scalar_chain_draws.csv", "convergence_diagnostics.csv",
  "checkpoint_summary.csv", "burnin_stability.csv", "chain_summaries.csv",
  "chain_agreement.csv", "method_recommendations.csv", "computational_summary.csv",
  "chain_status.csv", "seed_registry.csv", "simulation_truth.csv")
if (any(!file.exists(file.path(source_dir, result_files))))
  stop("Local Study 04 results are incomplete.", call. = FALSE)

dir.create(destination, recursive = TRUE)
file.copy(file.path(source_dir, result_files), destination)
study_files <- c("config.R", "diagnostic_registry.R", "methods.R", "chain_extraction.R",
  "diagnostics.R", "recommendations.R", "run_convergence_benchmark.R",
  "convergence_contract_smoke_test.R", "worked_convergence_example.R", "pilot.R", "targets.R")
file.copy(file.path("studies", "04_convergence", study_files), destination)
config <- source(file.path("studies", "04_convergence", "config.R"), local = TRUE)$value

data_manifest <- data.frame(repository = config$example_data$repository,
  commit = config$example_data$commit,
  path = file.path(config$example_data$subdirectory, config$example_data$files),
  filename = config$example_data$files,
  size_bytes = unname(config$example_data$size_bytes), md5 = unname(config$example_data$md5),
  download_url = paste0("https://raw.githubusercontent.com/", config$example_data$repository,
    "/", config$example_data$commit, "/", config$example_data$subdirectory, "/", config$example_data$files),
  role = c("PLINK genotype", "PLINK marker metadata", "PLINK sample metadata", "example phenotype", "example covariates"))
utils::write.csv(data_manifest, file.path(destination, "example_data_manifest.csv"), row.names = FALSE)

source_paths <- file.path("studies", "04_convergence", c(study_files, "promotion.R", "promote_convergence_results.R"))
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
source_files <- data.frame(path = source_paths,
  role = c("configuration", "diagnostic registry", "method execution", "chain extraction",
    "diagnostics", "recommendations", "entry point", "contract smoke test", "worked example",
    "study helpers", "targets pipeline", "capsule validation", "capsule promotion"),
  git_blob_sha = vapply(source_paths, function(x) system2("git", c("hash-object", x), stdout = TRUE), character(1)),
  md5 = unname(tools::md5sum(source_paths)), source_commit = commit,
  source_url = paste0("https://github.com/psoerensen/sblrbench/blob/", commit, "/", source_paths),
  stringsAsFactors = FALSE)
utils::write.csv(source_files, file.path(destination, "source_files.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(destination, "session_info.txt"))
writeLines(c("# Single-trait multichain convergence development benchmark", "",
  "This capsule contains one matched simulation per architecture, four methods, four identifiable",
  "chains per method, and 3,000 unthinned raw draws per chain. It assesses operational multichain",
  "diagnostics and provisional run-length recommendations. It does not assess prediction, parameter",
  "accuracy, model validity, or scientific superiority. Passing thresholds on these datasets does not",
  "prove universal convergence. Public qgdata inputs are pinned and checksum-validated; genotypes and",
  "native fits are not included. No explicit data licence was identified; reuse terms should be clarified."),
  file.path(destination, "README.md"))
files <- setdiff(list.files(destination, recursive = TRUE), "checksums.csv")
checks <- data.frame(file = files, size_bytes = file.info(file.path(destination, files))$size,
  md5 = unname(.study04_canonical_md5(file.path(destination, files))))
utils::write.csv(checks, file.path(destination, "checksums.csv"), row.names = FALSE)
.study04_validate_capsule(destination)
cat(normalizePath(destination), "\n")
