source(file.path("studies", "03_parameter_estimation", "promotion.R"))
source_dir <- file.path("results", "local", "03_parameter_estimation")
destination <- file.path("results", "reference", "03_parameter_estimation",
  "st-parameter-estimation-one-replicate-development-v1")
if (dir.exists(destination)) stop("Reference capsule already exists.", call. = FALSE)
dir.create(destination, recursive = TRUE)
result_files <- c("estimand_registry.csv", "simulation_truth.csv", "parameter_estimates.csv",
  "parameter_recovery_summary.csv", "paired_parameter_differences.csv",
  "paired_comparison_summary.csv", "computational_summary.csv", "replicate_status.csv", "seed_registry.csv")
if (any(!file.exists(file.path(source_dir, result_files)))) stop("Local Study 03 results are incomplete.", call. = FALSE)
file.copy(file.path(source_dir, result_files), destination)
study_files <- c("config.R", "estimands.R", "simulation.R", "methods.R", "metrics.R", "pilot.R", "targets.R",
  "run_parameter_estimation_benchmark.R", "parameter_estimation_contract_smoke_test.R", "worked_parameter_estimation_example.R")
file.copy(file.path("studies", "03_parameter_estimation", study_files), destination)
config <- source(file.path("studies", "03_parameter_estimation", "config.R"), local = TRUE)$value
data_manifest <- data.frame(repository = config$example_data$repository, commit = config$example_data$commit,
  path = file.path(config$example_data$subdirectory, config$example_data$files), filename = config$example_data$files,
  size_bytes = unname(config$example_data$size_bytes), md5 = unname(config$example_data$md5),
  download_url = paste0("https://raw.githubusercontent.com/", config$example_data$repository, "/",
    config$example_data$commit, "/", config$example_data$subdirectory, "/", config$example_data$files),
  role = c("PLINK genotype", "PLINK marker metadata", "PLINK sample metadata", "example phenotype", "example covariates"))
utils::write.csv(data_manifest, file.path(destination, "example_data_manifest.csv"), row.names = FALSE)
status <- read.csv(file.path(source_dir, "replicate_status.csv")); truth <- read.csv(file.path(source_dir, "simulation_truth.csv"))
manifest <- list(study_id = config$study, task = config$task, benchmark_scope = "one_replicate_development",
  benchmark_status = "complete", development_settings = TRUE,
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  repository_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  source_tree_clean = !length(system2("git", c("status", "--porcelain"), stdout = TRUE)),
  sblr_commit = system2("git", c("-C", "../sblr", "rev-parse", "HEAD"), stdout = TRUE),
  qgdata_commit = config$example_data$commit, active_methods = config$methods,
  architectures = names(config$simulation$architectures), replicate_count = 1L,
  expected_fit_count = 8L, successful_fit_count = sum(status$status == "ok"),
  failed_fit_count = sum(status$status != "ok"), analysis_sample_count = 5000L,
  canonical_marker_count = 37991L, causal_marker_count = config$simulation$n_causal,
  target_h2 = config$simulation$h2, realized_h2 = truth[truth$estimand_id == "heritability", c("architecture", "truth")],
  mcmc = config$profiles$development, seeds = read.csv(file.path(source_dir, "seed_registry.csv")),
  estimands = read.csv(file.path(source_dir, "estimand_registry.csv")),
  primary_estimands = read.csv(file.path(source_dir, "estimand_registry.csv"))$estimand_id[read.csv(file.path(source_dir, "estimand_registry.csv"))$primary],
  truth_definitions = unique(truth[c("estimand_id", "truth_type", "truth_definition")]),
  model_specification_rules = c("homogeneous BayesC matched", "mixture BayesR matched"),
  data_provenance = config$example_data)
jsonlite::write_json(manifest, file.path(destination, "benchmark_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
sources <- file.path("studies", "03_parameter_estimation", study_files)
source_inventory <- data.frame(path = sources, role = c("configuration", "estimand registry", "simulation and truth",
  "method and posterior extraction", "recovery metrics", "study helpers", "targets pipeline", "benchmark entry point",
  "contract smoke test", "worked example"), md5 = unname(tools::md5sum(sources)),
  source_commit = manifest$repository_commit, stringsAsFactors = FALSE)
utils::write.csv(source_inventory, file.path(destination, "source_files.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(destination, "session_info.txt"))
writeLines(c("# Single-trait parameter-estimation development benchmark", "",
  "This complete one-replicate development capsule validates eight fits, posterior extraction,",
  "truth-aware recovery summaries and paired BED/CSR comparisons. It does not establish bias,",
  "interval calibration, convergence, or method rankings. Public qgdata inputs are pinned and",
  "checksum-validated; genotype data, native fits and posterior draws are not included."),
  file.path(destination, "README.md"))
files <- setdiff(list.files(destination, recursive = TRUE), "checksums.csv")
checks <- data.frame(file = files, size_bytes = file.info(file.path(destination, files))$size,
  md5 = unname(.study03_canonical_md5(file.path(destination, files))))
utils::write.csv(checks, file.path(destination, "checksums.csv"), row.names = FALSE)
.study03_validate_capsule(destination)
cat(normalizePath(destination), "\n")
