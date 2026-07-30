# Promote a reviewed, complete Study 02 single-trait benchmark.
source(file.path("studies", "02_prediction", "promotion.R"), local = TRUE)

git_status <- system2("git", c("status", "--porcelain"), stdout = TRUE)
if (length(git_status)) stop("Final prediction promotion requires a clean Git working tree.", call. = FALSE)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
source_files <- c("studies/02_prediction/config.R",
  "studies/02_prediction/pilot.R", "studies/02_prediction/targets.R",
  "scripts/worked_prediction_example.R",
  "scripts/prediction_contract_smoke_test.R",
  "scripts/run_prediction_benchmark.R")
published <- vapply(source_files, function(path)
  system2("git", c("cat-file", "-e", paste0("HEAD:", path))) == 0L, logical(1))
if (!all(published)) stop("The clean source commit does not contain all revised Study 02 files.", call. = FALSE)

source_dir <- file.path("results", "local", "02_prediction")
destination <- file.path("results", "reference", "02_prediction",
  "st-bayesc-bayesr-development-v1")
result_files <- c("prediction_metrics.csv", "paired_method_differences.csv",
  "computational_summary.csv", "replicate_status.csv", "simulation_summary.csv",
  "prediction_manifest.json")
missing <- result_files[!file.exists(file.path(source_dir, result_files))]
if (length(missing)) stop("Missing local prediction outputs: ",
  paste(missing, collapse = ", "), call. = FALSE)
config <- source("studies/02_prediction/config.R", local = TRUE)$value
metrics <- read.csv(file.path(source_dir, "prediction_metrics.csv"))
paired <- read.csv(file.path(source_dir, "paired_method_differences.csv"))
computation <- read.csv(file.path(source_dir, "computational_summary.csv"))
status <- read.csv(file.path(source_dir, "replicate_status.csv"))
simulations <- read.csv(file.path(source_dir, "simulation_summary.csv"))
manifest <- jsonlite::read_json(file.path(source_dir, "prediction_manifest.json"),
  simplifyVector = TRUE)
.study02_validate_promotion_tables(config, metrics, paired, computation,
  status, simulations, manifest)
if (dir.exists(destination)) stop("Reference destination already exists: ", destination, call. = FALSE)
dir.create(destination, recursive = TRUE)
file.copy(file.path(source_dir, result_files), destination)
file.copy(source_files, destination)
writeLines(c("# Single-trait prediction development benchmark", "",
  paste("Source commit:", commit), "",
  "Four active ST BayesC/BayesR BED/CSR methods; multi-trait prediction is deferred.",
  "Development settings are not a final scientific comparison."),
  file.path(destination, "README.md"))
files <- setdiff(list.files(destination), "checksums.csv")
info <- file.info(file.path(destination, files))
checksums <- data.frame(file = files, size_bytes = info$size,
  md5 = unname(tools::md5sum(file.path(destination, files))))
write.csv(checksums, file.path(destination, "checksums.csv"), row.names = FALSE)
message("Created clean reference capsule at ", destination)
