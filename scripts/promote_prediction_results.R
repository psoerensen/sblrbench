# Promote reviewed Study 02 compact outputs after the benchmark source is committed.
status <- system2("git", c("status", "--porcelain"), stdout = TRUE)
if (length(status)) stop("Final prediction promotion requires a clean Git working tree.", call. = FALSE)
commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
source_dir <- file.path("results", "local", "02_prediction")
destination <- file.path("results", "reference", "02_prediction", "st-mt-bayesr-development-v1")
required <- c("prediction_metrics.csv", "paired_method_differences.csv", "computational_summary.csv",
  "replicate_status.csv", "simulation_summary.csv", "prediction_manifest.json",
  "method_summary.csv", "paired_method_summary.csv")
missing <- required[!file.exists(file.path(source_dir, required))]
if (length(missing)) stop("Missing local prediction outputs: ", paste(missing, collapse = ", "), call. = FALSE)
status_table <- utils::read.csv(file.path(source_dir, "replicate_status.csv"), stringsAsFactors = FALSE)
if (nrow(status_table) != 80L || any(status_table$status != "ok")) stop("Promotion requires 80 successful method branches.", call. = FALSE)
if (dir.exists(destination)) stop("Reference destination already exists: ", destination, call. = FALSE)
dir.create(destination, recursive = TRUE)
file.copy(file.path(source_dir, required), destination)
file.copy(c("studies/02_prediction/config.R", "scripts/worked_prediction_example.R",
  "scripts/prediction_contract_smoke_test.R", "scripts/run_prediction_benchmark.R"), destination)
writeLines(c("# ST/MT BayesR prediction development benchmark", "",
  paste("Source commit:", commit), "", "Development settings; not a final scientific comparison."),
  file.path(destination, "README.md"))
files <- setdiff(list.files(destination), "checksums.csv")
info <- file.info(file.path(destination, files))
checksums <- data.frame(file = files, size_bytes = info$size,
  md5 = unname(tools::md5sum(file.path(destination, files))))
utils::write.csv(checksums, file.path(destination, "checksums.csv"), row.names = FALSE)
message("Created clean reference capsule at ", destination)
