.study06v2_sha256 <- function(path) digest::digest(file = path,
  algo = "sha256", serialize = FALSE)

.study06v2_summary_statistics_hash <- function(store = file.path("results",
    "local", "study06_ld_operator", "_targets")) {
  x <- targets::tar_read_raw("study06_operator_stats", store = store)
  digest::digest(list(marker_names = x$marker_names, rows = x$rows,
    af = x$af, wy = x$wy, ww = x$ww, yy = x$yy, n = x$n),
    algo = "sha256")
}

.study06v2_repin_deterministic <- function(config,
    old_sha = "651d8637d8b04b6f6d049a2e8ab31c11d8318dc6",
    old_snapshot = file.path("results", "local",
      "study06_low_rank_operator_v2", "preflight", "sblr_source_651d8637")) {
  output <- file.path(config$local_dir, "deterministic")
  pre_path <- file.path(output, "pre_repin_hash_inventory.csv")
  if (!file.exists(pre_path)) stop("Preserved pre-repin inventory is absent.",
    call. = FALSE)
  pre <- read.csv(pre_path, stringsAsFactors = FALSE)
  relevant_source <- pre$path[pre$category == "low_rank_source"]
  rows <- list(); add <- function(category, path, old_hash, new_hash,
                                  relevance, note = "") {
    rows[[length(rows) + 1L]] <<- data.frame(category = category,
      path = path, old_sha256 = old_hash, new_sha256 = new_hash,
      unchanged = identical(old_hash, new_hash), relevance = relevance,
      note = note, stringsAsFactors = FALSE)
  }
  for (path in relevant_source) add("low_rank_source", path,
    pre$sha256[pre$category == "low_rank_source" & pre$path == path],
    .study06v2_sha256(file.path(config$source_snapshot, path)),
    "deterministic_operator_implementation")
  scripts <- pre$path[pre$category == "deterministic_script"]
  for (path in scripts) {
    old_hash <- pre$sha256[pre$category == "deterministic_script" &
      pre$path == path]
    if (basename(path) == "config.R") {
      raw <- readBin(path, "raw", n = file.info(path)$size)
      normalized <- gsub(config$required_sblr_sha, old_sha,
        rawToChar(raw), fixed = TRUE)
      normalized <- gsub("sblr_source_96487b3", "sblr_source_651d8637",
        normalized, fixed = TRUE)
      new_hash <- digest::digest(charToRaw(normalized), algo = "sha256",
        serialize = FALSE)
      note <- "byte-normalized for the two provenance-only pin substitutions"
    } else {
      new_hash <- .study06v2_sha256(path); note <- ""
    }
    add("deterministic_script", path, old_hash, new_hash,
      "deterministic_execution", note)
  }
  input_paths <- pre$path[pre$category %in% c("genotype_input",
    "block_definitions")]
  for (path in input_paths) add(pre$category[match(path, pre$path)], path,
    pre$sha256[match(path, pre$path)], .study06v2_sha256(path),
    "deterministic_input")
  old_stats <- pre$sha256[pre$category == "summary_statistics"]
  add("summary_statistics", "targets:study06_operator_stats", old_stats,
    .study06v2_summary_statistics_hash(), "deterministic_input")
  comparison <- do.call(rbind, rows)
  operator_paths <- c("src/st_block_eigen.cpp", "src/st_block_low_rank.cpp",
    "src/st_cpg_omp_csr.cpp", "src/st_cpg_omp_csr_bayesr.cpp",
    "R/stblr-block-eigen.R", "R/RcppExports.R")
  reuse <- all(comparison$unchanged) &&
    all(operator_paths %in% comparison$path)
  commit_files <- system2("git", c("-C", "../sblr", "diff-tree",
    "--no-commit-id", "--name-only", "-r", config$required_sblr_sha),
    stdout = TRUE)
  allowed_commit_files <- c("NEWS.md", "R/sparse_ld_bed_helper.R")
  reuse <- reuse && setequal(commit_files, allowed_commit_files)
  .study06v2_write_csv(comparison,
    file.path(output, "repin_hash_comparison.csv"))
  manifest <- list(
    validation_status = if (reuse) "passed_reuse_authorized" else
      "failed_rerun_required",
    old_source_sha = old_sha, new_source_sha = config$required_sblr_sha,
    package_version = config$required_sblr_version,
    commit_subject = "Correct BayesR prior variance calibration",
    changed_files = commit_files,
    retained_low_rank_sources_unchanged = all(comparison$unchanged[
      comparison$category == "low_rank_source"]),
    deterministic_scripts_semantically_unchanged = all(
      comparison$unchanged[comparison$category == "deterministic_script"]),
    deterministic_inputs_unchanged = all(comparison$unchanged[
      comparison$category %in% c("genotype_input", "block_definitions",
        "summary_statistics")]),
    deterministic_evidence_reused = reuse,
    numerical_summary_sha256 = .study06v2_sha256(file.path(output,
      "deterministic_identity_summary.csv")),
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  jsonlite::write_json(manifest,
    file.path(output, "deterministic_validation_manifest.json"),
    auto_unbox = TRUE, pretty = TRUE)
  if (!reuse) return(FALSE)
  checkpoint <- file.path(output, "near_full_operator_checkpoint.rds")
  object <- readRDS(checkpoint)
  if (!identical(object$source_sha, old_sha) ||
      !identical(object$package_version, config$required_sblr_version))
    stop("Preserved deterministic checkpoint has unexpected provenance.",
      call. = FALSE)
  object$source_sha <- config$required_sblr_sha
  object$repin_manifest <- "deterministic_validation_manifest.json"
  object$repin_reuse_authorized <- TRUE
  .study06v2_atomic_rds(object, checkpoint)
  TRUE
}
