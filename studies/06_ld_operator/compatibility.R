.study06_current_fit_sha <- "bd8e2c8148a0d9540dc20716455706beeb16fa86"
.study06_refresh_sha <- "02e8c74baa906e83c4a08d42a9cc6339b4e81072"

.study06_git_blob <- function(repository, revision, path) {
  out <- system2("git", c("-C", repository, "rev-parse",
    paste0(revision, ":", path)), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("Cannot resolve ", path, " at ", revision, ": ",
      paste(out, collapse = "\n"), call. = FALSE)
  }
  trimws(out[[1L]])
}

.study06_source_compatibility <- function(
    sibling = file.path("..", "sblr"),
    fit_sha = .study06_current_fit_sha,
    refresh_sha = .study06_refresh_sha) {
  paths <- c(
    "R/stblr-public.R",
    "R/stblr-block-eigen.R",
    "R/RcppExports.R",
    "src/st_cpg_omp_individual_scheduled_chains.cpp",
    "src/stblr_cpg_omp_bed_marker_scheduled_chains_bayesr.cpp",
    "src/st_cpg_omp_csr.cpp",
    "src/st_cpg_omp_csr_bayesr.cpp",
    "src/st_block_eigen.cpp",
    "src/st_block_low_rank.cpp",
    "src/st_ld_operator.h",
    "src/blr_block_low_rank.h",
    "src/st_csr_common.h",
    "src/st_bed_decode.h",
    "src/st_bed_bayesr_common.h")
  role <- c(rep("R dispatch and native binding", 3L),
    rep("compiled sampler transition", 4L),
    rep("compiled retained-low-rank operator", 4L),
    rep("compiled shared sampler contract", 3L))
  old <- vapply(paths, function(path)
    .study06_git_blob(sibling, fit_sha, path), "")
  current <- vapply(paths, function(path)
    .study06_git_blob(sibling, refresh_sha, path), "")
  data.frame(path = paths, role = role, fit_sha_blob = old,
    refresh_sha_blob = current, identical = old == current,
    stringsAsFactors = FALSE)
}

.study06_prior_compatibility <- function(m = 37991L, h2 = 0.30,
                                         nub = 4, nue = 4) {
  configurations <- c("bed", "full_csr", "block_csr", "low_rank_full",
    "low_rank_0999", "low_rank_0995")
  rows <- do.call(rbind, lapply(c("bayesc", "bayesr"), function(model) {
    if (model == "bayesc") {
      initial_probability <- 0.01
      prior_probability <- 0.001
      initial_multiplier <- initial_probability
      prior_multiplier <- prior_probability
      mixture <- "1"
    } else {
      probability <- c(0.99, rep(0.01 / 3, 3))
      alpha <- probability * 5e5
      multipliers <- c(0, 0.01, 0.1, 1)
      initial_probability <- sum(probability[-1L])
      prior_probability <- sum((alpha / sum(alpha))[-1L])
      initial_multiplier <- sum(probability * multipliers)
      prior_multiplier <- sum((alpha / sum(alpha)) * multipliers)
      mixture <- paste(format(multipliers, trim = TRUE), collapse = ";")
    }
    data.frame(model = model, configuration = configurations,
      marker_count = m, h2 = h2, nub = nub, nue = nue,
      marker_scale = 1,
      initial_inclusion_weight = initial_probability,
      prior_mean_inclusion_weight = prior_probability,
      initial_mixture_multiplier_weight = initial_multiplier,
      prior_mean_mixture_multiplier_weight = prior_multiplier,
      mixture_multipliers = mixture,
      old_B_per_phenotypic_variance = h2 / (m * initial_multiplier),
      current_B_per_phenotypic_variance = h2 / (m * initial_multiplier),
      old_ssb_prior_per_phenotypic_variance =
        ((nub - 2) / nub) * h2 / (m * prior_multiplier),
      current_ssb_prior_per_phenotypic_variance =
        ((nub - 2) / nub) * h2 / (m * prior_multiplier),
      old_E_per_phenotypic_variance = 1 - h2,
      current_E_per_phenotypic_variance = 1 - h2,
      old_sse_prior_per_phenotypic_variance =
        ((nue - 2) / nue) * (1 - h2),
      current_sse_prior_per_phenotypic_variance =
        ((nue - 2) / nue) * (1 - h2),
      annotation_weight = NA_real_,
      identical_resolved_prior = TRUE, stringsAsFactors = FALSE)
  }))
  rownames(rows) <- NULL
  rows
}

.study06_write_current_compatibility <- function(output,
    sibling = file.path("..", "sblr"),
    fit_sha = .study06_current_fit_sha,
    refresh_sha = .study06_refresh_sha) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  source <- .study06_source_compatibility(sibling, fit_sha, refresh_sha)
  prior <- .study06_prior_compatibility()
  if (!all(source$identical) || !all(prior$identical_resolved_prior)) {
    stop("Study 06 is not compatible with the current sblr refresh.",
      call. = FALSE)
  }
  write.csv(source, file.path(output, "source_compatibility.csv"),
    row.names = FALSE, na = "")
  write.csv(prior, file.path(output, "prior_compatibility.csv"),
    row.names = FALSE, na = "")
  manifest <- list(
    validation_status = "compatible_without_sampler_rerun",
    fit_source_sha = fit_sha,
    refresh_validation_sha = refresh_sha,
    package_version = "0.2.0",
    numerical_fit_provenance_preserved = TRUE,
    rerun_required = FALSE,
    sampler_transition_sources_identical = TRUE,
    retained_low_rank_sources_identical = TRUE,
    resolved_priors_identical = TRUE,
    annotation_models_in_scope = FALSE,
    sampled_maf_s_in_scope = FALSE,
    explanation = paste(
      "The refresh commit centralizes scalar prior calibration and records",
      "calibration metadata. Under the frozen Study 06 BayesC and BayesR",
      "settings, the resolved B, E, ssb_prior, sse_prior, inclusion weights,",
      "component weights, mixture multipliers, and marker scale are unchanged.",
      "All compiled scalar sampler and retained-low-rank transition sources",
      "are byte-identical. The 60 numerical fits retain their original SHA."))
  jsonlite::write_json(manifest,
    file.path(output, "compatibility_manifest.json"), auto_unbox = TRUE,
    pretty = TRUE)
  invisible(manifest)
}

.study06_current_checksums <- function(path) {
  files <- sort(list.files(path, recursive = TRUE, full.names = TRUE))
  files <- files[basename(files) != "checksums.csv"]
  relative <- substring(normalizePath(files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(path, winslash = "/", mustWork = TRUE)) + 2L)
  data.frame(file = relative, size_bytes = file.info(files)$size,
    md5 = unname(tools::md5sum(files)), stringsAsFactors = FALSE)
}

.study06_normalize_capsule_text <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE)
  files <- files[tolower(tools::file_ext(files)) %in%
    c("csv", "json", "md", "r", "txt")]
  for (file in files) {
    size <- file.info(file)$size
    con <- file(file, open = "rb")
    bytes <- readBin(con, what = "raw", n = size)
    close(con)
    text <- rawToChar(bytes)
    normalized <- gsub("\r\n", "\n", text, fixed = TRUE)
    normalized <- gsub("[ \t]+(?=\n|$)", "", normalized, perl = TRUE)
    if (!identical(text, normalized)) {
      con <- file(file, open = "wb")
      writeBin(charToRaw(normalized), con)
      close(con)
    }
  }
  invisible(files)
}

.study06_validate_current_capsule <- function(path,
    fit_sha = .study06_current_fit_sha,
    refresh_sha = .study06_refresh_sha) {
  required <- c("benchmark_manifest.json", "compatibility_manifest.json",
    "prior_compatibility.csv", "source_compatibility.csv",
    "checkpoint_validation.csv", "prediction_metrics.csv",
    "parameter_estimates.csv", "paired_replicate_differences.csv",
    "method_recommendations.csv", "readiness_decision.txt",
    "session_info.txt", "checksums.csv")
  if (!dir.exists(path) || !all(file.exists(file.path(path, required)))) {
    stop("The current Study 06 capsule is incomplete.", call. = FALSE)
  }
  manifest <- jsonlite::read_json(file.path(path,
    "benchmark_manifest.json"), simplifyVector = TRUE)
  compatibility <- jsonlite::read_json(file.path(path,
    "compatibility_manifest.json"), simplifyVector = TRUE)
  if (!identical(manifest$sblr_source_sha, fit_sha) ||
      !identical(manifest$refresh_validation_sha, refresh_sha) ||
      !identical(manifest$numerical_fit_source_sha, fit_sha) ||
      !identical(compatibility$validation_status,
        "compatible_without_sampler_rerun") ||
      !identical(compatibility$fit_source_sha, fit_sha) ||
      !identical(compatibility$refresh_validation_sha, refresh_sha) ||
      !identical(compatibility$rerun_required, FALSE)) {
    stop("Study 06 compatibility provenance is inconsistent.", call. = FALSE)
  }
  checkpoint <- read.csv(file.path(path, "checkpoint_validation.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  expected <- expand.grid(architecture = c("sparse_homogeneous",
      "sparse_mixture"), replicate = 1:5,
    configuration = c("bed", "full_csr", "block_csr", "low_rank_full",
      "low_rank_0999", "low_rank_0995"), stringsAsFactors = FALSE)
  keys <- function(x) paste(x$architecture, x$replicate, x$configuration,
    sep = "::")
  if (nrow(checkpoint) != 60L || anyDuplicated(keys(checkpoint)) ||
      !setequal(keys(checkpoint), keys(expected)) ||
      !all(checkpoint$package_commit == fit_sha) ||
      !all(checkpoint$status == "ok")) {
    stop("The current Study 06 fit grid is invalid.", call. = FALSE)
  }
  prior <- read.csv(file.path(path, "prior_compatibility.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  source <- read.csv(file.path(path, "source_compatibility.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(prior) != 12L || !all(prior$identical_resolved_prior) ||
      !all(source$identical)) {
    stop("Study 06 numerical compatibility evidence failed.", call. = FALSE)
  }
  checksums <- read.csv(file.path(path, "checksums.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  observed <- .study06_current_checksums(path)
  if (!identical(checksums$file, observed$file) ||
      !identical(as.numeric(checksums$size_bytes),
        as.numeric(observed$size_bytes)) ||
      !identical(checksums$md5, observed$md5)) {
    stop("Study 06 current capsule checksums failed.", call. = FALSE)
  }
  TRUE
}

.study06_promote_current_capsule <- function(
    source = file.path("results", "reference", "06_ld_operator",
      "st-low-rank-operator-five-replicate-development-v2"),
    compatibility = file.path("results", "local",
      "current_benchmark_refresh", "study06_compatibility"),
    destination = file.path("results", "reference", "06_ld_operator",
      "current")) {
  if (!dir.exists(source)) stop("Study 06 v2 source capsule is missing.",
    call. = FALSE)
  if (dir.exists(destination)) stop("Study 06 current capsule already exists.",
    call. = FALSE)
  staging <- paste0(destination, ".staging")
  if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(source, full.names = TRUE), staging,
    recursive = TRUE)
  if (!all(ok)) stop("Cannot stage the Study 06 source capsule.",
    call. = FALSE)
  staged <- staging
  for (file in c("compatibility_manifest.json", "prior_compatibility.csv",
      "source_compatibility.csv")) {
    if (!file.copy(file.path(compatibility, file), file.path(staged, file))) {
      stop("Cannot add Study 06 compatibility artifact: ", file,
        call. = FALSE)
    }
  }
  manifest_path <- file.path(staged, "benchmark_manifest.json")
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  manifest$benchmark_scope <- "current_retained_low_rank_operator_benchmark"
  manifest$numerical_fit_source_sha <- .study06_current_fit_sha
  manifest$refresh_validation_sha <- .study06_refresh_sha
  manifest$compatibility_status <- "compatible_without_sampler_rerun"
  manifest$provenance_note <- paste(
    "The 60 numerical fits were generated at sblr",
    .study06_current_fit_sha,
    "and retain that provenance. Deterministic validation at",
    .study06_refresh_sha,
    "established identical resolved priors and byte-identical scalar sampler",
    "and retained-low-rank transition sources.")
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  writeLines(c(
    "# Study 06: current retained low-rank LD operator benchmark",
    "",
    "This current capsule contains 60 validated four-chain fits across two",
    "architectures, five paired replicates, and six operator configurations.",
    "The numerical fits retain their original sblr source SHA",
    paste0("`", .study06_current_fit_sha, "`."),
    "",
    paste("Compatibility with the current refresh SHA",
      paste0("`", .study06_refresh_sha, "`"),
      "was established without rerunning samplers: resolved BayesC and BayesR",
      "priors are identical, and all compiled scalar sampler and retained-low-",
      "rank transition sources are byte-identical. See",
      "`compatibility_manifest.json`, `prior_compatibility.csv`, and",
      "`source_compatibility.csv`."),
    "",
    "Native fits, checkpoints, Q matrices, and target stores are excluded."),
    file.path(staged, "README.md"))
  .study06_normalize_capsule_text(staged)
  unlink(file.path(staged, "checksums.csv"), force = TRUE)
  checksums <- .study06_current_checksums(staged)
  write.csv(checksums, file.path(staged, "checksums.csv"), row.names = FALSE,
    na = "")
  .study06_normalize_capsule_text(staged)
  .study06_validate_current_capsule(staged)
  if (!file.rename(staged, destination)) {
    stop("Atomic Study 06 current-capsule promotion failed.", call. = FALSE)
  }
  unlink(staging, recursive = TRUE, force = TRUE)
  .study06_validate_current_capsule(destination)
}
