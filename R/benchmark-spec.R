# Study specifications are ordinary lists. Validation is deliberately scoped to
# fields required by the currently supported prediction task.

.benchmark_scalar_string <- function(x, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop(field, " must be one non-empty string.", call. = FALSE)
  x
}

#' Read a benchmark specification
#'
#' @param path Path to an R file that defines an ordinary list named `spec`.
#' @return A validated benchmark specification list.
#' @export
read_benchmark_spec <- function(path) {
  .benchmark_scalar_string(path, "path")
  if (!file.exists(path))
    stop("Benchmark specification does not exist: ", path, call. = FALSE)
  environment <- new.env(parent = baseenv())
  sys.source(path, envir = environment)
  if (!exists("spec", envir = environment, inherits = FALSE))
    stop("Benchmark specification must define an object named `spec`.",
      call. = FALSE)
  spec <- get("spec", envir = environment, inherits = FALSE)
  validate_benchmark_spec(spec)
  spec_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- dirname(spec_path)
  repeat {
    if (file.exists(file.path(root, "DESCRIPTION"))) break
    parent <- dirname(root)
    if (identical(parent, root)) {
      root <- NULL
      break
    }
    root <- parent
  }
  attr(spec, "repository_root") <- root
  spec
}

benchmark_spec_path <- function(spec, path) {
  if (file.exists(path) || grepl("^[A-Za-z]:[/\\\\]|^/", path)) return(path)
  root <- attr(spec, "repository_root", exact = TRUE)
  if (is.null(root)) path else file.path(root, path)
}

#' Validate a benchmark specification
#'
#' Validation currently supports only the Study 02 single-trait prediction
#' task. Fields for unrelated benchmark types are neither required nor allowed
#' to influence this contract.
#' @param spec An ordinary benchmark specification list.
#' @return `spec`, invisibly.
#' @export
validate_benchmark_spec <- function(spec) {
  if (!is.list(spec)) stop("spec must be an ordinary R list.", call. = FALSE)
  required <- c("study", "task", "supported_profiles", "data", "split",
    "markers", "replicate_count", "scenarios", "methods", "controls",
    "seeds", "metrics", "validation", "frozen_capsule", "packages")
  missing <- setdiff(required, names(spec))
  if (length(missing))
    stop("spec is missing required fields: ", paste(missing, collapse = ", "),
      ".", call. = FALSE)
  .benchmark_scalar_string(spec$study, "spec$study")
  .benchmark_scalar_string(spec$task, "spec$task")
  if (!identical(spec$task, "single_trait_prediction"))
    stop("Unsupported benchmark task `", spec$task,
      "`; only `single_trait_prediction` is implemented.", call. = FALSE)
  if (!identical(spec$study, "02_prediction"))
    stop("The prediction execution contract currently supports only study `02_prediction`.",
      call. = FALSE)
  profiles <- spec$supported_profiles
  if (!is.list(profiles) || !identical(sort(names(profiles)),
      sort(c("workshop", "benchmark"))))
    stop("spec$supported_profiles must define exactly `workshop` and `benchmark`.",
      call. = FALSE)
  for (profile in names(profiles)) {
    value <- profiles[[profile]]
    if (!is.list(value) || !is.numeric(value$replicate_count) ||
        length(value$replicate_count) != 1L || is.na(value$replicate_count) ||
        value$replicate_count < 1)
      stop("Profile `", profile,
        "` must define one positive replicate_count.", call. = FALSE)
  }
  if (!identical(as.integer(profiles$benchmark$replicate_count),
      as.integer(spec$replicate_count)))
    stop("The benchmark profile replicate_count must equal spec$replicate_count.",
      call. = FALSE)
  if (!is.list(spec$scenarios) || !length(spec$scenarios) ||
      is.null(names(spec$scenarios)) || anyNA(names(spec$scenarios)) ||
      any(!nzchar(names(spec$scenarios))) || anyDuplicated(names(spec$scenarios)))
    stop("spec$scenarios must be a uniquely named non-empty list.", call. = FALSE)
  distributions <- vapply(spec$scenarios, function(x)
    if (is.list(x) && is.character(x$effect_distribution) &&
        length(x$effect_distribution) == 1L) x$effect_distribution else NA_character_,
    character(1))
  if (anyNA(distributions) ||
      any(!distributions %in% c("single_normal", "variance_mixture")))
    stop("Each scenario must define a supported effect_distribution.",
      call. = FALSE)
  if (!is.list(spec$methods) || length(spec$methods) != 4L ||
      is.null(names(spec$methods)) || anyDuplicated(names(spec$methods)))
    stop("spec$methods must be a uniquely named four-method list.", call. = FALSE)
  expected_methods <- c("st_bed_bayesc", "st_bed_bayesr",
    "st_csr_sbayesc", "st_csr_sbayesr")
  if (!identical(names(spec$methods), expected_methods))
    stop("spec$methods must preserve the exact Study 02 method order: ",
      paste(expected_methods, collapse = ", "), ".", call. = FALSE)
  if (!is.numeric(spec$split$train_fraction) ||
      length(spec$split$train_fraction) != 1L ||
      spec$split$train_fraction <= 0 || spec$split$train_fraction >= 1 ||
      length(spec$split$seed) != 1L || is.na(spec$split$seed))
    stop("spec$split must define a train_fraction in (0, 1) and one seed.",
      call. = FALSE)
  if (!is.numeric(spec$seeds$simulation_base) ||
      !is.numeric(spec$seeds$fit_offset) ||
      !is.numeric(spec$seeds$chain_stride))
    stop("spec$seeds must define numeric simulation_base, fit_offset, and chain_stride.",
      call. = FALSE)
  if (!is.character(spec$metrics) || !length(spec$metrics) || anyNA(spec$metrics))
    stop("spec$metrics must be a non-empty character vector.", call. = FALSE)
  .benchmark_scalar_string(spec$frozen_capsule, "spec$frozen_capsule")
  .benchmark_scalar_string(spec$packages$sblr$version,
    "spec$packages$sblr$version")
  sha <- .benchmark_scalar_string(spec$packages$sblr$sha,
    "spec$packages$sblr$sha")
  if (!grepl("^[0-9a-f]{40}$", sha))
    stop("spec$packages$sblr$sha must be a 40-character Git SHA.",
      call. = FALSE)
  invisible(spec)
}

resolve_benchmark_profile <- function(spec, profile = "benchmark") {
  validate_benchmark_spec(spec)
  .benchmark_scalar_string(profile, "profile")
  if (!profile %in% names(spec$supported_profiles))
    stop("Unknown profile `", profile, "`; choose one of: ",
      paste(names(spec$supported_profiles), collapse = ", "), ".",
      call. = FALSE)
  profile_spec <- spec$supported_profiles[[profile]]
  profile_spec$id <- profile
  profile_spec$replicate_count <- as.integer(profile_spec$replicate_count)
  profile_spec
}

#' Build benchmark coordinates
#'
#' @param spec A validated prediction benchmark specification.
#' @param profile One of the profiles supported by `spec`.
#' @return A deterministic scenario-by-replicate-by-method data frame.
#' @export
benchmark_coordinates <- function(spec, profile = "benchmark") {
  resolved <- resolve_benchmark_profile(spec, profile)
  out <- expand.grid(scenario = names(spec$scenarios),
    replicate = seq_len(resolved$replicate_count), method = names(spec$methods),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  out <- out[order(match(out$scenario, names(spec$scenarios)), out$replicate,
    match(out$method, names(spec$methods))), , drop = FALSE]
  rownames(out) <- NULL
  out$profile <- profile
  out
}

#' Add deterministic Study 02 seeds to benchmark coordinates
#'
#' @inheritParams benchmark_coordinates
#' @return The coordinate grid with simulation, fit, and chain seeds.
#' @export
benchmark_seeds <- function(spec, profile = "benchmark") {
  coordinates <- benchmark_coordinates(spec, profile)
  architecture_index <- match(coordinates$scenario, names(spec$scenarios))
  method_index <- match(coordinates$method, names(spec$methods))
  coordinates$architecture_seed <- as.integer(spec$seeds$simulation_base +
    architecture_index * 1000L)
  coordinates$simulation_seed <- as.integer(coordinates$architecture_seed +
    coordinates$replicate)
  coordinates$fit_seed <- as.integer(spec$seeds$fit_offset +
    architecture_index * 10000L + coordinates$replicate * 100L + method_index)
  nchains <- if (identical(profile, "benchmark"))
    as.integer(spec$controls$benchmark$nchains) else
    as.integer(spec$controls$workshop$nchains)
  coordinates$chain_seeds <- I(lapply(coordinates$fit_seed, function(seed)
    as.integer(seed + seq_len(nchains) * spec$seeds$chain_stride)))
  coordinates
}
