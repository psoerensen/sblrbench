# Study specifications are ordinary lists. Validation is deliberately scoped to
# fields required by the currently supported prediction, parameter-estimation,
# and convergence tasks.

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
  if (identical(spec$task, "convergence"))
    return(.validate_convergence_spec(spec))
  required <- c("study", "task", "supported_profiles", "data",
    "markers", "replicate_count", "scenarios", "methods", "controls",
    "seeds", "metrics", "validation", "frozen_capsule", "packages")
  if (identical(spec$task, "single_trait_prediction"))
    required <- c(required, "split")
  if (identical(spec$task, "parameter_estimation"))
    required <- c(required, "estimands")
  missing <- setdiff(required, names(spec))
  if (length(missing))
    stop("spec is missing required fields: ", paste(missing, collapse = ", "),
      ".", call. = FALSE)
  .benchmark_scalar_string(spec$study, "spec$study")
  .benchmark_scalar_string(spec$task, "spec$task")
  supported <- c("single_trait_prediction", "parameter_estimation")
  if (!spec$task %in% supported)
    stop("Unsupported benchmark task `", spec$task,
      "`; supported tasks are: ", paste(supported, collapse = ", "), ".",
      call. = FALSE)
  expected_study <- if (identical(spec$task, "single_trait_prediction"))
    "02_prediction" else "03_parameter_estimation"
  if (!identical(spec$study, expected_study))
    stop("Task `", spec$task, "` requires study `", expected_study, "`.",
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
    stop("spec$methods must preserve the exact four-method order: ",
      paste(expected_methods, collapse = ", "), ".", call. = FALSE)
  if (identical(spec$task, "single_trait_prediction") &&
      (!is.numeric(spec$split$train_fraction) ||
      length(spec$split$train_fraction) != 1L ||
      spec$split$train_fraction <= 0 || spec$split$train_fraction >= 1 ||
      length(spec$split$seed) != 1L || is.na(spec$split$seed)))
    stop("spec$split must define a train_fraction in (0, 1) and one seed.",
      call. = FALSE)
  if (!is.numeric(spec$seeds$simulation_base) ||
      !is.numeric(spec$seeds$fit_offset) ||
      !is.numeric(spec$seeds$chain_stride))
    stop("spec$seeds must define numeric simulation_base, fit_offset, and chain_stride.",
      call. = FALSE)
  if (!is.character(spec$metrics) || !length(spec$metrics) || anyNA(spec$metrics))
    stop("spec$metrics must be a non-empty character vector.", call. = FALSE)
  if (identical(spec$task, "parameter_estimation") &&
      (!is.data.frame(spec$estimands) || !nrow(spec$estimands) ||
       !all(c("estimand_id", "posterior_source", "truth_source") %in%
         names(spec$estimands)) || anyDuplicated(spec$estimands$estimand_id)))
    stop("Parameter-estimation spec$estimands must be a non-empty unique registry.",
      call. = FALSE)
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

#' Resolve a benchmark profile
#'
#' @param spec A validated benchmark specification.
#' @param profile A supported profile name.
#' @return The resolved profile list.
#' @export
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
  if (identical(spec$task, "convergence"))
    return(benchmark_convergence_coordinates(spec, profile))
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

#' Add deterministic benchmark seeds to coordinates
#'
#' @inheritParams benchmark_coordinates
#' @return The coordinate grid with simulation, fit, and chain seeds.
#' @export
benchmark_seeds <- function(spec, profile = "benchmark") {
  if (identical(spec$task, "convergence"))
    return(benchmark_convergence_seeds(spec, profile))
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

benchmark_matched_spec <- function(spec) {
  path <- benchmark_spec_path(spec, spec$matched_study$spec_path)
  matched <- read_benchmark_spec(path)
  if (!identical(matched$study, spec$matched_study$study) ||
      !identical(matched$task, "parameter_estimation"))
    stop("Study 04 matched_study must resolve to the Study 03 parameter-estimation specification.",
      call. = FALSE)
  matched
}

.validate_convergence_spec <- function(spec) {
  required <- c("study", "task", "supported_profiles", "replicate_count",
    "matched_study", "matched_grid", "controls", "seeds", "diagnostics",
    "validation", "frozen_capsules", "packages")
  missing <- setdiff(required, names(spec))
  if (length(missing))
    stop("Convergence spec is missing required fields: ",
      paste(missing, collapse = ", "), ".", call. = FALSE)
  if (!identical(spec$study, "04_convergence") ||
      !identical(spec$task, "convergence"))
    stop("The convergence task requires study `04_convergence`.", call. = FALSE)
  profiles <- spec$supported_profiles
  if (!is.list(profiles) || !identical(sort(names(profiles)),
      sort(c("workshop", "benchmark"))))
    stop("spec$supported_profiles must define exactly `workshop` and `benchmark`.",
      call. = FALSE)
  for (id in names(profiles)) {
    profile <- profiles[[id]]
    if (!is.numeric(profile$replicate_count) || length(profile$replicate_count) != 1L ||
        is.na(profile$replicate_count) || profile$replicate_count < 1L ||
        !is.character(profile$stages) || !length(profile$stages) ||
        any(!profile$stages %in% c("selection", "validation")))
      stop("Convergence profile `", id,
        "` must define a positive replicate_count and valid stages.", call. = FALSE)
  }
  if (!identical(as.integer(profiles$benchmark$replicate_count),
      as.integer(spec$replicate_count)))
    stop("The convergence benchmark replicate_count must equal spec$replicate_count.",
      call. = FALSE)
  required_grid <- c("scenario", "method")
  if (!is.data.frame(spec$matched_grid) || nrow(spec$matched_grid) != 4L ||
      !all(required_grid %in% names(spec$matched_grid)) ||
      anyDuplicated(spec$matched_grid[required_grid]))
    stop("spec$matched_grid must contain four unique scenario-method rows.",
      call. = FALSE)
  diagnostic_fields <- c("registry", "burnin_candidates",
    "retained_draw_candidates", "thresholds", "recommendation_rules")
  absent <- setdiff(diagnostic_fields, names(spec$diagnostics))
  if (length(absent))
    stop("spec$diagnostics is missing: ", paste(absent, collapse = ", "),
      ".", call. = FALSE)
  registry <- spec$diagnostics$registry
  if (!is.data.frame(registry) || !nrow(registry) ||
      !all(c("quantity", "source", "classification", "required") %in%
        names(registry)) || anyDuplicated(registry$quantity))
    stop("Convergence diagnostic registry must be a unique non-empty data frame.",
      call. = FALSE)
  thresholds <- spec$diagnostics$thresholds
  needed_thresholds <- c("rhat", "ess_bulk", "ess_tail", "relative_mcse",
    "chain_count", "standardized_mean_shift")
  if (!is.list(thresholds) || !all(needed_thresholds %in% names(thresholds)) ||
      any(!is.finite(unlist(thresholds[needed_thresholds]))))
    stop("Convergence thresholds are incomplete or non-finite.", call. = FALSE)
  for (id in c("burnin_candidates", "retained_draw_candidates")) {
    value <- spec$diagnostics[[id]]
    if (!is.numeric(value) || !length(value) || anyNA(value) ||
        any(value < 0) || is.unsorted(value, strictly = TRUE))
      stop("spec$diagnostics$", id,
        " must be a strictly increasing non-negative numeric vector.",
        call. = FALSE)
  }
  if (!is.list(spec$controls$selection) ||
      !all(c("nit", "nburn", "nthin", "nchains", "ncores") %in%
        names(spec$controls$selection)))
    stop("spec$controls$selection is incomplete.", call. = FALSE)
  if (!is.numeric(spec$seeds$fit_base) ||
      !is.numeric(spec$seeds$chain_stride))
    stop("Convergence seeds must define fit_base and chain_stride.",
      call. = FALSE)
  if (!is.list(spec$frozen_capsules) ||
      !all(c("selection", "validation") %in% names(spec$frozen_capsules)))
    stop("Convergence spec must identify selection and validation capsules.",
      call. = FALSE)
  sha <- .benchmark_scalar_string(spec$packages$sblr$sha,
    "spec$packages$sblr$sha")
  if (!grepl("^[0-9a-f]{40}$", sha))
    stop("spec$packages$sblr$sha must be a 40-character Git SHA.",
      call. = FALSE)
  invisible(spec)
}
