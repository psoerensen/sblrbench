# Study specifications are ordinary lists. Validation is deliberately scoped to
# fields required by the currently supported prediction, parameter-estimation,
# convergence, fine-mapping, LD-operator, and annotation-model tasks.

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
  if (identical(spec$task, "finemapping"))
    return(.validate_finemapping_spec(spec))
  if (identical(spec$task, "ld_operator"))
    return(.validate_ld_operator_spec(spec))
  if (identical(spec$task, "annotation_models"))
    return(.validate_annotation_spec(spec))
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
  if (identical(spec$task, "ld_operator"))
    return(benchmark_ld_operator_coordinates(spec, profile))
  if (identical(spec$task, "annotation_models"))
    return(benchmark_annotation_coordinates(spec, profile, mode = "final"))
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
  if (identical(spec$task, "ld_operator")) {
    coordinates <- benchmark_ld_operator_coordinates(spec, profile)
    architecture_index <- match(coordinates$scenario, names(spec$scenarios))
    configuration_index <- match(coordinates$configuration,
      spec$operators$configurations)
    coordinates$simulation_seed <- as.integer(spec$seeds$simulation_base +
      (architecture_index - 1L) * spec$seeds$architecture_stride +
      coordinates$replicate * spec$seeds$replicate_stride)
    coordinates$effect_seed <- coordinates$simulation_seed +
      as.integer(spec$seeds$effect_offset)
    coordinates$residual_seed <- coordinates$simulation_seed +
      as.integer(spec$seeds$residual_offset)
    coordinates$fit_seed <- as.integer(coordinates$simulation_seed +
      spec$seeds$fit_base + configuration_index *
        spec$seeds$configuration_stride)
    coordinates$chain_seeds <- I(lapply(coordinates$fit_seed, function(seed)
      as.integer(seed + seq_len(4L) * spec$seeds$chain_stride)))
    return(coordinates)
  }
  if (identical(spec$task, "annotation_models"))
    return(benchmark_annotation_seeds(spec, profile, mode = "final"))
  coordinates <- benchmark_coordinates(spec, profile)
  if (identical(spec$task, "finemapping")) {
    method_index <- match(coordinates$method, names(spec$methods))
    coordinates$causal_seed <- as.integer(spec$seeds$simulation_base +
      coordinates$replicate)
    coordinates$simulation_seed <- as.integer(spec$seeds$simulation_base +
      spec$seeds$phenotype_offset + coordinates$replicate)
    coordinates$fit_seed <- as.integer(spec$seeds$fit_offset +
      coordinates$replicate * 100L + method_index)
    coordinates$chain_seeds <- I(lapply(coordinates$fit_seed, as.integer))
    return(coordinates)
  }
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

#' Build Study 05 LD-operator coordinates
#'
#' @inheritParams benchmark_coordinates
#' @return A deterministic architecture-by-replicate-by-configuration table.
#' @export
benchmark_ld_operator_coordinates <- function(spec, profile = "benchmark") {
  resolved <- resolve_benchmark_profile(spec, profile)
  out <- expand.grid(
    scenario = names(spec$scenarios),
    replicate = seq_len(resolved$replicate_count),
    configuration = spec$operators$configurations,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  out <- out[order(match(out$scenario, names(spec$scenarios)), out$replicate,
    match(out$configuration, spec$operators$configurations)), , drop = FALSE]
  rownames(out) <- NULL
  out$method <- paste(out$scenario, out$configuration, sep = "__")
  out$operator <- out$configuration
  out$profile <- profile
  out
}

.validate_ld_operator_spec <- function(spec) {
  required <- c("study", "task", "supported_profiles", "data", "split",
    "markers", "replicate_count", "scenarios", "methods", "controls",
    "seeds", "operators", "metrics", "validation", "frozen_capsule",
    "components", "packages")
  missing <- setdiff(required, names(spec))
  if (length(missing))
    stop("LD-operator spec is missing required fields: ",
      paste(missing, collapse = ", "), ".", call. = FALSE)
  if (!identical(spec$study, "05_ld_operator") ||
      !identical(spec$task, "ld_operator"))
    stop("The LD-operator task requires study `05_ld_operator`.", call. = FALSE)
  profiles <- spec$supported_profiles
  if (!is.list(profiles) || !identical(sort(names(profiles)),
      c("benchmark", "workshop")))
    stop("LD-operator profiles must define workshop and benchmark.", call. = FALSE)
  for (id in names(profiles)) {
    count <- profiles[[id]]$replicate_count
    if (!is.numeric(count) || length(count) != 1L || is.na(count) || count < 1L)
      stop("LD-operator profile `", id,
        "` must define one positive replicate_count.", call. = FALSE)
  }
  if (!identical(as.integer(profiles$benchmark$replicate_count),
      as.integer(spec$replicate_count)))
    stop("LD-operator benchmark replicate_count must equal spec$replicate_count.",
      call. = FALSE)
  configurations <- c("bed", "full_csr", "block_csr", "low_rank_full",
    "low_rank_0999", "low_rank_0995")
  if (!identical(spec$operators$configurations, configurations) ||
      !identical(names(spec$methods), configurations))
    stop("LD-operator configurations must preserve the historical six-operator order.",
      call. = FALSE)
  if (!identical(names(spec$scenarios),
      c("sparse_homogeneous", "sparse_mixture")))
    stop("LD-operator scenarios must preserve the two historical architectures.",
      call. = FALSE)
  recommendations <- spec$controls$benchmark$recommendations
  control_columns <- c("scenario", "configuration", "nit", "nburn",
    "nthin", "nchains", "ncores")
  if (!is.data.frame(recommendations) || nrow(recommendations) != 12L ||
      !all(control_columns %in% names(recommendations)) ||
      anyDuplicated(recommendations[c("scenario", "configuration")]) ||
      !setequal(recommendations$scenario, names(spec$scenarios)) ||
      !setequal(recommendations$configuration, configurations) ||
      any(recommendations$nchains != 4L) || any(recommendations$ncores != 4L) ||
      any(recommendations$nthin != 1L))
    stop("LD-operator benchmark MCMC recommendations are incomplete or changed.",
      call. = FALSE)
  block <- spec$operators$block
  if (!is.numeric(block$size) || length(block$size) != 1L || block$size != 1000L ||
      !identical(as.integer(block$sensitivity_sizes),
        c(250L, 500L, 1000L, 2000L)))
    stop("LD-operator block definitions do not match the frozen design.",
      call. = FALSE)
  eigen <- spec$operators$eigen
  expected_props <- c(low_rank_full = 1 - .Machine$double.eps,
    low_rank_0999 = 0.999, low_rank_0995 = 0.995)
  if (!identical(eigen$policy, "cumulative_positive_mass") ||
      !isTRUE(all.equal(eigen$proportions, expected_props)) ||
      !identical(eigen$tolerance, 1e-10))
    stop("LD-operator eigen-retention policy differs from the frozen design.",
      call. = FALSE)
  tolerances <- spec$operators$equivalence_tolerances
  needed <- c("absolute", "relative", "product_absolute",
    "quadratic_absolute", "probability")
  if (!is.list(tolerances) || !all(needed %in% names(tolerances)) ||
      any(!is.finite(unlist(tolerances[needed]))))
    stop("LD-operator equivalence tolerances are incomplete.", call. = FALSE)
  if (!identical(as.integer(spec$validation$expected_fit_count), 60L))
    stop("LD-operator validation must retain the 60-fit benchmark contract.",
      call. = FALSE)
  .benchmark_scalar_string(spec$frozen_capsule, "spec$frozen_capsule")
  if (!identical(sort(names(spec$components)),
      c("operator_validation", "sbayesr_ld_sensitivity")))
    stop("LD-operator spec must define the integrated operator-validation and ",
      "SBayesR LD-sensitivity components.", call. = FALSE)
  supplemental <- spec$components$sbayesr_ld_sensitivity
  if (!identical(as.integer(supplemental$marker_window$count), 1500L) ||
      !identical(as.integer(supplemental$block_starts),
        c(1L, 251L, 501L, 624L, 751L, 1001L, 1251L)) ||
      !identical(as.integer(supplemental$retained_block_ranks),
        c(248L, 248L, 123L, 127L, 248L, 248L, 248L)) ||
      !identical(as.integer(supplemental$retained_total_rank), 1490L) ||
      !identical(supplemental$retained_mass, 0.995))
    stop("Supplemental Study 05 window, blocks, or retained-rank policy changed.",
      call. = FALSE)
  sha <- .benchmark_scalar_string(spec$packages$sblr$sha,
    "spec$packages$sblr$sha")
  if (!grepl("^[0-9a-f]{40}$", sha))
    stop("spec$packages$sblr$sha must be a 40-character Git SHA.",
      call. = FALSE)
  invisible(spec)
}

.validate_finemapping_spec <- function(spec) {
  required <- c("study", "task", "supported_profiles", "data", "markers",
    "replicate_count", "scenarios", "causal_design", "locus_design",
    "methods", "controls", "seeds", "estimands", "metrics", "validation",
    "frozen_capsule", "packages")
  missing <- setdiff(required, names(spec))
  if (length(missing)) stop("Fine-mapping spec is missing required fields: ",
    paste(missing, collapse = ", "), ".", call. = FALSE)
  if (!identical(spec$study, "01_finemapping") ||
      !identical(spec$task, "finemapping"))
    stop("The fine-mapping task requires study `01_finemapping`.", call. = FALSE)
  profiles <- spec$supported_profiles
  if (!is.list(profiles) || !identical(sort(names(profiles)),
      c("benchmark", "workshop")))
    stop("Fine-mapping profiles must define workshop and benchmark.", call. = FALSE)
  if (!identical(as.integer(profiles$benchmark$replicate_count),
      as.integer(spec$replicate_count)))
    stop("Fine-mapping benchmark replicate_count must equal spec$replicate_count.",
      call. = FALSE)
  expected <- c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc",
    "st_csr_sbayesr")
  if (!identical(names(spec$methods), expected))
    stop("Fine-mapping methods must preserve the historical four-method order.",
      call. = FALSE)
  causal <- spec$causal_design
  if (!is.numeric(causal$min_distance_bp) || causal$min_distance_bp < 0 ||
      !is.numeric(causal$min_maf) || !is.numeric(causal$max_maf))
    stop("Fine-mapping causal-marker design is incomplete.", call. = FALSE)
  locus <- spec$locus_design
  needed <- c("credible_set_target", "min_r2", "pip_cutoff",
    "locus_pip_cutoff", "max_locus_distance", "algorithm")
  if (!all(needed %in% names(locus)) ||
      !identical(locus$algorithm, "sblr::make_credible_sets") ||
      locus$credible_set_target <= 0 || locus$credible_set_target > 1)
    stop("Fine-mapping locus and credible-set design is invalid.", call. = FALSE)
  for (profile in names(profiles)) {
    controls <- spec$controls[[profile]]
    if (!all(c("nit", "nburn", "nthin", "nchains", "ncores") %in%
        names(controls)) || controls$nchains != 1L)
      stop("Fine-mapping profile controls must preserve the one-chain policy.",
        call. = FALSE)
  }
  if (!is.data.frame(spec$estimands) || !nrow(spec$estimands) ||
      !all(c("estimand", "definition") %in% names(spec$estimands)))
    stop("Fine-mapping estimands must be a non-empty data frame.", call. = FALSE)
  sha <- .benchmark_scalar_string(spec$packages$sblr$sha,
    "spec$packages$sblr$sha")
  if (!grepl("^[0-9a-f]{40}$", sha))
    stop("spec$packages$sblr$sha must be a 40-character Git SHA.", call. = FALSE)
  invisible(spec)
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

#' Build Study 06 annotation-model coordinates
#'
#' @param spec A validated Study 06 specification.
#' @param profile A supported profile.
#' @param mode Either `"qualification"` or `"final"`.
#' @return A deterministic coordinate data frame.
#' @export
benchmark_annotation_coordinates <- function(spec, profile = "benchmark",
                                             mode = c("final",
                                               "qualification")) {
  mode <- match.arg(mode)
  validate_benchmark_spec(spec)
  resolved <- resolve_benchmark_profile(spec, profile)
  if (identical(mode, "qualification")) {
    out <- spec$qualification$entries
  } else {
    out <- expand.grid(scenario = names(spec$scenarios),
      replicate = seq_len(resolved$replicate_count),
      method = names(spec$methods), KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE)
    out <- out[order(match(out$scenario, names(spec$scenarios)),
      out$replicate, match(out$method, names(spec$methods))), , drop = FALSE]
  }
  rownames(out) <- NULL
  out$mode <- mode
  out$profile <- profile
  out
}

#' Add Study 06 deterministic seeds
#'
#' @inheritParams benchmark_annotation_coordinates
#' @return Study 06 coordinates with simulation, fit, and chain seeds.
#' @export
benchmark_annotation_seeds <- function(spec, profile = "benchmark",
                                       mode = c("final", "qualification")) {
  mode <- match.arg(mode)
  coordinates <- benchmark_annotation_coordinates(spec, profile, mode)
  scenario_index <- match(coordinates$scenario, names(spec$scenarios))
  method_index <- match(coordinates$method, names(spec$methods))
  base <- as.integer(spec$seeds$scenario_base +
    scenario_index * spec$seeds$scenario_stride +
    coordinates$replicate * spec$seeds$replicate_stride)
  coordinates$component_seed <- base + spec$seeds$component_offset
  coordinates$effect_seed <- base + spec$seeds$effect_offset
  coordinates$residual_seed <- base + spec$seeds$residual_offset
  coordinates$fit_seed <- as.integer(spec$seeds$fit_base +
    scenario_index * spec$seeds$scenario_stride +
    coordinates$replicate * 1000L +
    method_index * spec$seeds$method_stride)
  chain_count <- if (identical(mode, "qualification") ||
      identical(profile, "benchmark")) 4L else
    as.integer(spec$controls$workshop$nchains)
  coordinates$chain_seeds <- I(lapply(coordinates$fit_seed, function(seed)
    as.integer(seed + seq_len(chain_count) * spec$seeds$chain_stride)))
  coordinates
}

.validate_annotation_spec <- function(spec) {
  required <- c("study", "task", "status", "supported_profiles", "data",
    "split", "markers", "replicate_count", "scenarios",
    "annotation_design", "methods", "controls", "qualification", "seeds",
    "estimands", "metrics", "validation", "frozen_capsule", "packages")
  missing <- setdiff(required, names(spec))
  if (length(missing))
    stop("Annotation-model spec is missing required fields: ",
      paste(missing, collapse = ", "), ".", call. = FALSE)
  if (!identical(spec$study, "06_annotation_models") ||
      !identical(spec$task, "annotation_models"))
    stop("The annotation-model task requires study `06_annotation_models`.",
      call. = FALSE)
  if (!identical(spec$status, "qualification_pending"))
    stop("Study 06 must remain `qualification_pending` until fitting passes.",
      call. = FALSE)
  profiles <- spec$supported_profiles
  if (!is.list(profiles) || !identical(sort(names(profiles)),
      c("benchmark", "workshop")))
    stop("Annotation-model profiles must define workshop and benchmark.",
      call. = FALSE)
  if (!identical(as.integer(profiles$benchmark$replicate_count), 5L) ||
      !identical(as.integer(spec$replicate_count), 5L))
    stop("Study 06 benchmark must retain five replicates.", call. = FALSE)
  if (!identical(as.integer(spec$data$expected_sample_count), 2000L) ||
      !identical(as.integer(spec$data$expected_marker_count), 37991L) ||
      !identical(as.integer(spec$data$chromosome), 1L))
    stop("Study 06 data counts or chromosome differ from the audited design.",
      call. = FALSE)
  if (!identical(spec$data$sparse_ld, list(max_distance_bp = 0,
      max_distance_variants = 1000L, r2_threshold = .001,
      block_size = 1024L, nthreads = 1L)))
    stop("Study 06 sparse-LD controls changed.", call. = FALSE)
  if (!identical(spec$data$example_data$commit,
      "6cca5819e711d326cfb2614d7e9d9f34942612cd") ||
      !identical(spec$packages$qgdata$sha,
        "6cca5819e711d326cfb2614d7e9d9f34942612cd"))
    stop("Study 06 must use the validated qgdata SHA.", call. = FALSE)
  if (!identical(spec$split[c("train_fraction", "seed")],
      list(train_fraction = 0.70, seed = 3101L)))
    stop("Study 06 train/test split differs from the audited design.",
      call. = FALSE)
  expected_scenarios <- c("informative_annotations",
    "uninformative_annotations")
  if (!identical(names(spec$scenarios), expected_scenarios))
    stop("Study 06 scenarios must preserve informative then uninformative order.",
      call. = FALSE)
  expected_columns <- c("Intercept", "enriched_binary",
    "continuous_signal", "null_annotation")
  if (!identical(spec$annotation_design$columns, expected_columns))
    stop("Study 06 annotation columns must preserve this exact order: ",
      paste(expected_columns, collapse = ", "), ".", call. = FALSE)
  alpha <- spec$annotation_design$informative_nonintercept_alpha
  if (!is.matrix(alpha) || !identical(dim(alpha), c(3L, 3L)) ||
      !identical(rownames(alpha), expected_columns[-1L]) ||
      !identical(colnames(alpha), paste0("step_", 1:3)) ||
      !isTRUE(all.equal(unname(alpha), matrix(c(.98, .20, 0, .25, 0,
        .20, 0, 0, 0), 3L, byrow = TRUE))))
    stop("Study 06 informative annotation coefficients changed.",
      call. = FALSE)
  expected_methods <- c("st_bed_bayesr", "st_bed_bayesrc",
    "st_csr_sbayesr", "st_csr_sbayesrc")
  if (!identical(names(spec$methods), expected_methods))
    stop("Study 06 methods must preserve the audited four-method order.",
      call. = FALSE)
  observed_methods <- vapply(spec$methods, function(x) paste(x$interface,
    x$native_method, x$representation, x$prior_class,
    if (is.null(x$annotation_model)) "" else x$annotation_model,
    x$computational_policy, sep = "|"), character(1))
  expected_method_contracts <- c(
    st_bed_bayesr = "stblr_bed|bayesr|BED|BayesR||full unscheduled BED sweep",
    st_bed_bayesrc = "stblr_bed|bayesrc|BED|BayesRC||full unscheduled BED sweep",
    st_csr_sbayesr = "stblr_csr|sbayesr|CSR|BayesR||training-only sparse LD",
    st_csr_sbayesrc = paste("stblr_csr_annot|sbayesrc|CSR|BayesRC",
      "annotation_probit_stick|training-only sparse LD", sep = "|"))
  if (!identical(observed_methods, expected_method_contracts))
    stop("Study 06 method interfaces, models, or computational policies changed.",
      call. = FALSE)
  annotation_methods <- names(spec$methods)[vapply(spec$methods,
    function(x) isTRUE(x$annotation_aware), logical(1))]
  if (!identical(annotation_methods,
      c("st_bed_bayesrc", "st_csr_sbayesrc")))
    stop("Study 06 baseline/annotation-aware method pairing is invalid.",
      call. = FALSE)
  simulation <- spec$controls$simulation
  if (!identical(simulation$h2, .30) ||
      !identical(as.integer(simulation$target_expected_nonnull), 50L) ||
      !identical(simulation$mixture_var, c(0, .01, .1, 1)) ||
      !identical(simulation$active_component_weights, c(.60, .30, .10)))
    stop("Study 06 simulation or mixture controls changed.", call. = FALSE)
  priors <- spec$controls$priors
  if (!identical(priors$h2, .30) ||
      !identical(priors$bayesr_active_probability, 50 / 37991) ||
      !identical(priors$bayesr_mixture_var, c(0, .01, .1, 1)) ||
      !identical(priors$intercept_flat, TRUE) ||
      !identical(priors$sigmaSqAlpha_init, c(1, 1, 1)) ||
      !identical(priors[c("sigmaSqAlpha_a", "sigmaSqAlpha_b")],
        list(sigmaSqAlpha_a = 2, sigmaSqAlpha_b = 2)) ||
      !identical(priors$alpha_update_every, 1L) ||
      !identical(priors$pi_floor, 1e-12) ||
      !all(vapply(priors[c("updateAlpha", "updateB", "updateE")],
        isTRUE, logical(1))))
    stop("Study 06 priors or annotation update controls changed.",
      call. = FALSE)
  entries <- spec$qualification$entries
  expected_entries <- data.frame(
    scenario = rep(expected_scenarios, each = 2L), replicate = 1L,
    method = rep(c("st_bed_bayesrc", "st_csr_sbayesrc"), 2L),
    stringsAsFactors = FALSE)
  if (!is.data.frame(entries) || !identical(entries, expected_entries))
    stop("Study 06 must define exactly the four audited qualification entries.",
      call. = FALSE)
  qualification <- spec$qualification
  if (!identical(as.integer(qualification$maximum_history), 9000L) ||
      !identical(as.integer(qualification$burnin_candidates),
        c(1000L, 2000L, 3000L)) ||
      !identical(as.integer(qualification$retained_candidates),
        c(2000L, 4000L, 6000L)) ||
      !identical(as.integer(qualification$nchains), 4L))
    stop("Study 06 qualification history or candidate windows changed.",
      call. = FALSE)
  thresholds <- qualification$thresholds
  if (!identical(thresholds[c("rhat", "ess_bulk", "ess_tail",
      "relative_mcse", "chain_count")], list(rhat = 1.01,
      ess_bulk = 400, ess_tail = 400, relative_mcse = .05,
      chain_count = 4L)))
    stop("Study 06 convergence thresholds changed.", call. = FALSE)
  candidates <- expand.grid(burnin = qualification$burnin_candidates,
    retained = qualification$retained_candidates)
  if (!any(candidates$burnin + candidates$retained <=
      qualification$maximum_history))
    stop("Study 06 qualification has no feasible candidate window.",
      call. = FALSE)
  if (!is.data.frame(spec$estimands) || !nrow(spec$estimands) ||
      anyDuplicated(spec$estimands$estimand_id) ||
      !all(c("estimand_id", "definition", "required") %in%
        names(spec$estimands)))
    stop("Study 06 estimand registry is invalid.", call. = FALSE)
  if (!is.character(spec$metrics) || !length(spec$metrics) ||
      anyNA(spec$metrics) || anyDuplicated(spec$metrics))
    stop("Study 06 metrics must be a unique non-empty character vector.",
      call. = FALSE)
  if (!identical(spec$frozen_capsule$status, "partial_evidence_only") ||
      !is.null(spec$frozen_capsule$final))
    stop("Study 06 must not claim a final capsule before qualification.",
      call. = FALSE)
  sha <- .benchmark_scalar_string(spec$packages$sblr$sha,
    "spec$packages$sblr$sha")
  if (!identical(sha, "02e8c74baa906e83c4a08d42a9cc6339b4e81072"))
    stop("Study 06 must use the validated sblr SHA.", call. = FALSE)
  # Validate the complete seed space without recursively calling public helpers.
  grid <- expand.grid(scenario = expected_scenarios, replicate = 1:5,
    method = expected_methods, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  s <- match(grid$scenario, expected_scenarios)
  m <- match(grid$method, expected_methods)
  fit <- spec$seeds$fit_base + s * spec$seeds$scenario_stride +
    grid$replicate * 1000L + m * spec$seeds$method_stride
  chains <- unlist(lapply(fit, function(seed)
    seed + seq_len(4L) * spec$seeds$chain_stride), use.names = FALSE)
  if (anyDuplicated(fit) || anyDuplicated(chains))
    stop("Study 06 fit or chain seeds are not unique.", call. = FALSE)
  invisible(spec)
}
