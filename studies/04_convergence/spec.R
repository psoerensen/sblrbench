# Exact scientific specification for Study 04. The benchmark profile preserves
# both historical stages: four maximum-history selection fits followed by
# twenty fixed-setting validation fits. Workshop exposes the four-coordinate
# selection design only and is unsuitable for performance claims.

spec <- list(
  study = "04_convergence",
  task = "convergence",
  supported_profiles = list(
    workshop = list(
      replicate_count = 1L,
      stages = "selection",
      suitable_for_method_performance_claims = FALSE),
    benchmark = list(
      replicate_count = 5L,
      stages = c("selection", "validation"),
      suitable_for_method_performance_claims = FALSE)),
  replicate_count = 5L,
  matched_study = list(
    study = "03_parameter_estimation",
    spec_path = file.path("studies", "03_parameter_estimation", "spec.R")),
  matched_grid = data.frame(
    scenario = c("sparse_homogeneous", "sparse_homogeneous",
      "sparse_mixture", "sparse_mixture"),
    method = c("st_bed_bayesc", "st_csr_sbayesc", "st_bed_bayesr",
      "st_csr_sbayesr"), stringsAsFactors = FALSE),
  controls = list(
    selection = list(nit = 3000L, nburn = 0L, nthin = 1L,
      nchains = 4L, ncores = 4L, convergence = "core",
      keep_chains = TRUE,
      convergence_control = list(warn = FALSE, keep_traces = TRUE)),
    validation = list(nthin = 1L, nchains = 4L, ncores = 4L,
      convergence = "core", keep_chains = TRUE,
      convergence_control = list(warn = FALSE, keep_traces = TRUE),
      recommendation_source = file.path("results", "reference",
        "04_convergence", "current-selection",
        "method_recommendations.csv"))),
  seeds = list(fit_base = 40000L, chain_stride = 100L),
  diagnostics = list(
    registry = data.frame(
      quantity = c("effect_variance", "genetic_variance",
        "residual_variance", "heritability", "causal_proportion",
        "total_marker_effect_variance"),
      source = c("vbs", "vgs", "ves", "vgs/(vgs+ves)", "pi_trace",
        "vbs*pi_trace*marker_count"),
      classification = c(rep("core", 4L), "unavailable", "unavailable"),
      required = c(rep(TRUE, 4L), FALSE, FALSE),
      lower_bound = 0,
      upper_bound = c(Inf, Inf, Inf, 1, 1, Inf),
      stringsAsFactors = FALSE),
    burnin_candidates = c(250L, 500L, 1000L),
    retained_draw_candidates = c(250L, 500L, 1000L, 2000L),
    thresholds = list(rhat = 1.01, ess_bulk = 400, ess_tail = 400,
      relative_mcse = 0.05, chain_count = 4L,
      standardized_mean_shift = 0.10),
    recommendation_rules = list(
      reference_burnin = 1000L,
      reference_retained_draws = 2000L,
      selection = "earliest stable burn-in and retained-draw candidate for which all current and longer windows pass every core threshold",
      unavailable = "no stable checkpoint passed all core thresholds")),
  validation = list(
    oracle_tolerance = 1e-10,
    required_successful_replicates = 5L,
    status_rule = "indeterminate if fewer than five successful replicates; otherwise all=5, most=3-4, mixed=1-2, not_supported=0",
    expected_selection_fit_count = 4L,
    expected_validation_fit_count = 20L,
    expected_selection_chain_count = 16L,
    expected_validation_chain_count = 80L),
  frozen_capsules = list(
    selection = file.path("results", "reference", "04_convergence",
      "current-selection"),
    validation = file.path("results", "reference", "04_convergence",
      "current-validation")),
  packages = list(
    sblr = list(version = "0.2.0",
      sha = "02e8c74baa906e83c4a08d42a9cc6339b4e81072"),
    qgg = list(version = "1.1.6"),
    sblrbench = list(version = "0.0.0.9000"),
    qgdata = list(sha = "6cca5819e711d326cfb2614d7e9d9f34942612cd")))
