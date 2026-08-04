# Study 05: integrated LD-operator validation.
# Scientific values are transcribed from the completed main and SBayesR capsules.

.study05_configurations <- c(
  "bed", "full_csr", "block_csr", "low_rank_full",
  "low_rank_0999", "low_rank_0995"
)

.study05_architectures <- list(
  sparse_homogeneous = list(
    target_h2 = 0.30, n_causal = 50L,
    effect_distribution = "single_normal",
    mixture_prob = NA, mixture_var = NA,
    mixture_probabilities = NA, mixture_multipliers = NA
  ),
  sparse_mixture = list(
    target_h2 = 0.30, n_causal = 50L,
    effect_distribution = "variance_mixture",
    mixture_prob = c(0.60, 0.30, 0.10),
    mixture_var = c(0.01, 0.1, 1),
    mixture_probabilities = c(0.60, 0.30, 0.10),
    mixture_multipliers = c(0.01, 0.1, 1)
  )
)

.study05_methods <- lapply(.study05_configurations, function(configuration) {
  list(
    id = configuration,
    configuration = configuration,
    operator = configuration,
    interface = if (configuration == "bed") "stblr_bed" else if (
      startsWith(configuration, "low_rank_")) "stblr_block_eigen" else "stblr_csr",
    representation = if (startsWith(configuration, "low_rank_"))
      "low_rank" else NA_character_
  )
})
names(.study05_methods) <- .study05_configurations

spec <- list(
  study = "05_ld_operator",
  title = "Study 05 — LD-operator validation",
  task = "ld_operator",
  supported_profiles = list(
    workshop = list(replicate_count = 1L),
    benchmark = list(replicate_count = 5L)
  ),
  data = list(
    source = "qgg_example",
    chromosome = 1L,
    trait = "trait1",
    sample_limit = NULL,
    sample_count = 5000L,
    marker_count = 37991L,
    sample_policy = "all samples; fixed 70/30 train/test split",
    marker_policy = "all chromosome-1 HapMap3 markers in source order",
    alignment = "strict sample and marker identifiers; no silent reordering",
    summary_statistics = "training samples only; scaling inherited from qgdata Glist",
    qgdata_sha = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
    example_data = list(
      repository = "psoerensen/qgdata",
      commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
      subdirectory = "simulated_human_data",
      files = c("human.bed", "human.bim", "human.fam", "human.pheno",
        "human.covar"),
      size_bytes = c(human.bed = 62500003, human.bim = 1882359,
        human.fam = 117786, human.pheno = 92786, human.covar = 641513),
      md5 = c(human.bed = "e89bea9a6cedd9eeef3fd0a5c807db81",
        human.bim = "0105119b04c67b7ac7f66cc5e6680963",
        human.fam = "3c5db3d9eb7f3fc893c75f6f2b89836d",
        human.pheno = "6a9e7cb1162e43999c170a363863176d",
        human.covar = "d06002aa2b1b79bdc4c0e92c21f27ae5")
    ),
    sparse_ld = list(max_distance_bp = 0, max_distance_variants = 1000L,
      r2_threshold = 0.001, block_size = 1024L, nthreads = 1L)
  ),
  split = list(train_fraction = 0.70, seed = 3101L),
  markers = list(selection = "all_qgdata_chr1_hapmap3", count = 37991L,
    qc = list(excludeMAF = 0.05, excludeMISS = 0.05, excludeCGAT = TRUE,
      excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
      excludeMHC = FALSE)),
  replicate_count = 5L,
  scenarios = .study05_architectures,
  methods = .study05_methods,
  controls = list(
    simulation = list(n_causal = 50L, h2 = 0.30),
    workshop = list(nit = 50L, nburn = 25L, nthin = 1L,
      nchains = 4L, ncores = 4L),
    benchmark = list(nchains = 4L, ncores = 4L,
      recommendations = data.frame(
        scenario = rep(c("sparse_homogeneous", "sparse_mixture"), each = 6L),
        configuration = rep(.study05_configurations, 2L),
        nit = c(rep(250L, 6L), 2000L, 2000L, 1000L, 1000L, 1000L, 1000L),
        nburn = c(250L, 250L, 250L, 250L, 500L, 250L,
          250L, 250L, 250L, 250L, 250L, 500L),
        nthin = 1L, nchains = 4L, ncores = 4L,
        stringsAsFactors = FALSE
      )),
    convergence = list(
      maximum_nit = 3000L, nburn = 0L, nthin = 1L,
      candidate_burnin = c(250L, 500L, 1000L),
      candidate_retained = c(250L, 500L, 1000L, 2000L),
      thresholds = list(rhat = 1.01, ess_bulk = 400, ess_tail = 400,
        relative_mcse = 0.05, chains = 4L)
    ),
    residual_rebuild_every = 100L,
    bayesc = list(pi_init = 0.01),
    bayesr = list(
      pi = c(0.99, 0.01 / 3, 0.01 / 3, 0.01 / 3),
      mixture_var = c(0, 0.01, 0.1, 1)
    )
  ),
  seeds = list(
    data_selection = 3101L,
    simulation_base = 61000L,
    architecture_stride = 10000L,
    replicate_stride = 100L,
    effect_offset = 11L,
    residual_offset = 37L,
    fit_base = 600000L,
    configuration_stride = 1000L,
    chain_stride = 101L,
    operator_probe = 6601L
  ),
  operators = list(
    contract = "block_low_rank_v1",
    configurations = .study05_configurations,
    block = list(size = 1000L, sensitivity_sizes = c(250L, 500L, 1000L, 2000L),
      policy = "contiguous marker order"),
    sparse = list(threshold = 0, window = 0L, source = "full CSR"),
    eigen = list(
      policy = "cumulative_positive_mass",
      tolerance = 1e-10,
      proportions = c(low_rank_full = 1 - .Machine$double.eps,
        low_rank_0999 = 0.999, low_rank_0995 = 0.995),
      full_rank_rule = "retain every numerically positive eigenmode"
    ),
    equivalence_tolerances = list(
      absolute = 2e-3, relative = 2e-6,
      product_absolute = 2e-2, quadratic_absolute = 5e-2,
      probability = 1e-10
    )
  ),
  metrics = c(
    "prediction_correlation", "genetic_value_correlation", "prediction_nmse",
    "genetic_value_rmse", "heritability_recovery", "marker_effect_agreement",
    "variable_selection", "operator_error", "runtime"
  ),
  validation = list(
    expected_fit_count = 60L,
    successful_fit_count = 60L,
    benchmark_sample_count = 5000L,
    benchmark_training_count = 3500L,
    benchmark_test_count = 1500L,
    benchmark_marker_count = 37991L,
    oracle_tolerance = 1e-10,
    operator_probe_count = 5L,
    projected_variance_contract =
      "yy - sum_b crossprod(w_b) + sum_b crossprod(r_b)",
    readiness = list(max_prediction_difference = 0.05,
      max_h2_difference = 0.05, min_marker_correlation = 0.95)
  ),
  frozen_capsule = "results/reference/05_ld_operator/current",
  components = list(
    operator_validation = list(
      scenarios = names(.study05_architectures),
      configurations = .study05_configurations,
      replicate_count = 5L,
      retained_mass = c(near_full = 1 - .Machine$double.eps,
        high_retention = 0.999, canonical = 0.995)
    ),
    sbayesr_ld_sensitivity = list(
      component = "sbayesr_ld_sensitivity",
      scheduler_variants = c("bed_scheduled", "bed_full_sweep",
        "csr_current"),
      operator_variants = c("bed_exact_data", "csr_exact_ld",
        "csr_sparse_ld", "block_eigen_full_rank",
        "block_eigen_retained"),
      marker_window = list(count = 1500L, source_rows = 13392:14891),
      sample_count = 5000L,
      simulation_seed = 17002L,
      causal_markers = 15L,
      chain_seeds = c(150104L, 250104L, 350104L, 450104L),
      sparse = list(r2 = 0.001, window = 1000L),
      exact = list(r2 = 0, window = 0L),
      block_size = 250L,
      block_starts = c(1L, 251L, 501L, 624L, 751L, 1001L, 1251L),
      retained_mass = 0.995,
      retained_total_rank = 1490L,
      retained_block_ranks = c(248L, 248L, 123L, 127L, 248L, 248L, 248L),
      audits = c("corrected_marker_score", "component_probability",
        "conditional_effect", "quadratic_form", "residual_expression",
        "spectral_projection"),
      controls = list(method = "sbayesr", h2 = 0.30, adjE = 0.9,
        pi = c(0.99, rep(0.01 / 3, 3L)),
        alpha = c(0.99, rep(0.01 / 3, 3L)) * 5e5,
        mixture_var = c(0, 0.01, 0.1, 1),
        updateB = TRUE, updateE = TRUE, updatePi = TRUE,
        updateLDswap = FALSE, nburn = 250L, nit = 1000L, nthin = 1L,
        nchains = 4L, ncores = 4L, seed = 50104L,
        convergence = "extended", extended_groups = "probability")
    )
  ),
  packages = list(
    sblr = list(version = "0.2.0",
      sha = "02e8c74baa906e83c4a08d42a9cc6339b4e81072"),
    numerical_fit_sblr_sha = "bd8e2c8148a0d9540dc20716455706beeb16fa86"
  )
)

rm(.study05_configurations, .study05_architectures, .study05_methods)
