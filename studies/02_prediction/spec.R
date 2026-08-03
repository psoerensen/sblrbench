# Exact scientific specification for the completed Study 02 prediction study.
# The `benchmark` profile reproduces the frozen five-replicate design. The
# shorter `workshop` profile is instructional and unsuitable for performance
# claims.

spec <- list(
  study = "02_prediction",
  task = "single_trait_prediction",
  supported_profiles = list(
    workshop = list(replicate_count = 1L,
      suitable_for_method_performance_claims = FALSE),
    benchmark = list(replicate_count = 5L,
      suitable_for_method_performance_claims = FALSE)
  ),
  data = list(
    source = "qgg_example",
    chromosome = 1L,
    sample_limit = NULL,
    trait = "trait1",
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
    sparse_ld = list(max_distance_bp = 0,
      max_distance_variants = 1000L, r2_threshold = 0.001,
      block_size = 1024L, nthreads = 1L)
  ),
  split = list(train_fraction = 0.70, seed = 3101L),
  markers = list(qc = list(excludeMAF = 0.05, excludeMISS = 0.05,
    excludeCGAT = TRUE, excludeINDEL = TRUE, excludeDUPS = TRUE,
    excludeHWE = 1e-12, excludeMHC = FALSE)),
  replicate_count = 5L,
  scenarios = list(
    sparse_homogeneous = list(effect_distribution = "single_normal"),
    sparse_mixture = list(effect_distribution = "variance_mixture",
      mixture_var = c(0.01, 0.1, 1),
      mixture_prob = c(0.60, 0.30, 0.10))
  ),
  methods = list(
    st_bed_bayesc = list(label = "ST-BED BayesC", interface = "stblr_bed",
      native_method = "bayesc", representation = "BED",
      prior_class = "BayesC"),
    st_bed_bayesr = list(label = "ST-BED BayesR", interface = "stblr_bed",
      native_method = "bayesr", representation = "BED",
      prior_class = "BayesR"),
    st_csr_sbayesc = list(label = "ST-CSR SBayesC", interface = "stblr_csr",
      native_method = "sbayesc", representation = "CSR",
      prior_class = "BayesC"),
    st_csr_sbayesr = list(label = "ST-CSR SBayesR", interface = "stblr_csr",
      native_method = "sbayesr", representation = "CSR",
      prior_class = "BayesR")
  ),
  controls = list(
    simulation = list(h2 = 0.30, n_causal = 50L),
    priors = list(h2 = 0.30, bayesc_inclusion_probability = 0.01,
      bayesr_active_probability = 0.01,
      bayesr_mixture_var = c(0, 0.01, 0.1, 1)),
    benchmark = list(nchains = 4L,
      recommendation_source = file.path("results", "reference",
        "04_convergence", "current-selection", "method_recommendations.csv")),
    workshop = list(nit = 100L, nburn = 50L, nthin = 1L, nchains = 1L,
      ncores = 1L, convergence = "none", keep_chains = TRUE,
      convergence_control = list(warn = FALSE, keep_traces = TRUE))
  ),
  seeds = list(simulation_base = 4001L, fit_offset = 10000L,
    chain_stride = 100000L),
  metrics = c("prediction_correlation", "prediction_mse",
    "prediction_nmse", "phenotype_prediction_correlation",
    "prediction_calibration", "effect_rmse"),
  validation = list(oracle_tolerance = 1e-10,
    benchmark_sample_count = 5000L, benchmark_training_count = 3500L,
    benchmark_test_count = 1500L, benchmark_marker_count = 37991L,
    expected_fit_count = 40L),
  frozen_capsule = file.path("results", "reference", "02_prediction",
    "current"),
  packages = list(
    sblr = list(version = "0.2.0",
      sha = "02e8c74baa906e83c4a08d42a9cc6339b4e81072"),
    qgg = list(version = "1.1.6"),
    sblrbench = list(version = "0.0.0.9000")
  )
)
