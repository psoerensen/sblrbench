spec <- list(
  study = "01_finemapping",
  task = "finemapping",
  supported_profiles = list(
    workshop = list(replicate_count = 1L),
    benchmark = list(replicate_count = 10L)
  ),
  data = list(
    source = "qgg_example", chromosome = 1L, sample_limit = NULL,
    example_data = list(
      repository = "psoerensen/qgdata",
      commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
      subdirectory = "simulated_human_data",
      files = c("human.bed", "human.bim", "human.fam", "human.pheno",
        "human.covar")
    ),
    sparse_ld = list(max_distance_bp = 0, max_distance_variants = 1000L,
      r2_threshold = 0.001, block_size = 1024L, nthreads = 1L)
  ),
  markers = list(qc = list(excludeMAF = 0.05, excludeMISS = 0.05,
    excludeCGAT = TRUE, excludeINDEL = TRUE, excludeDUPS = TRUE,
    excludeHWE = 1e-12, excludeMHC = FALSE)),
  replicate_count = 10L,
  scenarios = list(separated = list(effect_distribution = "sblr_mtsim",
    target_heritability = 0.2, causal_markers = 10L)),
  causal_design = list(min_distance_bp = 1e6, min_maf = 0.05,
    max_maf = 0.5, selection = "seeded greedy selection in marker order"),
  locus_design = list(
    implementation = "studies/01_finemapping/locus-design.R",
    credible_set_target = 0.95, min_r2 = 0.5, pip_cutoff = 0.001,
    locus_pip_cutoff = 0.01, max_locus_distance = 1e6,
    algorithm = "sblr::make_credible_sets"
  ),
  methods = list(
    st_bed_bayesc = list(label = "ST-BED BayesC", interface = "stblr_bed",
      native_method = "bayesc", representation = "BED", prior_class = "BayesC"),
    st_bed_bayesr = list(label = "ST-BED BayesR", interface = "stblr_bed",
      native_method = "bayesr", representation = "BED", prior_class = "BayesR"),
    st_csr_sbayesc = list(label = "ST-CSR SBayesC", interface = "stblr_csr",
      native_method = "sbayesc", representation = "CSR", prior_class = "BayesC"),
    st_csr_sbayesr = list(label = "ST-CSR SBayesR", interface = "stblr_csr",
      native_method = "sbayesr", representation = "CSR", prior_class = "BayesR")
  ),
  controls = list(
    simulation = list(h2 = 0.2, n_causal = 10L),
    workshop = list(nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L,
      ncores = 1L, convergence = "core", keep_chains = FALSE,
      convergence_control = list(warn = FALSE)),
    benchmark = list(nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L,
      ncores = 1L, convergence = "core", keep_chains = FALSE,
      convergence_control = list(warn = FALSE))
  ),
  seeds = list(simulation_base = 2001L, causal_offset = 0L,
    phenotype_offset = 1000L, fit_offset = 10000L, chain_stride = 1L),
  estimands = data.frame(
    estimand = c("marker_pip", "causal_marker_rank", "credible_set_size",
      "credible_set_coverage"),
    definition = c("posterior marginal inclusion probability",
      "descending PIP rank of each causal marker",
      "number of markers in the native marginal-PIP credible set",
      "exact and r-squared proxy coverage of the matched causal marker"),
    stringsAsFactors = FALSE),
  metrics = c("pip_brier", "effect_rmse", "average_precision",
    "causal_ranks", "credible_set_summary"),
  validation = list(oracle_tolerance = 1e-10,
    benchmark_sample_count = 5000L, benchmark_marker_count = 37991L,
    expected_fit_count = 40L),
  frozen_capsule = "results/reference/01_finemapping/current",
  packages = list(
    sblr = list(version = "0.2.0",
      sha = "02e8c74baa906e83c4a08d42a9cc6339b4e81072"),
    qgdata = list(sha = "6cca5819e711d326cfb2614d7e9d9f34942612cd")
  )
)
