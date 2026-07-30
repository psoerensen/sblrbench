list(
  study = "02_prediction",
  task = "single_trait_prediction",
  development_settings = TRUE,
  data_source = "qgg_example",
  example_data = list(
    repository = "psoerensen/qgdata",
    commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
    subdirectory = "simulated_human_data",
    files = c("human.bed", "human.bim", "human.fam", "human.pheno", "human.covar")
  ),
  chr = 1L,
  sample_limit = NULL,
  trait = "trait1",
  replicate_counts = c(development = 1L, intermediate = 5L, pilot = 10L),
  split = list(train_fraction = 0.70, seed = 3101L),
  qc = list(excludeMAF = 0.05, excludeMISS = 0.05, excludeCGAT = TRUE,
    excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
    excludeMHC = FALSE),
  sparse_ld = list(max_distance_bp = 0, max_distance_variants = 1000L,
    r2_threshold = 0.001, block_size = 1024L, nthreads = 1L),
  simulation = list(
    h2 = 0.30,
    n_causal = 50L,
    base_seed = 4001L,
    architectures = list(
      sparse_homogeneous = list(effect_distribution = "single_normal"),
      sparse_mixture = list(effect_distribution = "variance_mixture",
        mixture_var = c(0.01, 0.1, 1), mixture_prob = c(0.60, 0.30, 0.10))
    )
  ),
  methods = c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"),
  priors = list(
    h2 = 0.30,
    bayesc_inclusion_probability = 0.01,
    bayesr_active_probability = 0.01,
    bayesr_mixture_var = c(0, 0.01, 0.1, 1)
  ),
  multitrait = list(
    enabled = FALSE,
    status = "deferred",
    reason = paste("Multi-trait BED and CSR models are excluded from the current",
      "benchmark pending computational profiling and optimization."),
    methods = c("mt_bed_bayesr", "mt_csr_sbayesr")
  ),
  mcmc = list(nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L,
    ncores = 1L, convergence = "none", seed_offset = 10000L),
  oracle_tolerance = 1e-10
)
