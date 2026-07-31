list(
  study = "02_prediction",
  task = "single_trait_prediction",
  profile = Sys.getenv("SBLR_BENCH_PROFILE", "one_replicate_development"),
  development_settings = TRUE,
  data_source = "qgg_example",
  example_data = list(
    repository = "psoerensen/qgdata",
    commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
    subdirectory = "simulated_human_data",
    files = c("human.bed", "human.bim", "human.fam", "human.pheno", "human.covar"),
    size_bytes = c(human.bed = 62500003, human.bim = 1882359,
      human.fam = 117786, human.pheno = 92786, human.covar = 641513),
    md5 = c(human.bed = "e89bea9a6cedd9eeef3fd0a5c807db81",
      human.bim = "0105119b04c67b7ac7f66cc5e6680963",
      human.fam = "3c5db3d9eb7f3fc893c75f6f2b89836d",
      human.pheno = "6a9e7cb1162e43999c170a363863176d",
      human.covar = "d06002aa2b1b79bdc4c0e92c21f27ae5")
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
  reference_profiles = list(
    one_replicate_development = list(
      replicate_count = 1L,
      capsule_id = "st-bayesc-bayesr-one-replicate-development-v1",
      status = "complete_development_benchmark"
    ),
    five_replicate_development = list(
      replicate_count = 5L,
      capsule_id = "st-bayesc-bayesr-five-replicate-development-v1",
      status = "complete_five_replicate_development_benchmark"
    )
  ),
  priors = list(
    h2 = 0.30,
    bayesc_inclusion_probability = 0.01,
    bayesr_active_probability = 0.01,
    bayesr_mixture_var = c(0, 0.01, 0.1, 1)
  ),
  mcmc = list(nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L,
    ncores = 1L, convergence = "none", seed_offset = 10000L),
  oracle_tolerance = 1e-10
)
