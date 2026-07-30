list(
  study = "02_prediction",
  task = "prediction",
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
  traits = c("trait1", "trait2", "trait3"),
  replicate_counts = c(development = 1L, intermediate = 5L, pilot = 10L),
  split = list(train_fraction = 0.70, seed = 3101L),
  qc = list(excludeMAF = 0.05, excludeMISS = 0.05, excludeCGAT = TRUE,
    excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
    excludeMHC = FALSE),
  sparse_ld = list(max_distance_bp = 0, max_distance_variants = 1000L,
    r2_threshold = 0.001, block_size = 1024L, nthreads = 1L),
  simulation = list(
    h2 = c(0.30, 0.30, 0.30),
    shared_effect_correlation = 0.60,
    residual_correlation = 0,
    base_seed = 4001L,
    architectures = list(
      mostly_shared = list(n_shared = 20L, n_specific = 5L),
      mostly_trait_specific = list(n_shared = 5L, n_specific = 20L)
    )
  ),
  methods = c("st_bed_bayesr", "mt_bed_bayesr", "st_csr_sbayesr", "mt_csr_sbayesr"),
  bayesr = list(mixture_var = c(0, 0.01, 0.1, 1), h2 = c(0.30, 0.30, 0.30),
    active_probability = 0.001),
  mcmc = list(nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L,
    ncores = 1L, convergence = "none", seed_offset = 10000L),
  oracle_tolerance = 1e-10
)
