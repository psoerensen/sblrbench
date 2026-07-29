development_mcmc <- list(
  nit = 500L, nburn = 250L, nthin = 1L, nchains = 1L, ncores = 1L,
  seed_offset = 10000L, convergence_control = list(warn = FALSE)
)

list(
  study = "01_finemapping", architecture = "separated", task = "pilot",
  data_source = "qgg_example", chr = 1L, sample_limit = NULL,
  example_data = list(
    repository = "psoerensen/qgdata",
    commit = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
    subdirectory = "simulated_human_data",
    files = c("human.bed", "human.bim", "human.fam", "human.pheno", "human.covar")
  ),
  replicate_counts = c(development = 1L, intermediate = 5L, pilot = 10L),
  qc = list(excludeMAF = 0.05, excludeMISS = 0.05, excludeCGAT = TRUE,
    excludeINDEL = TRUE, excludeDUPS = TRUE, excludeHWE = 1e-12,
    excludeMHC = FALSE),
  sparse_ld = list(max_distance_bp = 0, max_distance_variants = 1000L,
    r2_threshold = 0.001, block_size = 1024L, nthreads = 1L),
  simulation = list(h2 = 0.2, n_causal = 10L, base_seed = 2001L),
  causal_selection = list(min_distance_bp = 1e6, min_maf = 0.05, max_maf = 0.5),
  methods = c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"),
  mcmc = development_mcmc,
  credible_sets = list(coverage = 0.95, min_r2 = 0.5, pip_cutoff = 0.001,
    locus_pip_cutoff = 0.01, max_locus_distance = 1e6),
  oracle_tolerance = 1e-10
)
