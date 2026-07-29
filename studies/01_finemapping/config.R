list(
  study = "finemapping",
  task = "genotype_setup",
  data_source = "qgg_example",
  chr = 1L,
  sample_limit = NULL,
  qc = list(
    excludeMAF = 0.05,
    excludeMISS = 0.05,
    excludeCGAT = TRUE,
    excludeINDEL = TRUE,
    excludeDUPS = TRUE,
    excludeHWE = 1e-12,
    excludeMHC = FALSE
  ),
  sparse_ld = list(
    max_distance_bp = 0,
    max_distance_variants = 1000L,
    r2_threshold = 0.001,
    block_size = 1024L,
    nthreads = 1L
  ),
  simulation = list(
    nt = 1L,
    h2 = 0.2,
    n_shared = 10L,
    n_specific = 0L,
    seed = 1001L
  ),
  oracle_tolerance = 1e-10
)
