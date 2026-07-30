.study04_project_root <- if (file.exists(file.path("studies", "03_parameter_estimation", "config.R"))) "." else file.path("..", "..")
.study04_study03_config <- source(file.path(.study04_project_root, "studies", "03_parameter_estimation", "config.R"), local = TRUE)$value

list(
  study = "04_convergence", task = "single_trait_multichain_convergence",
  development_settings = TRUE, chr = 1L, trait = "trait1", sample_limit = NULL,
  example_data = .study04_study03_config$example_data,
  qc = .study04_study03_config$qc,
  sparse_ld = .study04_study03_config$sparse_ld,
  simulation = .study04_study03_config$simulation,
  methods = c("st_bed_bayesc", "st_bed_bayesr", "st_csr_sbayesc", "st_csr_sbayesr"),
  matched_grid = data.frame(architecture = c("sparse_homogeneous", "sparse_homogeneous",
    "sparse_mixture", "sparse_mixture"), method = c("st_bed_bayesc", "st_csr_sbayesc",
    "st_bed_bayesr", "st_csr_sbayesr"), stringsAsFactors = FALSE),
  profiles = list(development = list(nit = 3000L, nburn = 0L, nthin = 1L,
    nchains = 4L, ncores = 4L, convergence = "core", keep_chains = TRUE,
    convergence_control = list(warn = FALSE, keep_traces = TRUE)),
    robustness_unexecuted = list(enabled = FALSE, design = "all eight architecture-method combinations")),
  burnin_candidates = c(250L, 500L, 1000L),
  retained_draw_candidates = c(250L, 500L, 1000L, 2000L),
  thresholds = list(rhat = 1.01, ess_bulk = 400, ess_tail = 400,
    relative_mcse = 0.05, chain_count = 4L, standardized_mean_shift = 0.10),
  seeds = list(base_fit = 40000L, chain_stride = 100L),
  priors = .study04_study03_config$priors,
  oracle_tolerance = 1e-10)
