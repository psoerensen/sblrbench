# Developer convergence contract smoke test
#
# This script probes the verified native two-chain retention contract with a
# deliberately tiny public-data subset. It is not benchmark evidence and its
# short chains are not convergence recommendations.

library(qgg)
library(sblr)
library(sblrbench)

source(file.path("studies", "01_finemapping", "setup_example_data.R"))
source(file.path("studies", "04_convergence", "chain_extraction.R"))

config <- read_benchmark_spec(file.path("studies","03_parameter_estimation","spec.R"))
paths <- list(glist_path = "", data_dir = file.path("results", "local", "03_parameter_estimation", "data"))
files <- .study01_example_files(paths$data_dir, config$data$example_data)
Glist <- .study01_load_glist(paths, files)
sample_ids <- Glist$ids[seq_len(min(120L, length(Glist$ids)))]
marker_ids <- Glist$rsids[[1L]][seq_len(min(120L, length(Glist$rsids[[1L]])))]
Glist <- .study01_set_rsids_ld(Glist, 1L, marker_ids)
Z <- qgg::getG(Glist = Glist, chr = 1L, ids = sample_ids, rsids = marker_ids,
  impute = TRUE, scale = TRUE)

spec <- list(architecture = "sparse_homogeneous", replicate = 1L, simulation_seed = 991L)
config$controls$simulation$n_causal <- 8L
simulation <- sblrbench:::simulate_prediction_architecture(list(
  scenario=spec$architecture,replicate=spec$replicate,
  simulation_seed=spec$simulation_seed),Z,config)
stopifnot(check_oracle_genetic_values(simulation)$ok)

fit <- sblr::stblr_bed(
  Glist = Glist, y = simulation$truth$phenotypes, chr = 1L,
  rows = seq_along(sample_ids), method = "bayesc", h2 = 0.30,
  pi_init = 0.05, nit = 30L, nburn = 0L, nthin = 1L,
  nchains = 2L, ncores = 2L, seed = 8001L,
  chain_seeds = c(8101L, 8201L), keep_chains = TRUE,
  convergence = "core", convergence_control = list(warn = FALSE, keep_traces = TRUE),
  verbose = FALSE
)

draws <- .study04_extract_chain_draws(fit, spec$architecture, "st_bed_bayesc", expected_chains = 2L)
stopifnot(identical(sort(unique(draws$chain)), 1:2), all(table(draws$chain) == 30L))
cat("Tiny native multichain contract: 2 chains x 30 retained raw draws; chain identity preserved.\n")
