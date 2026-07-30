# ============================================================
# Four-chain convergence diagnostics: compact worked example
# ============================================================
# These short settings demonstrate the workflow. They are not recommendations.

library(qgg)
library(sblr)
library(sblrbench)
library(posterior)

source(file.path("studies", "01_finemapping", "setup_example_data.R"))
source(file.path("studies", "03_parameter_estimation", "simulation.R"))
source(file.path("studies", "04_convergence", "chain_extraction.R"))
study03_config <- source(file.path("studies", "03_parameter_estimation", "config.R"), local = TRUE)$value
data_directory <- Sys.getenv("SBLR_BENCH_DATA_DIR", file.path("results", "local", "example_data"))
example_files <- .study01_example_files(data_directory, study03_config$example_data)
Glist <- .study01_load_glist(list(glist_path = "", data_dir = data_directory), example_files)

sample_ids <- Glist$ids[seq_len(200L)]
marker_ids <- Glist$rsids[[1L]][seq_len(200L)]
Glist <- .study01_set_rsids_ld(Glist, 1L, marker_ids)
genotypes <- qgg::getG(Glist = Glist, chr = 1L, ids = sample_ids,
  rsids = marker_ids, impute = TRUE, scale = TRUE)
study03_config$simulation$n_causal <- 10L
simulation <- .study03_simulate(list(architecture = "sparse_homogeneous",
  replicate = 1L, simulation_seed = 9401L), genotypes, study03_config)
stopifnot(check_oracle_genetic_values(simulation)$ok)

fit <- sblr::stblr_bed(Glist = Glist, y = simulation$truth$phenotypes,
  chr = 1L, rows = seq_along(sample_ids), method = "bayesc", h2 = 0.30,
  pi_init = 0.05, nit = 100L, nburn = 0L, nthin = 1L,
  nchains = 4L, ncores = 4L, seed = 9501L,
  chain_seeds = c(9511L, 9521L, 9531L, 9541L), keep_chains = TRUE,
  convergence = "core", convergence_control = list(warn = FALSE, keep_traces = TRUE),
  verbose = FALSE)

chain_draws <- .study04_extract_chain_draws(fit, "sparse_homogeneous", "st_bed_bayesc")
effect_wide <- reshape(chain_draws[c("raw_iteration", "chain", "effect_variance")],
  idvar = "raw_iteration", timevar = "chain", direction = "wide")
effect_matrix <- as.matrix(effect_wide[-1L])
diagnostic_result <- data.frame(rhat = posterior::rhat(effect_matrix),
  ess_bulk = posterior::ess_bulk(effect_matrix), ess_tail = posterior::ess_tail(effect_matrix),
  mcse_mean = posterior::mcse_mean(effect_matrix))
print(diagnostic_result)

chain_draws$running_mean <- ave(chain_draws$effect_variance, chain_draws$chain,
  FUN = function(x) cumsum(x) / seq_along(x))
print(utils::head(chain_draws[c("chain", "raw_iteration", "effect_variance", "running_mean")]))
