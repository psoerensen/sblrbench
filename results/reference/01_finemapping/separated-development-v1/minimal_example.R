# Minimal, fast contract demonstration; no large Glist or sampler fit is required.
library(sblr)
library(sblrbench)

set.seed(42)
W <- matrix(sample(0:2, 240, replace = TRUE), nrow = 30, ncol = 8,
  dimnames = list(paste0("sample", 1:30), paste0("marker", 1:8)))
sim_raw <- sblr::mtsim(W = W, nt = 1L, n_shared = 2L, n_specific = 0L,
  h2 = 0.2, seed = 2001L)
simulation <- as_sblrbench_simulation(sim_raw, study = "minimal_example",
  architecture = "synthetic", replicate = 1L, seed = 2001L)
check_oracle_genetic_values(simulation)

# Lightweight standard-result demonstration using known truth.
causal <- matrix(as.numeric(rownames(simulation$truth$effects) %in% sim_raw$causal$all),
  ncol = 1L, dimnames = dimnames(simulation$truth$effects))
result <- new_sblrbench_result("truth_demo", effects = simulation$truth$effects, pip = causal)
metric_effect_rmse(simulation, result)
metric_pip_brier(simulation, result)
# A full public sblr sampler requires its corresponding BED/Glist or CSR/LD input;
