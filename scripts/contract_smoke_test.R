# ============================================================
# Developer contract smoke test
# ============================================================
# This fast truth-based check exercises simulation, oracle, result-object, and
# metric contracts without running a sampler. It is not the user tutorial.

library(sblr)
library(sblrbench)

set.seed(42)
genotype_values <- sample(0:2, 240, replace = TRUE)
genotype_matrix <- matrix(genotype_values, nrow = 30L, ncol = 8L)
rownames(genotype_matrix) <- paste0("sample", seq_len(nrow(genotype_matrix)))
colnames(genotype_matrix) <- paste0("marker", seq_len(ncol(genotype_matrix)))

simulation_raw <- sblr::mtsim(
  W = genotype_matrix,
  nt = 1L,
  n_shared = 2L,
  n_specific = 0L,
  h2 = 0.2,
  seed = 2001L
)

simulation <- sblrbench::as_sblrbench_simulation(
  simulation_raw,
  study = "contract_smoke_test",
  architecture = "synthetic",
  replicate = 1L,
  seed = 2001L
)

oracle <- sblrbench::check_oracle_genetic_values(simulation)
stopifnot(oracle$ok)

causal_indicator <- as.numeric(
  rownames(simulation$truth$effects) %in% simulation_raw$causal$all
)
causal_pip <- matrix(causal_indicator, ncol = 1L)
dimnames(causal_pip) <- dimnames(simulation$truth$effects)

truth_result <- sblrbench::new_sblrbench_result(
  method_id = "truth_contract_demo",
  effects = simulation$truth$effects,
  pip = causal_pip
)

metrics <- sblrbench::evaluate_metrics(
  simulation,
  truth_result,
  metrics = c("effect_rmse", "pip_brier")
)

stopifnot(all(metrics$status == "ok"))
stopifnot(all(metrics$value == 0))
print(metrics)
