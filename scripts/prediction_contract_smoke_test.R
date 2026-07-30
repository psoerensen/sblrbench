# Developer prediction contract smoke test
# No sampler is fitted. This checks sample subsetting, prediction attachment,
# strict alignment, and hand-verifiable prediction metrics.
library(sblrbench)

sample_ids <- paste0("sample", 1:6)
marker_ids <- paste0("marker", 1:3)
trait_names <- c("trait1", "trait2")
genotypes <- matrix(seq_len(18) / 10, 6, 3,
  dimnames = list(sample_ids, marker_ids))
effects <- matrix(c(.2, -.1, 0, .1, .05, -.2), 3, 2,
  dimnames = list(marker_ids, trait_names))
genetic_values <- genotypes %*% effects
phenotypes <- genetic_values + matrix(seq_len(12) / 100, 6, 2,
  dimnames = list(sample_ids, trait_names))

simulation <- list(schema_version = 1L,
  data = list(marker_ids = marker_ids, sample_ids = sample_ids,
    trait_names = trait_names, train_ids = sample_ids[1:4],
    test_ids = sample_ids[5:6], reference_ids = NULL, genotypes = genotypes),
  truth = list(effects = effects, genetic_values = genetic_values,
    phenotypes = phenotypes, residuals = phenotypes - genetic_values,
    causal = list(shared = marker_ids[1], specific = list(trait1 = marker_ids[2], trait2 = marker_ids[3])),
    parameters = list()), scenario = list(study = "contract_smoke",
      architecture = "fabricated", replicate = 1L),
  provenance = list(seed = 1L), extras = list())
class(simulation) <- c("sblrbench_simulation", "list")
validate_sblrbench_simulation(simulation)

test_simulation <- subset_sblrbench_simulation_samples(simulation, sample_ids[5:6])
result <- new_sblrbench_result("oracle", effects = effects)
result <- add_sblrbench_predictions(result,
  test_simulation$truth$genetic_values, test_simulation)
metrics <- evaluate_metrics(test_simulation, result, c("prediction_correlation",
  "prediction_mse", "prediction_nmse", "phenotype_prediction_correlation",
  "prediction_calibration"))
stopifnot(all(metrics$status == "ok"), all(is.finite(metrics$value)))
print(metrics)
