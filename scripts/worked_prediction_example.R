# ============================================================
# Predicting three traits with public simulated genotype data
# ============================================================
# Demonstration settings: 100 iterations, 50 burn-in, one chain, one core.
# These settings illustrate the workflow; they are not convergence guidance.

# 1. Load packages
library(qgg)
library(sblr)
library(sblrbench)

set.seed(5201)
data_directory <- Sys.getenv(
  "SBLRBENCH_EXAMPLE_DATA_DIR",
  file.path("results", "local", "worked_prediction_example")
)

# 2. Download the public genotype data
example_files <- download_sblrbench_example_data(data_directory)

# 3. Create the Glist
glist_cache <- file.path(data_directory, "human_glist.rds")
if (file.exists(glist_cache)) {
  genotype_list <- readRDS(glist_cache)
} else {
  genotype_list <- qgg::gprep(
    study = "sblrbench prediction example",
    bedfiles = unname(example_files["human.bed"]),
    bimfiles = unname(example_files["human.bim"]),
    famfiles = unname(example_files["human.fam"])
  )
  saveRDS(genotype_list, glist_cache)
}

# 4. Select samples and markers
qc_markers <- qgg::gfilter(
  Glist = genotype_list,
  excludeMAF = 0.05,
  excludeMISS = 0.05,
  excludeCGAT = TRUE,
  excludeINDEL = TRUE,
  excludeDUPS = TRUE,
  excludeHWE = 1e-12,
  excludeMHC = FALSE
)

chromosome <- 1L
sample_ids <- genotype_list$ids[seq_len(min(600L, length(genotype_list$ids)))]
chromosome_markers <- genotype_list$rsids[[chromosome]]
marker_ids <- chromosome_markers[!is.na(match(chromosome_markers, qc_markers))]
marker_ids <- marker_ids[seq_len(min(600L, length(marker_ids)))]

# 5. Create the train/test split
split <- make_prediction_split(
  sample_ids = sample_ids,
  train_fraction = 0.70,
  seed = 3101L
)

# 6. Learn genotype scaling from training data
raw_genotypes <- qgg::getG(
  Glist = genotype_list,
  chr = chromosome,
  ids = sample_ids,
  rsids = marker_ids,
  impute = FALSE,
  scale = FALSE
)

scaled <- training_scaled_genotypes(
  raw_genotypes = raw_genotypes,
  train_rows = split$train_rows
)

working_glist <- genotype_list
marker_columns <- match(marker_ids, working_glist$rsids[[chromosome]])
working_glist$af[[chromosome]][marker_columns] <- scaled$allele_frequency
working_glist$maf[[chromosome]][marker_columns] <- pmin(
  scaled$allele_frequency,
  1 - scaled$allele_frequency
)
working_glist$rsidsLD[[chromosome]] <- marker_ids

# 7. Simulate three correlated traits
trait_names <- c("trait1", "trait2", "trait3")
genetic_correlation <- matrix(0.60, 3L, 3L)
diag(genetic_correlation) <- 1

simulation_raw <- sblr::mtsim(
  W = scaled$all,
  standardize_W = FALSE,
  nt = 3L,
  n_shared = 20L,
  n_specific = 5L,
  h2 = rep(0.30, 3L),
  rg = genetic_correlation,
  re = 0,
  seed = 5201L
)

colnames(simulation_raw$y) <- trait_names
colnames(simulation_raw$B) <- trait_names
colnames(simulation_raw$G) <- trait_names
colnames(simulation_raw$E) <- trait_names

simulation <- as_sblrbench_simulation(
  simulation_raw,
  study = "worked_prediction_example",
  architecture = "mostly_shared",
  replicate = 1L,
  seed = 5201L
)
simulation$data$train_ids <- split$train_ids
simulation$data$test_ids <- split$test_ids
validate_sblrbench_simulation(simulation)

# 8. Validate the simulation
oracle <- check_oracle_genetic_values(simulation)
print(oracle)

training_phenotypes <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
test_simulation <- subset_sblrbench_simulation_samples(simulation, split$test_ids)

# 9. Fit the single-trait model
fit_controls <- list(
  method = "bayesr",
  h2 = 0.30,
  pi = c(0.999, 0.001 / 3, 0.001 / 3, 0.001 / 3),
  mixture_var = c(0, 0.01, 0.1, 1),
  nit = 100L,
  nburn = 50L,
  nthin = 1L,
  nchains = 1L,
  ncores = 1L,
  convergence = "none",
  verbose = FALSE
)

single_start <- proc.time()[["elapsed"]]
single_fits <- lapply(seq_along(trait_names), function(j) {
  do.call(sblr::stblr_bed, c(list(
    y = training_phenotypes[, j, drop = FALSE],
    Glist = working_glist,
    rows = split$train_rows,
    seed = 6100L + j
  ), fit_controls))
})
single_runtime <- proc.time()[["elapsed"]] - single_start

single_effects <- do.call(cbind, lapply(single_fits, function(fit) fit$bm[, 1L]))
rownames(single_effects) <- marker_ids
colnames(single_effects) <- trait_names

# 10. Fit the multi-trait model
multi_start <- proc.time()[["elapsed"]]
multi_fit <- sblr::mtblr_bed(
  y = training_phenotypes,
  Glist = working_glist,
  rows = split$train_rows,
  scale = TRUE,
  residual_covariance = "diagonal",
  method = "bayesr",
  h2 = rep(0.30, 3L),
  pi = 0.001,
  mixture_var = c(0, 0.01, 0.1, 1),
  nit = 100L,
  nburn = 50L,
  nthin = 1L,
  seed = 6201L,
  nchains = 1L,
  ncores = 1L,
  convergence = "none",
  verbose = FALSE
)
multi_runtime <- proc.time()[["elapsed"]] - multi_start
multi_effects <- multi_fit$bm

# 11. Predict genetic values in the test set
make_result <- function(method_id, effects, runtime) {
  result <- new_sblrbench_result(
    method_id = method_id,
    effects = effects,
    elapsed_seconds = runtime,
    keep_native_fit = FALSE
  )
  predictions <- scaled$test %*% effects
  add_sblrbench_predictions(result, predictions, test_simulation)
}

single_result <- make_result("st_bed_bayesr", single_effects, single_runtime)
multi_result <- make_result("mt_bed_bayesr", multi_effects, multi_runtime)

# 12. Compare predictive performance
metric_names <- c(
  "prediction_correlation",
  "prediction_mse",
  "prediction_nmse",
  "phenotype_prediction_correlation",
  "prediction_calibration"
)
single_metrics <- evaluate_metrics(test_simulation, single_result, metric_names)
multi_metrics <- evaluate_metrics(test_simulation, multi_result, metric_names)
single_metrics$model <- "ST-BED BayesR"
multi_metrics$model <- "MT-BED BayesR"
comparison <- rbind(single_metrics, multi_metrics)

print(comparison[, c("model", "trait", "metric", "value", "status")], row.names = FALSE)

paired <- merge(
  single_metrics[, c("trait", "metric", "value")],
  multi_metrics[, c("trait", "metric", "value")],
  by = c("trait", "metric"),
  suffixes = c("_st", "_mt")
)
error_metric <- paired$metric %in% c("prediction_mse", "prediction_nmse")
paired$mt_advantage <- ifelse(error_metric, paired$value_st - paired$value_mt,
                              paired$value_mt - paired$value_st)
print(paired, row.names = FALSE)
cat(sprintf("\nST runtime: %.2f seconds; MT runtime: %.2f seconds\n",
            single_runtime, multi_runtime))
