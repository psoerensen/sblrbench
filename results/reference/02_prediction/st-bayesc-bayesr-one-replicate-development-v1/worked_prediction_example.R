# ============================================================
# Single-trait prediction with public simulated genotype data
# ============================================================
# Demonstration settings: 100 iterations, 50 burn-in, one chain, one core.
# They illustrate the workflow and are not convergence recommendations.

# 1. Load packages
library(qgg)
library(sblr)
library(sblrbench)

data_directory <- Sys.getenv("SBLRBENCH_EXAMPLE_DATA_DIR",
  file.path("results", "local", "worked_prediction_example"))
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)

# 2. Download public genotype data
example_files <- download_sblrbench_example_data(data_directory)

# 3. Create the Glist
glist_cache <- file.path(data_directory, "human_glist.rds")
if (file.exists(glist_cache)) {
  genotype_list <- readRDS(glist_cache)
} else {
  genotype_list <- qgg::gprep(study = "sblrbench prediction example",
    bedfiles = unname(example_files["human.bed"]),
    bimfiles = unname(example_files["human.bim"]),
    famfiles = unname(example_files["human.fam"]))
  saveRDS(genotype_list, glist_cache)
}

# 4. Select samples and markers
qc_markers <- qgg::gfilter(Glist = genotype_list, excludeMAF = 0.05,
  excludeMISS = 0.05, excludeCGAT = TRUE, excludeINDEL = TRUE,
  excludeDUPS = TRUE, excludeHWE = 1e-12, excludeMHC = FALSE)
chromosome <- 1L
sample_ids <- genotype_list$ids[seq_len(600L)]
marker_ids <- genotype_list$rsids[[chromosome]]
marker_ids <- marker_ids[!is.na(match(marker_ids, qc_markers))][seq_len(600L)]

# 5. Create the train/test split
split <- make_prediction_split(sample_ids, train_fraction = 0.70, seed = 3101L)

# 6. Learn genotype scaling from training data
raw_genotypes <- qgg::getG(Glist = genotype_list, chr = chromosome,
  ids = sample_ids, rsids = marker_ids, impute = FALSE, scale = FALSE)
scaled <- training_scaled_genotypes(raw_genotypes, split$train_rows)
working_glist <- genotype_list
marker_columns <- match(marker_ids, working_glist$rsids[[chromosome]])
working_glist$af[[chromosome]][marker_columns] <- scaled$allele_frequency
working_glist$maf[[chromosome]][marker_columns] <- pmin(scaled$allele_frequency,
  1 - scaled$allele_frequency)
working_glist$rsidsLD[[chromosome]] <- marker_ids

# 7. Simulate a sparse phenotype
# This is a variance-mixture architecture followed by one global rescaling to
# exact realized h2. It is approximately BayesR-like, not an exact prior draw.
set.seed(5201L)
trait_name <- "trait1"
causal_index <- sort(sample.int(length(marker_ids), 30L))
mixture_var <- c(0.01, 0.1, 1)
component <- sample.int(3L, 30L, replace = TRUE, prob = c(.60, .30, .10))
effects <- matrix(0, length(marker_ids), 1L,
  dimnames = list(marker_ids, trait_name))
effects[causal_index, 1L] <- rnorm(30L, sd = sqrt(mixture_var[component]))
genetic_values <- scaled$all %*% effects
effects <- effects * sqrt((.30 / .70) / var(genetic_values[, 1L]))
genetic_values <- scaled$all %*% effects
residual <- rnorm(length(sample_ids)); residual <- (residual - mean(residual)) / sd(residual)
residuals <- matrix(residual, ncol = 1L, dimnames = list(sample_ids, trait_name))
phenotypes <- genetic_values + residuals
causal_ids <- marker_ids[causal_index]
simulation_raw <- list(y = phenotypes, W = scaled$all, B = effects,
  G = genetic_values, E = residuals,
  causal = list(shared = causal_ids,
    specific = setNames(list(character()), trait_name), all = causal_ids),
  rsids = marker_ids, ids = sample_ids, h2_target = .30,
  h2_observed = var(genetic_values[, 1L]) /
    (var(genetic_values[, 1L]) + var(residuals[, 1L])))
simulation <- as_sblrbench_simulation(simulation_raw,
  study = "worked_prediction_example", architecture = "sparse_mixture",
  replicate = 1L, seed = 5201L)
simulation$data$train_ids <- split$train_ids
simulation$data$test_ids <- split$test_ids
validate_sblrbench_simulation(simulation)

# 8. Validate the simulation oracle
print(check_oracle_genetic_values(simulation))
training_y <- simulation$truth$phenotypes[split$train_ids, , drop = FALSE]
test_simulation <- subset_sblrbench_simulation_samples(simulation, split$test_ids)

development_controls <- list(nit = 100L, nburn = 50L, nthin = 1L,
  nchains = 1L, ncores = 1L, convergence = "none", verbose = FALSE,
  h2 = .30)
fit_timed <- function(fun, args) {
  start <- proc.time()[["elapsed"]]
  fit <- do.call(fun, args)
  list(fit = fit, runtime = unname(proc.time()[["elapsed"]] - start))
}

# 9. Fit BED BayesC
bed_bayesc <- fit_timed(sblr::stblr_bed, c(list(y = training_y,
  Glist = working_glist, rows = split$train_rows, method = "bayesc",
  pi_init = .02, seed = 6101L), development_controls))

# 10. Fit BED BayesR
bayesr_pi <- c(.98, .02 / 3, .02 / 3, .02 / 3)
bed_bayesr <- fit_timed(sblr::stblr_bed, c(list(y = training_y,
  Glist = working_glist, rows = split$train_rows, method = "bayesr",
  pi = bayesr_pi, mixture_var = c(0, .01, .1, 1), seed = 6102L),
  development_controls))

# 11. Build training-only sparse LD and summary statistics
ld_prefix <- file.path(data_directory, "worked_prediction_training_ld")
ld_cache <- paste0(ld_prefix, "_glist.rds")
if (file.exists(ld_cache)) {
  ld_glist <- readRDS(ld_cache)
} else {
  ld_glist <- sblr::make_sparse_ld(Glist = working_glist,
    rows = split$train_rows, out_prefix = ld_prefix, chr = chromosome,
    max_distance_variants = 200L, r2_threshold = 0.001,
    block_size = 256L, nthreads = 1L)
  saveRDS(ld_glist, ld_cache)
}
stopifnot(identical(ld_glist$sparseLD$rows, split$train_rows),
  identical(as.integer(ld_glist$sparseLD$reference_n), length(split$train_rows)))
summary_stats <- sblr::make_summary_stats(Glist = ld_glist, y = training_y,
  chr = chromosome, rows = split$train_rows, scale = TRUE, nthreads = 1L)
stopifnot(identical(summary_stats$marker_names, marker_ids),
  identical(as.integer(summary_stats$n), length(split$train_rows)))

# 12. Fit CSR SBayesC
csr_sbayesc <- fit_timed(sblr::stblr_csr, c(list(stats = summary_stats,
  Glist = ld_glist, method = "sbayesc", pi_init = .02, seed = 6103L),
  development_controls))

# 13. Fit CSR SBayesR
csr_sbayesr <- fit_timed(sblr::stblr_csr, c(list(stats = summary_stats,
  Glist = ld_glist, method = "sbayesr", pi = bayesr_pi,
  mixture_var = c(0, .01, .1, 1), seed = 6104L), development_controls))

# 14. Predict genetic values in the test set
native_runs <- list(st_bed_bayesc = bed_bayesc, st_bed_bayesr = bed_bayesr,
  st_csr_sbayesc = csr_sbayesc, st_csr_sbayesr = csr_sbayesr)
results <- lapply(names(native_runs), function(method_id) {
  run <- native_runs[[method_id]]
  result <- extract_sblr_result(list(native_fit = run$fit,
    elapsed_seconds = run$runtime), method_id = method_id,
    keep_native_fit = FALSE)
  effects_aligned <- align_markers(result$estimates$effects, marker_ids)
  predictions <- scaled$test %*% effects_aligned
  add_sblrbench_predictions(result, predictions, test_simulation)
})
names(results) <- names(native_runs)

# 15. Compare prediction metrics
metrics <- do.call(rbind, lapply(names(results), function(method_id) {
  z <- evaluate_metrics(test_simulation, results[[method_id]], c(
    "prediction_correlation", "prediction_nmse",
    "phenotype_prediction_correlation", "prediction_calibration", "effect_rmse"))
  z$method <- method_id; z
}))
summary <- reshape(metrics[, c("method", "metric", "value")],
  idvar = "method", timevar = "metric", direction = "wide")
summary$runtime <- vapply(native_runs, `[[`, numeric(1), "runtime")
print(summary, row.names = FALSE)
