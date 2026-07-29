# ============================================================
# Fine-mapping with public simulated genotype data
# ============================================================
# Tested with qgg 1.1.6, sblr 0.1.2, and sblrbench 0.0.0.9000.
# The short MCMC settings demonstrate workflow use. They are not scientific
# convergence recommendations.

# 1. Load packages ---------------------------------------------------------
library(qgg)
library(sblr)
library(sblrbench)

# 2. Download the example genotype data -----------------------------------
data_directory <- Sys.getenv(
  "SBLRBENCH_EXAMPLE_DATA_DIR",
  file.path("results", "local", "example_data")
)

example_files <- download_sblrbench_example_data(
  destination = data_directory,
  overwrite = FALSE,
  quiet = FALSE
)

# 3. Create and inspect the Glist ------------------------------------------
genotype_list <- qgg::gprep(
  study = "sblrbench worked example",
  bedfiles = unname(example_files["human.bed"]),
  bimfiles = unname(example_files["human.bim"]),
  famfiles = unname(example_files["human.fam"])
)

stopifnot(length(genotype_list$ids) >= 400L)
stopifnot(length(genotype_list$rsids[[1L]]) >= 500L)

# 4. Select samples and markers --------------------------------------------
retained_marker_ids <- qgg::gfilter(
  Glist = genotype_list,
  excludeMAF = 0.05,
  excludeMISS = 0.05,
  excludeCGAT = TRUE,
  excludeINDEL = TRUE,
  excludeDUPS = TRUE,
  excludeHWE = 1e-12,
  excludeMHC = FALSE
)

chromosome_marker_ids <- genotype_list$rsids[[1L]]
qc_marker_ids <- chromosome_marker_ids[
  !is.na(match(chromosome_marker_ids, retained_marker_ids))
]

sample_ids <- genotype_list$ids[seq_len(400L)]
marker_ids <- qc_marker_ids[seq_len(500L)]

stopifnot(!anyDuplicated(sample_ids))
stopifnot(!anyDuplicated(marker_ids))

# Tell the public BED model which chromosome-1 markers to analyse.
analysis_glist <- genotype_list
analysis_glist$rsidsLD <- vector("list", length(genotype_list$rsids))
analysis_glist$rsidsLD[[1L]] <- marker_ids

# 5. Extract standardized genotypes ----------------------------------------
genotype_matrix <- qgg::getG(
  Glist = analysis_glist,
  chr = 1L,
  ids = sample_ids,
  rsids = marker_ids,
  impute = TRUE,
  scale = TRUE
)

stopifnot(identical(rownames(genotype_matrix), sample_ids))
stopifnot(identical(colnames(genotype_matrix), marker_ids))
stopifnot(all(is.finite(genotype_matrix)))

# 6. Choose causal markers and simulate a phenotype -----------------------
causal_marker_ids <- marker_ids[c(50L, 150L, 250L, 350L, 450L)]
simulation_seed <- 20260729L

simulation_raw <- sblr::mtsim(
  W = genotype_matrix,
  standardize_W = FALSE,
  nt = 1L,
  n_shared = length(causal_marker_ids),
  n_specific = 0L,
  causal_rsids = causal_marker_ids,
  h2 = 0.2,
  seed = simulation_seed
)

stopifnot(identical(sort(simulation_raw$causal$all), sort(causal_marker_ids)))

# 7. Validate the simulated truth ------------------------------------------
simulation <- sblrbench::as_sblrbench_simulation(
  simulation_raw,
  study = "worked_finemapping_example",
  architecture = "separated",
  replicate = 1L,
  seed = simulation_seed
)

oracle <- sblrbench::check_oracle_genetic_values(
  simulation,
  tolerance = 1e-10,
  stop_on_failure = TRUE
)
stopifnot(oracle$ok)

# 8. Prepare model input ----------------------------------------------------
phenotype <- matrix(
  simulation_raw$y,
  ncol = 1L,
  dimnames = list(sample_ids, simulation$data$trait_names)
)

# 9. Fit one real public sblr model ----------------------------------------
fit_start <- proc.time()[["elapsed"]]

native_fit <- sblr::stblr_bed(
  y = phenotype,
  Glist = analysis_glist,
  method = "bayesc",
  nit = 100L,
  nburn = 50L,
  nthin = 1L,
  seed = 12001L,
  nchains = 1L,
  ncores = 1L,
  convergence = "core",
  convergence_control = list(warn = FALSE),
  verbose = FALSE
)

elapsed_seconds <- unname(proc.time()[["elapsed"]] - fit_start)

# 10. Extract posterior effects and PIPs -----------------------------------
result <- sblrbench::extract_sblr_result(
  run = list(
    native_fit = native_fit,
    elapsed_seconds = elapsed_seconds
  ),
  method_id = "worked_st_bed_bayesc",
  keep_native_fit = FALSE
)

sblrbench::validate_sblrbench_result(result, simulation)
stopifnot(identical(rownames(result$estimates$effects), marker_ids))
stopifnot(all(is.finite(result$estimates$effects)))
stopifnot(all(is.finite(result$estimates$pip)))
stopifnot(all(result$estimates$pip >= 0 & result$estimates$pip <= 1))

# 11. Evaluate recovery of the simulated truth ----------------------------
metric_rows <- sblrbench::evaluate_metrics(
  simulation = simulation,
  result = result,
  metrics = c(
    "effect_rmse",
    "pip_brier",
    "average_precision",
    "causal_ranks"
  )
)

requested_metrics <- c(
  "effect_rmse",
  "pip_brier",
  "average_precision",
  "causal_top_10_recall",
  "causal_rank_median"
)

result_table <- metric_rows[
  match(requested_metrics, metric_rows$metric),
  c("metric", "value")
]
result_table <- rbind(
  result_table,
  data.frame(metric = "runtime_seconds", value = elapsed_seconds)
)

stopifnot(all(is.finite(result_table$value)))
print(result_table, row.names = FALSE)

# 12. Inspect top-ranked markers -------------------------------------------
pip <- result$estimates$pip[, 1L]
estimated_effect <- result$estimates$effects[, 1L]
true_effect <- simulation$truth$effects[, 1L]
ranking <- order(-pip, seq_along(pip))
top_indices <- ranking[seq_len(10L)]

top_markers <- data.frame(
  marker = marker_ids[top_indices],
  PIP = pip[top_indices],
  estimated_effect = estimated_effect[top_indices],
  true_effect = true_effect[top_indices],
  causal = marker_ids[top_indices] %in% causal_marker_ids,
  row.names = NULL
)

print(top_markers, row.names = FALSE)

invisible(list(
  oracle = oracle,
  metrics = result_table,
  top_markers = top_markers,
  result = result
))
