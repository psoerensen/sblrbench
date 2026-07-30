#' Create a deterministic prediction train/test split
#'
#' Sampling determines membership, while the original sample order is retained
#' within both returned sets.
#' @param sample_ids Unique sample identifiers.
#' @param train_fraction Fraction assigned to training.
#' @param seed Deterministic sampling seed.
#' @export
make_prediction_split <- function(sample_ids, train_fraction = 0.7, seed = 1L) {
  sample_ids <- .assert_ids(as.character(sample_ids), "sample_ids")
  if (!is.numeric(train_fraction) || length(train_fraction) != 1L ||
      !is.finite(train_fraction) || train_fraction <= 0 || train_fraction >= 1) {
    stop("train_fraction must be a finite scalar strictly between zero and one.", call. = FALSE)
  }
  n_train <- floor(length(sample_ids) * train_fraction)
  if (n_train < 1L || n_train >= length(sample_ids)) stop("The split must contain training and test samples.", call. = FALSE)
  if (length(seed) != 1L || is.na(seed)) stop("seed must be a non-missing scalar.", call. = FALSE)
  set.seed(as.integer(seed))
  chosen <- sort(sample.int(length(sample_ids), n_train, replace = FALSE))
  test_rows <- setdiff(seq_along(sample_ids), chosen)
  list(train_ids = sample_ids[chosen], test_ids = sample_ids[test_rows],
       train_rows = chosen, test_rows = test_rows, split_seed = as.integer(seed),
       train_fraction = train_fraction,
       split_id = paste0("train", n_train, "-test", length(test_rows), "-seed", as.integer(seed)))
}

#' Learn and apply training-only genotype scaling
#'
#' Missing genotypes are imputed to the training marker mean. Genotypes are
#' centered by `2p` and divided by `sqrt(2p(1-p))`, where `p` is learned only
#' from training rows.
#' @param raw_genotypes Numeric sample-by-marker dosage matrix.
#' @param train_rows Integer training row indices.
#' @return A list containing training allele frequencies and consistently
#' scaled all/training/test matrices.
#' @export
training_scaled_genotypes <- function(raw_genotypes, train_rows) {
  if (!is.matrix(raw_genotypes) || !is.numeric(raw_genotypes) ||
      is.null(rownames(raw_genotypes)) || is.null(colnames(raw_genotypes))) {
    stop("raw_genotypes must be a named numeric matrix.", call. = FALSE)
  }
  .assert_ids(rownames(raw_genotypes), "genotype sample IDs")
  .assert_ids(colnames(raw_genotypes), "genotype marker IDs")
  if (any(!is.finite(raw_genotypes[!is.na(raw_genotypes)]))) stop("Observed genotypes must be finite.", call. = FALSE)
  train_rows <- as.integer(train_rows)
  if (!length(train_rows) || anyNA(train_rows) || anyDuplicated(train_rows) ||
      any(train_rows < 1L | train_rows > nrow(raw_genotypes))) stop("train_rows must be unique valid row indices.", call. = FALSE)
  train <- raw_genotypes[train_rows, , drop = FALSE]
  means <- colMeans(train, na.rm = TRUE)
  if (any(!is.finite(means))) stop("Every marker must have an observed training genotype.", call. = FALSE)
  af <- means / 2
  scale <- sqrt(2 * af * (1 - af))
  if (any(!is.finite(af)) || any(af <= 0 | af >= 1) || any(scale <= 0)) stop("Training allele frequencies must lie strictly between zero and one.", call. = FALSE)
  filled <- raw_genotypes
  missing <- which(is.na(filled), arr.ind = TRUE)
  if (nrow(missing)) filled[missing] <- means[missing[, 2L]]
  scaled <- sweep(sweep(filled, 2L, means, "-"), 2L, scale, "/")
  test_rows <- setdiff(seq_len(nrow(raw_genotypes)), train_rows)
  list(allele_frequency = stats::setNames(af, colnames(raw_genotypes)),
       center = stats::setNames(means, colnames(raw_genotypes)),
       scale = stats::setNames(scale, colnames(raw_genotypes)),
       all = scaled, train = scaled[train_rows, , drop = FALSE],
       test = scaled[test_rows, , drop = FALSE], train_rows = train_rows,
       test_rows = test_rows)
}

#' Subset a benchmark simulation to selected samples
#' @param simulation A validated simulation.
#' @param sample_ids Requested sample IDs in desired order.
#' @export
subset_sblrbench_simulation_samples <- function(simulation, sample_ids) {
  validate_sblrbench_simulation(simulation)
  sample_ids <- .assert_ids(as.character(sample_ids), "sample_ids")
  idx <- match(sample_ids, simulation$data$sample_ids)
  if (anyNA(idx)) stop("Requested sample IDs are missing from the simulation.", call. = FALSE)
  out <- simulation
  out$data$sample_ids <- sample_ids
  out$data$genotypes <- if (is.null(out$data$genotypes)) NULL else out$data$genotypes[idx, , drop = FALSE]
  subset_ids <- function(x) { z <- intersect(x %||% character(), sample_ids); if (length(z)) z else NULL }
  out$data$train_ids <- subset_ids(out$data$train_ids)
  out$data$test_ids <- subset_ids(out$data$test_ids)
  out$data$reference_ids <- subset_ids(out$data$reference_ids)
  for (nm in c("genetic_values", "phenotypes", "residuals")) out$truth[[nm]] <- out$truth[[nm]][idx, , drop = FALSE]
  validate_sblrbench_simulation(out)
  out
}

#' Attach aligned genetic-value predictions to a benchmark result
#' @param result A benchmark result.
#' @param genetic_value Named sample-by-trait prediction matrix.
#' @param simulation Optional simulation used for strict validation.
#' @export
add_sblrbench_predictions <- function(result, genetic_value, simulation = NULL) {
  validate_sblrbench_result(result)
  if (!is.matrix(genetic_value) || !is.numeric(genetic_value) ||
      is.null(rownames(genetic_value)) || is.null(colnames(genetic_value)) || any(!is.finite(genetic_value))) {
    stop("genetic_value must be a finite named numeric matrix.", call. = FALSE)
  }
  .assert_ids(rownames(genetic_value), "prediction sample IDs")
  .assert_ids(colnames(genetic_value), "prediction trait names")
  out <- result
  out$predictions$genetic_value <- genetic_value
  validate_sblrbench_result(out, simulation)
  out
}

#' Create paired method advantages
#' @param metrics Long metric data containing architecture, replicate, method,
#' trait, metric, value, and status.
#' @param comparisons Data frame with unique `comparison_id`, `focal_method`,
#' and `comparison_method` rows.
#' @return Paired rows with incomplete pairs retained. Positive advantages
#' consistently indicate better focal-method performance.
#' @export
paired_method_advantages <- function(metrics, comparisons) {
  required <- c("architecture", "replicate", "method", "trait", "metric", "value", "status")
  if (!is.data.frame(metrics) || length(setdiff(required, names(metrics)))) stop("metrics lacks required pairing columns.", call. = FALSE)
  needed <- c("comparison_id", "focal_method", "comparison_method")
  if (!is.data.frame(comparisons) || length(setdiff(needed, names(comparisons)))) stop("comparisons lacks required columns.", call. = FALSE)
  if (anyNA(comparisons[, needed]) || anyDuplicated(comparisons$comparison_id)) stop("comparison_id values must be unique and complete.", call. = FALSE)
  z <- metrics[, required, drop = FALSE]
  z$value[z$status != "ok"] <- NA_real_
  z$paired_metric <- z$metric
  intercept <- z$metric == "prediction_calibration_intercept"
  slope <- z$metric == "prediction_calibration_slope"
  z$value[intercept] <- abs(z$value[intercept])
  z$value[slope] <- abs(z$value[slope] - 1)
  z$paired_metric[intercept] <- "absolute_calibration_intercept_error"
  z$paired_metric[slope] <- "absolute_calibration_slope_error"
  rows <- lapply(seq_len(nrow(comparisons)), function(i) {
    cmp <- comparisons[i, ]
    focal <- z[z$method == cmp$focal_method, c("architecture", "replicate", "trait", "paired_metric", "value")]
    reference <- z[z$method == cmp$comparison_method, c("architecture", "replicate", "trait", "paired_metric", "value")]
    names(focal)[names(focal) == "value"] <- "focal_value"
    names(reference)[names(reference) == "value"] <- "comparison_value"
    out <- merge(focal, reference, by = c("architecture", "replicate", "trait", "paired_metric"), all = TRUE, sort = TRUE)
    out$comparison_id <- cmp$comparison_id
    out$focal_method <- cmp$focal_method
    out$comparison_method <- cmp$comparison_method
    out$complete_pair <- is.finite(out$focal_value) & is.finite(out$comparison_value)
    lower <- out$paired_metric %in% c("prediction_mse", "prediction_nmse", "effect_rmse",
      "absolute_calibration_intercept_error", "absolute_calibration_slope_error")
    out$advantage <- ifelse(lower, out$comparison_value - out$focal_value,
      out$focal_value - out$comparison_value)
    out
  })
  do.call(rbind, rows)
}
