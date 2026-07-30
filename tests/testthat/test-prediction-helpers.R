test_that("prediction split is deterministic, complete, and ordered", {
  ids <- paste0("id", 1:10)
  a <- make_prediction_split(ids, 0.7, 31)
  b <- make_prediction_split(ids, 0.7, 31)
  expect_identical(a, b)
  expect_length(a$train_ids, 7L)
  expect_length(intersect(a$train_ids, a$test_ids), 0L)
  expect_setequal(c(a$train_ids, a$test_ids), ids)
  expect_true(all(diff(a$train_rows) > 0))
  expect_true(all(diff(a$test_rows) > 0))
  expect_error(make_prediction_split(ids, 0), "strictly")
  expect_error(make_prediction_split(ids, 1), "strictly")
})

test_that("genotype scaling is learned only from training rows", {
  x <- matrix(c(0, 0, 1, 1, 2, 2, 0, 1, 2, 1, 0, 2), 4, 3,
              dimnames = list(paste0("s", 1:4), paste0("m", 1:3)))
  a <- training_scaled_genotypes(x, 1:3)
  changed <- x
  changed[4, ] <- c(2, 0, 1)
  b <- training_scaled_genotypes(changed, 1:3)
  expect_equal(a$allele_frequency, colMeans(x[1:3, , drop = FALSE]) / 2)
  expect_identical(a$allele_frequency, b$allele_frequency)
  expect_identical(a$train, b$train)
  expect_identical(colnames(a$train), colnames(a$test))
  expect_equal(a$all, sweep(sweep(x, 2, a$center, "-"), 2, a$scale, "/"))
})

test_that("simulation views and prediction attachment preserve contracts", {
  s <- bench_fixture()
  view <- subset_sblrbench_simulation_samples(s, rev(s$data$sample_ids[1:2]))
  expect_identical(view$data$sample_ids, rev(s$data$sample_ids[1:2]))
  expect_identical(view$truth$effects, s$truth$effects)
  expect_error(subset_sblrbench_simulation_samples(s, c("absent")), "missing")
  expect_error(subset_sblrbench_simulation_samples(s, rep(s$data$sample_ids[[1]], 2)), "unique")
  base <- new_sblrbench_result("method", effects = s$truth$effects,
                               provenance = list(note = "kept"))
  p <- s$truth$genetic_values[view$data$sample_ids, , drop = FALSE]
  out <- add_sblrbench_predictions(base, p, view)
  expect_identical(out$predictions$genetic_value, p)
  expect_identical(out$provenance, base$provenance)
})

test_that("prediction metrics are hand calculable and reject zero variance", {
  s <- bench_fixture()
  delta <- matrix(seq_len(length(s$truth$genetic_values)) / 10,
                  nrow(s$truth$genetic_values), ncol(s$truth$genetic_values))
  predicted <- s$truth$genetic_values + delta
  r <- new_sblrbench_result("method", genetic_value = predicted)
  mse <- metric_prediction_mse(s, r)
  nmse <- metric_prediction_nmse(s, r)
  phen <- metric_phenotype_prediction_correlation(s, r)
  cal <- metric_prediction_calibration(s, r)
  expect_equal(mse$value, unname(colMeans((predicted - s$truth$genetic_values)^2)))
  expect_equal(nmse$value, unname(mse$value / apply(s$truth$genetic_values, 2, var)))
  expect_equal(phen$value, vapply(seq_len(ncol(predicted)), function(j)
    cor(predicted[, j], s$truth$phenotypes[, j]), numeric(1)))
  for (j in seq_along(s$data$trait_names)) {
    expected <- coef(lm(s$truth$genetic_values[, j] ~ predicted[, j]))
    z <- cal[cal$trait == s$data$trait_names[[j]], ]
    expect_equal(z$value, unname(expected))
  }
  flat <- predicted; flat[, 1] <- 1
  failed <- metric_prediction_calibration(s, new_sblrbench_result("method", genetic_value = flat))
  expect_true(all(failed$status[failed$trait == s$data$trait_names[[1]]] == "failed"))
  zero <- s; zero$truth$genetic_values[, 1] <- 1
  expect_identical(metric_prediction_nmse(zero, r)$status[[1]], "failed")
})

test_that("paired MT advantages orient benefit positively and retain missing pairs", {
  metrics <- data.frame(architecture = "a", replicate = 1L,
    method = c("st", "mt", "st"), trait = "t",
    metric = c("prediction_correlation", "prediction_correlation", "prediction_mse"),
    value = c(.4, .6, 2), status = "ok")
  map <- data.frame(method = c("st", "mt"), representation = "BED", scope = c("ST", "MT"))
  z <- paired_mt_advantages(metrics, map)
  expect_equal(z$advantage[z$metric == "prediction_correlation"], .2)
  expect_false(z$complete_pair[z$metric == "prediction_mse"])
  expect_true(is.na(z$advantage[z$metric == "prediction_mse"]))
})
