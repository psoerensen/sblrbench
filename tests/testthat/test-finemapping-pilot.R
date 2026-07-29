pilot_helper <- testthat::test_path("..", "..", "studies", "01_finemapping", "pilot.R")
.pilot_available <- file.exists(pilot_helper)
if (.pilot_available) source(pilot_helper, local = TRUE)

test_that("separated causal selection is deterministic and respects distance", {
  skip_if_not(.pilot_available, "repository-only pilot helpers are excluded from package builds")
  ids <- paste0("m", 1:8)
  gl <- list(rsids = list(ids), pos = list(seq(0, by = 10, length.out = 8)), maf = list(rep(.2, 8)))
  a <- select_separated_causal_markers(gl, 1, ids, 3, 20, 42, .05, .5)
  b <- select_separated_causal_markers(gl, 1, ids, 3, 20, 42, .05, .5)
  expect_identical(a, b); expect_length(a$marker_ids, 3); expect_true(all(a$marker_ids %in% ids))
  expect_gte(a$pairwise_min_distance_bp, 20)
  expect_error(select_separated_causal_markers(gl, 1, ids, 5, 20, 42), "Insufficient separated")
})

test_that("average precision is hand calculated with canonical ties", {
  s <- bench_fixture(); s$data$trait_names <- "t1"; s$truth$effects <- s$truth$effects[, "t1", drop = FALSE]
  s$truth$genetic_values <- s$truth$genetic_values[, "t1", drop = FALSE]; s$truth$phenotypes <- s$truth$phenotypes[, "t1", drop = FALSE]
  s$truth$residuals <- s$truth$residuals[, "t1", drop = FALSE]; s$truth$causal <- list(shared = c("m1", "m3"), specific = list(t1 = character()), all = c("m1", "m3"))
  p <- matrix(c(.8, .8, .2), 3, 1, dimnames = list(s$data$marker_ids, "t1"))
  expect_equal(metric_average_precision(s, new_sblrbench_result("x", pip = p))$value, (1 + 2/3)/2)
})

test_that("causal rank summaries and recalls are exact", {
  s <- bench_fixture(); p <- matrix(c(.2, .9, .1, .8, .1, .2), 3, 2, dimnames = list(s$data$marker_ids, s$data$trait_names))
  z <- metric_causal_ranks(s, new_sblrbench_result("x", pip = p)); t1 <- z[z$trait == "t1", ]
  expect_equal(t1$value[t1$metric == "causal_rank_mean"], 2)
  expect_equal(t1$value[t1$metric == "causal_rank_median"], 2)
  expect_equal(t1$value[t1$metric == "causal_rank_best"], 2)
  expect_equal(t1$value[t1$metric == "causal_rank_worst"], 2)
  expect_equal(t1$value[t1$metric == "causal_top_1_recall"], 0)
})

test_that("credible-set scoring distinguishes exact, proxy, and uncovered", {
  skip_if_not(.pilot_available, "repository-only pilot helpers are excluded from package builds")
  s <- bench_fixture(); s$data$trait_names <- "t1"; s$truth$effects <- s$truth$effects[,1,drop=FALSE]
  s$truth$genetic_values <- s$truth$genetic_values[,1,drop=FALSE]; s$truth$phenotypes <- s$truth$phenotypes[,1,drop=FALSE]; s$truth$residuals <- s$truth$residuals[,1,drop=FALSE]
  s$truth$causal <- list(shared = c("m1", "m3"), specific = list(t1=character()), all=c("m1","m3"))
  native <- list(sets = list(L1 = list(S1 = c("m1")), L2 = list(S2 = c("m2"))),
    loci = data.frame(locus=c("L1","L2"), start=c(0,90), end=c(10,110)), summary=data.frame())
  LD <- diag(3); dimnames(LD) <- list(paste0("m",1:3),paste0("m",1:3)); LD["m3","m2"] <- LD["m2","m3"] <- .8
  z <- evaluate_credible_sets(list(status="ok",reason="",native=native), s, c(m1=5,m2=100,m3=105), LD, "x", .5, 20)
  expect_true(z$covered[z$set_id == "L1:S1" & z$coverage_type == "exact"])
  expect_false(z$covered[z$set_id == "L2:S2" & z$coverage_type == "exact"])
  expect_true(z$covered[z$set_id == "L2:S2" & z$coverage_type == "ld_proxy"])
  expect_true(all(z$set_size == 1L))
})
