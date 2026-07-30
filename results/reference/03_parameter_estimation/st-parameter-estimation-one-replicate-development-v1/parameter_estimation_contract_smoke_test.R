# Developer parameter-estimation contract smoke test (no sampler fitting).
for (f in c("estimands.R", "methods.R", "metrics.R"))
  source(file.path("studies", "03_parameter_estimation", f))
registry <- .study03_validate_registry(.study03_estimand_registry())
fake <- list(input = list(nburn = 2L, nit = 4L, nthin = 1L),
  pi_trace = matrix(c(.1, .2, .2, .3, .4, .5)),
  vbs = matrix(c(1, 2, 2, 4, 3, 5)), vgs = matrix(c(1, 1, 1, 2, 3, 4)),
  ves = matrix(c(2, 2, 2, 2, 2, 2)))
draws <- .study03_extract_draws(fake, "st_bed_bayesc", registry, 10L)
stopifnot(nrow(draws) == 24L, all(is.finite(draws$value)),
  all(draws$value[draws$estimand_id == "heritability"] >= 0),
  all(draws$value[draws$estimand_id == "heritability"] <= 1))
print(.study03_summarise_draws(draws))
