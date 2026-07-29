# sblrbench

`sblrbench` provides small, explicit contracts for simulation-based evaluation of scalable Bayesian linear regression methods. Statistical models, samplers, data preparation, output semantics, convergence, structural diagnostics, and credible sets remain the responsibility of `sblr`; larger simulations, strict alignment, oracle checks, adapters, truth-aware metrics, and provenance belong here. This package calls only exported functions from the installed `sblr` package.

Install with `remotes::install_github("psoerensen/sblrbench")` and attach with `library(sblrbench)`. For local development use `devtools::load_all()`. To install the clone and run the default contract study:

```r
devtools::install(upgrade = "never")
targets::tar_make()
```

Version 0.1 offers `sblrbench_simulation`, strict marker/sample/trait alignment, an oracle `Z %*% B` check, lightweight method and result lists, native public-API adapters, four metrics in stable long format, and compact JSON manifests.

```r
simulation <- as_sblrbench_simulation(sblr::mtsim(...))
oracle <- check_oracle_genetic_values(simulation)

method <- new_sblrbench_method(
  id = "example", label = "Example", capabilities = "posterior_effects",
  fit = function(...) list(...), extract = function(x) x
)
result <- run_sblrbench_method(method)
metrics <- evaluate_metrics(simulation, result,
                            metrics = c("effect_rmse", "pip_brier"))
```

Alignment is identity-based, restores canonical order, and rejects duplicates, missing IDs, unnamed axes, and extras by default. Prediction studies must pass the oracle genotype-scale check before scoring. Native adapters retain the unmodified fit optionally and map only compatible common fields.

Add a future method by constructing an ordinary `sblrbench_method`; no registry is required. Start a future study under `studies/<name>/`, keeping study-specific code there until it becomes a stable shared contract. External adapters, full fine-mapping/prediction benchmarks, workflow engines, and method-specific semantic normalization are deliberately deferred.

Study 01 Task 1 validates the qgg-backed genotype path: `qgg::gprep()` creates or loads a Glist, `qgg::gfilter()` performs QC, chromosome order defines the canonical markers, `sblr::make_sparse_ld()` prepares LD, and `qgg::getG(impute = TRUE, scale = TRUE)` produces the matrix passed unchanged to `sblr::mtsim(standardize_W = FALSE)`. `sblrbench` then proves `Z %*% B = G` with its strict oracle. qgg owns genotype reading, imputation, and scaling; no model fitting occurs in this task. Select it with `SBLR_BENCH_STUDY=01_finemapping` before `targets::tar_make()`.
