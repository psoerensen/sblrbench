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

Study 01 implements a separated-locus fine-mapping pilot. Ten exact causal markers are selected by a seeded randomized-greedy scan after MAF filtering, retained in chromosome order, and separated by at least 1 Mb. The four scalar-trait methods are BED BayesC/BayesR and CSR SBayesC/SBayesR. Set `SBLR_BENCH_STUDY=01_finemapping` and `SBLR_BENCH_REPLICATES` to `1`, `5`, or `10` before `targets::tar_make()`; genotype and sparse-LD targets are cached across stages.

Outputs cover effect RMSE, PIP Brier score, average precision, causal ranks/top-K recall, and exact/LD-proxy credible-set coverage. The configured 500-iteration, 250-burn-in, one-chain controls are development settings only. Structural success does not prove convergence or support rankings. Cumulative marginal-PIP credible sets are not per-effect causal configurations; scientific runs require longer chains and sensitivity analyses.

## Website

Reusable benchmark reports will be published at <https://psoerensen.github.io/sblrbench/>. The repository owner may need to select **GitHub Actions** as the Pages source in repository settings.

## Running a study

```r
Sys.setenv(SBLR_BENCH_STUDY = "01_finemapping")
targets::tar_make()
```

## Published benchmark capsules

Reviewed, compact, versioned benchmark snapshots live under `results/reference/`. Working outputs and full local study artifacts remain ignored under `results/local/`.

## Reproducing reports

```bash
quarto render
```

Report rendering reads only tracked reference capsules and does not run the targets pipeline.

Study 01 uses public simulated PLINK data from a pinned `psoerensen/qgdata` revision. Download and validate the five example files with:

```r
download_sblrbench_example_data()
```

The improved `separated-development-v1.1` capsule contains the same numerical benchmark results as v1 with a real worked sampler example, a separate developer smoke test, pinned data checksums, and clearer full-workflow reproduction guidance.
