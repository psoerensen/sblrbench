# Study 01: separated-locus fine-mapping pilot

This scalar-trait pilot reuses the qgg-prepared chromosome-1 Glist, QC-retained canonical markers, sparse LD, and scaled genotype matrix. It excludes high-LD, multi-trait, annotation, and block-eigen models.

For each replicate, candidates have MAF in `[0.05, 0.5]`. A deterministic seed permutes them; a greedy scan accepts markers at least 1 Mb from all previously accepted markers. The ten selected markers are restored to chromosome order and passed as the complete `causal_rsids` pool to `sblr::mtsim()`. The exact set and `Z %*% B = G` oracle are verified.

The same phenotype, individuals, markers, and LD feed BED BayesC, BED BayesR, CSR SBayesC, and CSR SBayesR. Metrics are effect RMSE, PIP Brier score, average precision, causal rank mean/median/best/worst, top-10/20/50 recall, and exact/reference-LD proxy credible-set coverage. PIP ties use canonical marker order. Credible sets come only from public `sblr::make_credible_sets()`.

```r
Sys.setenv(SBLR_BENCH_STUDY = "01_finemapping", SBLR_BENCH_REPLICATES = "1")
targets::tar_make()
Sys.setenv(SBLR_BENCH_REPLICATES = "5"); targets::tar_make()
Sys.setenv(SBLR_BENCH_REPLICATES = "10"); targets::tar_make()
```

Changing counts creates only missing branches; shared genotype/LD targets are reused. Compact files are written below `results/local/01_finemapping/separated/`; full fits remain in `_targets/`.

The controls (`nit=500`, `nburn=250`, `nthin=1`, one chain/core) are development-only. Structural checks do not prove MCMC convergence or justify method rankings. Cumulative marginal-PIP credible sets are not SuSiE-style per-effect configurations. Scientific runs need longer MCMC, multiple chains, convergence review, and sensitivity analyses.
