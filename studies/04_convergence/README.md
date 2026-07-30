# Study 04: single-trait multichain convergence

Study 04 applies rank-normalized split R-hat, bulk and tail effective sample
size, Monte Carlo standard error, burn-in stability, and stable-checkpoint rules
to four matched single-trait method/architecture combinations. Four native
chains are retained for each method. The maximum development run stores 3,000
raw draws per chain with no sampler thinning.

The study evaluates convergence diagnostics, not parameter accuracy,
prediction, or method rankings. Its operational thresholds support provisional
development recommendations only.

Run the tiny contract probe first:

```r
source("studies/04_convergence/convergence_contract_smoke_test.R")
```

Then run the targets pipeline with `SBLR_BENCH_STUDY=04_convergence`. Public
qgdata inputs and the sparse-LD cache use the pinned Study 03 provenance.
