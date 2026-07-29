# Study 01: fine-mapping genotype setup

Task 1 is development infrastructure, not a scientific benchmark:

```text
PLINK BED/BIM/FAM
    ↓
qgg::gprep()
    ↓
qgg::gfilter()
    ↓
canonical chromosome marker order
    ↓
sblr::make_sparse_ld()
    ↓
qgg::getG(impute = TRUE, scale = TRUE)
    ↓
sblr::mtsim(standardize_W = FALSE)
    ↓
sblrbench oracle validation
```

qgg owns Glist preparation, marker QC, genotype extraction, missing-genotype imputation, and scaling. `sblrbench` does not reimplement these operations. sblr owns sparse-LD construction and simulation. No model is fitted in Task 1. A passing oracle proves the internal identity `Z %*% B = G` after strict alignment; it does not establish model correctness.

By default the pipeline downloads the five public `human.*` qgg example files only when absent, under ignored `results/local/01_finemapping/data/`. Override that directory with `SBLR_BENCH_DATA_DIR`, or set `SBLR_BENCH_GLIST` to an existing Glist RDS. The original RDS is never overwritten. Generated Glist caches, sparse LD, matrices, and simulations remain local and ignored.

From the repository root:

```powershell
$env:SBLR_BENCH_STUDY = "01_finemapping"
Rscript -e "targets::tar_make()"
```

Future work remains the documented one-trait pilot with separated and high-LD causal architectures, 10 replicates, ST-BED BayesC/BayesR and ST-CSR SBayesC/SBayesR, initial PIP Brier scoring, and later ranking and exported-sblr credible-set metrics. None of that fitting or scoring is implemented here. No study work may modify `sblr`.

In the development environment recorded in `docs/sblr_inventory.md`, installed `sblr` 0.1.0 could not complete the full example matrix's one-trait `mtsim(h2 = 0.2)` residual-covariance step. The pipeline stops with the native error and writes no oracle success summary. Run the command above with a compatible installed public `sblr`; do not install or patch the sibling source as a workaround.
