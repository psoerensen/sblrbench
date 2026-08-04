# Study 06 LD-operator development run

> **Retired internal runbook.** This records the historical execution that
> produced the frozen evidence. The launchers and per-study targets graph named
> below were removed by the shared-framework migration. New validation uses
> `scripts/run_benchmark.R --study 06_ld_operator`; see
> [study06_migration.md](study06_migration.md). Commands below are provenance,
> not active instructions.

## Scope and initial provenance

Study 06 is a package-specific development benchmark of the installed `sblr`
implementation. It compares one individual-level BED fit and five
summary-statistics/operator configurations for each matched BayesC or BayesR
simulation. `S` in SBayesC and SBayesR means summary statistics; it is not the
MAF-dependent `S` parameter. Annotation-informed models and MTBLR are outside
the benchmark grid.

The clean-tree gate passed on branch `master` at
`0d5b7d854e655c88aac69cef59279be513f4b37d`. The initial
`git status --short --untracked-files=all` was empty. Installed versions were
`sblr` 0.1.2 (source commit
`92ff3f6e7a0b1228f9f04b693d91a36d86934b0f`) and `sblrbench`
0.0.0.9000. Study 05 was tracked and clean. No active R benchmark process was
present.

The cached qgdata source commit is
`6cca5819e711d326cfb2614d7e9d9f34942612cd`. Both existing local cache copies
matched these pinned values:

| file | bytes | MD5 |
|---|---:|---|
| `human.bed` | 62,500,003 | `e89bea9a6cedd9eeef3fd0a5c807db81` |
| `human.bim` | 1,882,359 | `0105119b04c67b7ac7f66cc5e6680963` |
| `human.fam` | 117,786 | `3c5db3d9eb7f3fc893c75f6f2b89836d` |
| `human.pheno` | 92,786 | `6a9e7cb1162e43999c170a363863176d` |
| `human.covar` | 641,513 | `d06002aa2b1b79bdc4c0e92c21f27ae5` |

## Verified public interfaces

| Configuration | Public entry point | method |
|---|---|---|
| BED BayesC/BayesR | `sblr::stblr_bed()` | `bayesc`, `bayesr` |
| full or block CSR | `sblr::stblr_csr()` | `sbayesc`, `sbayesr` |
| block-eigen | `sblr::stblr_block_eigen()` | `sbayesc`, `sbayesr` |

`stblr_bed()` requires a phenotype and BED-backed `Glist`.
`stblr_csr()` requires `stats` plus either a CSR `ld_prefix` or matching sparse
LD provenance in `Glist`. `stblr_block_eigen()` requires `stats`, the
BED-backed reference `Glist`, and one-based contiguous `block_start` values.
Public block starts must begin at one; the native builder uses zero-based
starts.

All summary statistics come from `sblr::make_summary_stats(scale = TRUE)`.
For each trait, `wy` is \(X^\top y\), `ww` is the marker diagonal
\(x_j^\top x_j\), `yy` is \(y^\top y\), and `n` is the analysis sample count.
Markers are standardized with the training-reference allele frequency as
\((g-2p)/\sqrt{2p(1-p)}\). Missing BED genotypes decode to zero on that scale.
The sparse-LD CSR stores the standardized correlation upper triangle with an
implicit unit diagonal; the sampler combines it with `ww`. Marker order must
be identical across `stats`, CSR metadata, BED columns, allele frequencies,
and block definitions.

The block-eigen builder reads the same selected BED rows and columns. In each
block it forms the double-precision cross-product \(A=Z^\top Z\), obtains
\(C=D^{-1/2}AD^{-1/2}\), applies the requested filter, reconstructs a dense
runtime cross-product, and stores its upper triangle as float32. Cross-block
entries are zero.

### Filtering contract

- `hard_truncate`: retains eigenvalues of `C` at least
  `max(eigen_tau, 0.01)` and projects `wy` into the retained eigenspace.
  Consequently, `eigen_tau = 0` is **not** unfiltered.
- `ridge_fixed`: uses shrinkage weight
  \(a=\eta/(1+\eta)\), reconstructing
  \((1-a)A+a\,\mathrm{diag}(A)\). Therefore
  `eigen_filter = "ridge_fixed", eigen_eta = 0` is the verified exact
  no-shrink path.
- `ridge_lw`: computes a block-specific Ledoit--Wolf-style weight, shrinks
  off-diagonal entries toward zero, preserves the cross-product diagonal, and
  does not project `wy`.

For `ridge_lw`, the installed source computes
`a = min(bbar, d2) / d2`, clipped to `[0, 1]`, and reconstructs
`(1-a)A + a diag(diag(A))`. Thus `a = 0` means no shrinkage and `a = 1`
means complete removal of off-diagonal entries. The target is the diagonal
matrix on the standardized-genotype cross-product scale, not identity. The
training-row count is the sample size supplied to the formula. The qgdata
pilot produced weights near one, so Ledoit--Wolf is retained as deterministic
sensitivity evidence and the main sixth configuration uses a frozen 1%
fixed-ridge weight (`eigen_eta = 0.01 / 0.99`).

Negative eigenvalues are retained by the exact no-shrink ridge path because it
reconstructs `A` directly. Hard truncation drops eigenvalues below the
effective threshold. The runtime operator is the packed reconstructed matrix,
not a retained low-rank factorization, so runtime storage is
\(O(\sum_b m_b(m_b+1)/2)\).

Spectral filtering in the current implementation changes or regularizes the
numerical LD operator but does not reduce the packed dense runtime storage
requirement. Study 06 does not describe this filtering as compression.

Three diagonals are kept distinct:

1. `stats$ww`, the summary-statistic/reference diagonal;
2. the diagonal of the unfiltered block cross-product;
3. the filtered packed runtime diagonal.

The deterministic gate constructs runtime-matched block CSR from the exact
float32 packed no-shrink operator. It does not compare matrices with different
cleanup or reconstruction histories.

## Design revision after the first deterministic pilot

The initial implementation conservatively required a public hard threshold to
reduce rank. The qgdata blocks were full rank at thresholds 0.02, 0.05, and
0.10, so that rule stopped the first sampler-free run. That evidence is
preserved in the local log. The rule was superseded because Study 06 tests
operator and downstream BLR agreement, not compression.

The frozen real-data public hard route is now `hard_truncate` with
`eigen_tau = 0.10`. A full-rank result is recorded as `effective_no_op`, not a
failure. Separate deterministic matrices with controlled spectra verify that
the hard route removes eigenvalues when values below the effective threshold
exist. The primary block design remains 1,000 markers; the global fraction of
squared LD retained is descriptive and is not itself a validity gate.

The scientific interpretation order is:

1. runtime-matched block CSR versus zero-ridge block eigen;
2. full CSR versus runtime-matched block CSR;
3. hard and fixed-ridge routes versus zero-ridge block eigen;
4. summary-statistics configurations versus BED.

### MCMC and outputs

For all three entry points, `nit` is retained post-warm-up draws, `nburn` is
additional warm-up, and `nthin = 1` means no thinning. Native multichain
execution uses four explicit distinct `chain_seeds`, `nchains = 4`,
`ncores = 4`, and `keep_chains = TRUE`. Extended trace retention is used for
applicable low-dimensional quantities.

Canonical fits include marker posterior means (`bm`), marker non-null/activity
probabilities (`dm`), variance traces such as `vgs`, `ves`, and `vbs`,
mixture summaries such as `pis` for BayesR, and chain/convergence metadata.
Summary-statistics fits may additionally expose `vld`, `vle`, `log_cpo`, and
`mean_log_cpo`. Similar names are not treated as interchangeable: every
extracted field is tagged with its package path and availability. Operator
metadata includes filter controls and per-block diagnostics (start, size,
retained rank, minimum eigenvalue, and shrinkage weight); the deterministic
gate independently records reconstructed matrices, diagonals, ranks, spectral
mass, products, and quadratic forms.

## Read-only source audit

The installed package is the only package used for fits. The sibling checkout
was inspected read-only with `git show` at the installed source commit and was
not loaded or modified. Consulted files are:

- `R/stblr-public.R`
- `R/stblr-block-eigen.R`
- `R/stblr-csr-annot.R`
- `R/mtblr-bed.R`
- `R/mtblr-csr.R`
- `R/mtblr-block-eigen.R`
- `R/mtblr-bayesr.R`
- `R/mtblr-bayesrc.R`
- `R/interface_mtblr.R`
- `R/sparse_ld_bed_helper.R`
- `src/st_bed_decode.h`
- `src/st_block_eigen.cpp`
- `src/blr_block_eigen.h`
- `tests/testthat/helper-blr-fixtures.R`
- `tests/testthat/test-stblr-block-eigen.R`
- `tests/testthat/test-blr-operator-reductions.R`
- `docs/methods/block_eigen_operator.qmd`
- `docs/dev/blr_output_schema.md`
- `docs/dev/blr_model_capability_matrix.md`
- `docs/methods/mt_bayesr_sbayesr.qmd`
- `docs/methods/mt_bayesrc_sbayesrc.qmd`
- `examples/workflows/operator_comparison_workflow.R`

## Execution

Sampler-free validation:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_study06_ld_operator.ps1 --validate-only
```

Detached resumable execution:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_study06_ld_operator.ps1 --phase all --resume
```

Foreground fallback:

```powershell
Rscript scripts/run_study06_ld_operator.R --phase all --resume
```

Progress is written atomically to
`results/local/study06_ld_operator/status.csv`. Persistent phase logs are in
the same directory. No renderer loads qgdata, constructs an operator, or calls
a sampler.
