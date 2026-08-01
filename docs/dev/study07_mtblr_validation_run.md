# Study 07: multivariate implementation validation

## Paused execution status

Study 07 is paused pending implementation and validation of the retained
low-rank MT operator in `sblr`.

The current `mtblr_block_eigen()` backend uses the historical
reconstructed-dense block representation and must not be used for future
primary low-rank analyses. MT-BED, MT-CSR, alignment, state-mapping,
covariance, runtime-scaling, checkpoint, and reporting infrastructure remain
valid and are preserved.

This is a pause, not abandonment. Deterministic contracts and bounded MT-BED,
full MT-CSR, and runtime-matched block MT-CSR interface checks remain
available. The primary MT block-eigen benchmark will resume only after Study
06 v2 has been reviewed and the retained low-rank MT implementation has been
implemented and validated in `sblr`. Existing local block-eigen timing evidence
describes the reconstructed-dense backend and must not be presented as
retained-low-rank evidence.

The unattended runner conservatively pauses every long phase (`runtime`,
`convergence`, `benchmark`, `stress`, `aggregate`, and `all`). Matching guards
at the start of the long target bodies prevent direct `tar_make()` calls from
dispatching fits. The fit wrapper independently rejects the historical MT
block-eigen implementation ID. Audit, validation-only, sampler-free contracts,
and bounded MT-BED/MT-CSR contract work remain available.

## Provenance and boundaries

Development began on `master` at commit
`a4bd3fdaedbe6e36400c97318299f1965fba72bf` with a clean working tree.
The installed packages are `sblr` 0.1.2 from source commit
`92ff3f6e7a0b1228f9f04b693d91a36d86934b0f` and `sblrbench`
0.0.0.9000. Both cached copies of qgdata commit
`6cca5819e711d326cfb2614d7e9d9f34942612cd` passed pinned size and MD5
validation. No network or package installation is permitted.

The sibling `../sblr` tree is used only for read-only contract audit. It is
never loaded or compiled. Consulted files are recorded in `config.R` and in
each promoted capsule's `interface_audit_sources.csv`.

## Narrow scientific scope

The first benchmark uses two traits, BayesC/SBayesC, one common genotype/LD
resource, identical individuals, diagonal residual covariance, and zero
generating residual covariance. It excludes BayesR, annotations, sampled MAF-S,
partial sample overlap, correlated residuals, and more than two traits.

## Public interface map

| Function | Method | Inputs | Residual policy | Summary overlap policy |
|---|---|---|---|---|
| `mtblr_bed()` | `bayesc` | phenotype matrix `n x 2`, one BED-backed Glist, selected rows/markers | `full` or `diagonal`; Study 07 explicitly uses `diagonal` | not applicable |
| `mtblr_csr()` | `sbayesc` | two-trait sufficient statistics plus shared or trait-specific CSR LD | diagonal | `sample_overlap = "not_modeled"`; marginal `yy` only |
| `mtblr_block_eigen()` | `sbayesc` | same-BED sufficient statistics, Glist and contiguous block starts | diagonal | `sample_overlap = "not_modeled"`; marginal `yy` only |

Phenotypes are sample-by-trait. Marker effects and states are marker-by-trait.
`wy` and `ww` are trait lists in canonical marker order; `yy` and `n` have one
entry per trait. Summary statistics and LD use standardized genotypes.
The historical reconstructed-dense block-eigen Study 07 route used
`ridge_fixed` with `eigen_eta = 0` and the Study 06 1,000-marker contiguous
block convention. That execution route is now guarded and is not the future
retained low-rank benchmark.

`nit` is the retained post-warm-up draw count, `nburn` is additional warm-up,
and `nthin` affects returned posterior summaries but not convergence traces.
Each native multichain call fits one complete joint model per chain. Explicit
signed 32-bit chain seeds preserve supplied order.

## Joint-state contract

The exact two-trait BayesC model matrix is:

| Internal ID | Label | Trait 1 | Trait 2 | Internal key |
|---:|---|---:|---:|---|
| 0 | neither | 0 | 0 | `0_0` |
| 1 | trait 1 only | 1 | 0 | `1_0` |
| 2 | trait 2 only | 0 | 1 | `0_1` |
| 3 | both traits | 1 | 1 | `1_1` |

This order follows the public default `expand.grid(rep(list(0:1), 2))` and is
also passed explicitly. `pimodels` and the returned `pi_mean` are global
probabilities over these four patterns. Trait-1 PIP is the sum of states 1 and
3; trait-2 PIP is the sum of states 2 and 3. Trait permutation swaps states 1
and 2 and leaves states 0 and 3 unchanged.

For BayesC, `dm` is a marker-by-trait marginal inclusion probability. The
installed schema does not return marker-specific four-state posterior
probabilities. Such metrics are marked unavailable; they are not reconstructed
from unrelated objects.

## Output and chain semantics

- `bm` and `dm`: pooled retained-draw posterior marker means and marginal PIPs.
- `cov_g_mean`, `cov_e_mean`, `cov_b_mean`: pooled retained-draw covariance
  means.
- `b`, `d`, final covariance matrices and `pi_final`: primary-chain final state.
- `vgs`, `ves`, `vbs`, `vld`, `vle`: iterationwise trace summaries.
- extended convergence covariance and probability groups: actual unpooled,
  post-burn per-chain histories used for R-hat, ESS, and MCSE.
- `bm_chain_mean_sd/min/max` and `dm_chain_mean_sd/min/max`: variation among
  per-chain posterior means, not posterior SD or credible intervals.

Compact chains omit shared BED data, phenotypes and genetic values. Prediction
is reconstructed explicitly from the frozen genotype scale and `bm`.

The installed MT-BED formatter exposes `nchains` at top level. MT-CSR and
MT-block-eigen retain the same identifiable chain list and record the count as
`input$nchains`, but do not populate top-level `nchains`. Study 07 validates the
top-level value when present and otherwise uses the explicit input metadata.

`make_summary_stats()` returns the per-chromosome `cls` list named (`chr1` in
this benchmark) while its aligned `af` list is unnamed. MT-block-eigen compares
the named `lengths()` vectors exactly. Study 07 assigns the existing `cls` list
names to `af`; values, marker order, and allele frequencies are unchanged.

## Deterministic gates

Permanent tests cover state round trips and probability sums, trait/marker/
sample permutations, effect-to-genetic-value reconstruction, covariance and
genetic-correlation identities, positive-semidefinite checks, deterministic
seeds, summary-statistic alignment, runtime-matched CSR/block-eigen actions,
chain identity, runtime-table selection, capsule checksums, and safe optional
manifest fields.

The runtime-matched block CSR and unfiltered block eigen must use the same
reconstructed runtime blocks. Matrix, diagonal, matrix-vector, matrix-matrix,
quadratic-form, summary-score, and true-effect actions must pass the frozen
Study 06 tolerances before any scientific sampler phase.

The installed MT-CSR alignment code recognizes BED-by-construction orientation
when the Glist descriptor uses its canonical `make_sparse_ld` source token and
the marker ID, chromosome/file, and BED-column provenance are identical. The
runtime-matched CSR is reconstructed from the same selected BED Glist, retains
those exact fields, and records the additional reference ID
`study07_runtime_matched_block_reconstruction`.

## Runtime gate

Nested 1,000, 2,000, 4,000 and conditionally 8,000-marker subsets are timed.
The largest common feasible count is frozen using these development limits:

- projected four-chain fit no longer than two hours;
- projected 30-fit grid no longer than 36 hours;
- estimated memory no more than 70% of detected physical RAM, or a conservative
  16 GiB configured limit when detection is unavailable;
- one fit object no larger than 5 GiB.

These are operational limits, not claims about universal implementation scale.

## Known limitations

This design does not validate arbitrary overlapping GWAS samples, correlated
residuals, trait-specific LD panels, missing phenotypes, BayesR/SBayesR,
annotation models, sampled MAF-S, or more than two traits. Five replicates are
descriptive development evidence rather than a universal method ranking or
convergence result.
