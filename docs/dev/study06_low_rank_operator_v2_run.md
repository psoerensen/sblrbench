# Historical pre-renumbering Study 06 v2 retained-low-rank validation

> **Retired numbering:** this completed evidence is now integrated into Study
> 05. Commands and paths below document the historical run only.

> **Retired internal runbook.** This file preserves execution provenance for
> the immutable capsule. The v2 orchestration and launchers were removed after
> migration to the shared framework. See [study06_migration.md](study06_migration.md)
> for the active specification, CLI, checkpoint policy, and report contract.

## Finalization control gate (2026-08-02)

Final checkpoint validation found that all five `sparse_mixture / full_csr`
benchmark fits retained 1,000 iterations. The frozen v2 convergence evidence
requires 2,000 retained iterations after 250 burn-in iterations. This was an
execution-dispatch defect: `full_csr` incorrectly inherited the Study 04
baseline controls instead of its Study 06 v2 convergence recommendation.

The five checkpoints, all other checkpoints, aggregates, and generated report
were preserved. The rejected benchmark capsule and report remain in ignored
local quarantine and are not valid frozen reference evidence. Benchmark
dispatch now reads the v2 recommendation for `full_csr`, and capsule validation
requires exact agreement with every non-BED recommendation.

Only the five rejected `sparse_mixture / full_csr` coordinates were rerun,
using their original deterministic fit and chain seeds with `nburn = 250`,
`nit = 2000`, `nthin = 1`, `nchains = 4`, and `ncores = 4`. All 55 unaffected
checkpoint SHA-256 hashes, sizes, and modification times were unchanged. The
corrected 60-fit grid passed final validation and was reaggregated into a new
benchmark capsule; quarantined evidence was excluded.

## Status and provenance

Study 06 v2 is a new validation study. The completed Study 06 v1 capsules and
their reconstructed-dense results are immutable historical evidence. V2 uses
the scalar retained low-rank operator from `sblr` 0.2.0 at source revision
`bd8e2c8148a0d9540dc20716455706beeb16fa86` (Optimize retained low-rank
scalar sampling).

The ordinary shared user library still contains `sblr` 0.1.2 because a
long-lived unrelated RStudio process holds that DLL. The optimized package is
built and installed normally from the clean local sibling into the isolated
Study 06 library `results/local/study06_low_rank_operator_v2/rlib`. Each fresh
runner prepends that library and verifies version, `RemoteSha`, API, and
namespace path before work begins. `pkgload::load_all()` is not used for
`sblr`, and the sibling repository remains read-only.

The pre-optimization v2 run is preserved at
`results/local/study06_low_rank_operator_v2_pre_optimization_bd8e2c8` and is
ineligible for combination with optimized checkpoints. The optimized backend
uses the unchanged statistical model and retained-rank policy, explicitly
records `low_rank_residual_rebuild_every = 100`, and returns rebuild-count and
maximum-drift diagnostics.

## Historical reconstructed-dense Study 06 v1

V1 used `sblr` 0.1.2 at
`92ff3f6e7a0b1228f9f04b693d91a36d86934b0f`. Its convergence capsule contains
8 successful fits and its five-replicate capsule contains all 60 successful
fits. V1 evaluated reconstructed-dense unfiltered, hard-thresholded, and
fixed-ridge routes. Those policies are not scientific configurations in v2,
and v1 values are never labelled retained low rank.

## V1-to-v2 design crosswalk

| Field | V1 | V2 disposition |
|---|---|---|
| Architectures | sparse homogeneous; sparse mixture | identical |
| Replicates | five paired replicates | identical |
| Samples and markers | 5,000 samples; canonical 37,991 markers | identical |
| Split | deterministic 70/30, seed 3101 | identical |
| Truth | 50 causal markers; realized h2 0.30 | identical |
| Simulation seeds | base 61000 with frozen strides | identical |
| BED/full CSR/block CSR | three reference configurations | retained |
| Block design | 38 contiguous blocks, target 1,000 markers | identical |
| Dense block eigen | unfiltered, hard, fixed ridge | historical only; removed |
| Retained low rank | unavailable | near-full, 0.999, and 0.995 |
| Operator contract | reconstructed packed dense | `block_low_rank_v1` |
| Residual state | marker-space dense reconstruction | retained Q/w residual state |
| Residual variance | historical reconstructed operator | global projected SSE contract |
| Primary grid | 2 x 5 x 6 = 60 | 2 x 5 x 6 = 60 |

## V2 execution policy

Every low-rank fit must pass `representation = "low_rank"` and an explicit
`eigen_prop`. Canonical model specifications containing `dense_reconstructed`,
`eigen_filter`, `eigen_tau`, or `eigen_eta` fail before dispatch. A single
unfiltered reconstructed-dense configuration is permitted only in the
post-optimization two-architecture pilot as a labelled implementation
comparator; it cannot enter convergence or the primary 60-fit grid. Direct
`crossprod(Q)` and `crossprod(Q, w)` calculations remain the deterministic
dense eigenspace oracle.

The six primary configurations are packed BED, full CSR, block CSR, retained
near-full positive rank, retained 0.999, and retained 0.995. The exact v1 block
boundaries are reused. The canonical 0.995 route is the primary future
low-rank configuration.

## Projected residual variance

The retained SBayesRC-style eigenspace likelihood is expressed in `sblr`
cross-product units. It uses the global projected contract

`SSE = yy - sum_b crossprod(w_b) + sum_b crossprod(r_b)`.

It does not reproduce GCTB's block-specific residual-variance procedure.
Deterministic and fitted-state validation must establish the equivalent
quadratic form and nonnegative finite projected SSE before scientific claims.
