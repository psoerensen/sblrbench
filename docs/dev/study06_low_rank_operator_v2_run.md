# Study 06 v2: retained low-rank LD operator validation

## Status and provenance

Study 06 v2 is a new validation study. The completed Study 06 v1 capsules and
their reconstructed-dense results are immutable historical evidence. V2 uses
the scalar retained low-rank operator from `sblr` 0.2.0 at source revision
`96487b3194fc1f8c6789060da5f2e2a0eea89974` (Correct BayesR prior
variance calibration).

The ordinary user library still contains `sblr` 0.1.2. V2 does not install or
update it. Each v2 R session loads a compiled, generated `git archive` snapshot
of the pinned clean sibling source from the v2 local preflight directory and
verifies its version, API, and source path before work begins. The sibling
repository remains read-only.

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
`eigen_prop`. Model specifications containing `dense_reconstructed`,
`eigen_filter`, `eigen_tau`, or `eigen_eta` fail before dispatch. Direct
`crossprod(Q)` and `crossprod(Q, w)` calculations are the only dense
eigenspace oracle. The reconstructed-dense sampler is never called.

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
