# Current `sblr` implementation map

This map describes the installed `sblr` 0.1.2 implementation at source commit
`92ff3f6e7a0b1228f9f04b693d91a36d86934b0f`. It separates confirmed public
contracts from benchmark observations and unresolved concerns. It is not a
claim that every public path has received a scientific benchmark.

## Public implementation matrix

| Trait structure | Public function | Data/operator | Public methods | annotations | residual covariance / overlap | MAF-S | multichain | benchmark status |
|---|---|---|---|---|---|---|---|---|
| single | `stblr_bed()` | individual-level packed BED | `bayesc`, `bayesr`, `bayesrc` | BayesRC | scalar residual variance | fixed or sampled where supported | native; optional compact chains | BayesC/R Studies 02--04; BayesRC Study 05 interface only |
| single | `stblr_csr()` | summary statistics + sparse-LD CSR | `sbayesc`, `sbayesr` | no | scalar residual variance | fixed or sampled where supported | native; optional compact chains | SBayesC/R Studies 02--04 |
| single | `stblr_csr_annot()` | summary statistics + sparse-LD CSR | fixed-marker, group, learned-logistic, or `sbayesrc` annotation route | yes | scalar residual variance | model-dependent | native; optional compact chains | SBayesRC Study 05 interface only |
| single | `stblr_block_eigen()` | summary statistics + reconstructed block-eigen operator | `sbayesc`, `sbayesr`, `sbayesrc` | SBayesRC | scalar residual variance | fixed or sampled where supported | native; optional compact chains | Study 06 complete five-replicate development benchmark |
| multivariate | `mtblr_bed()` | joint individual-level packed BED | `bayesc`, `bayesr`, `bayesrc` | BayesRC | full or diagonal residual covariance | fixed only for BayesR/RC; sampled MT MAF-S unsupported | native joint chains | Study 07 contract scaffold preserved; scientific benchmark paused |
| multivariate | `mtblr_csr()` | joint summary statistics + sparse-LD CSR | `sbayesc`, `sbayesr`, `sbayesrc` | SBayesRC | diagonal residual policy; `sample_overlap = "not_modeled"` | fixed only for BayesR/RC | native joint chains | Study 07 contract scaffold preserved; scientific benchmark paused |
| multivariate | `mtblr_block_eigen()` | joint summary statistics + reconstructed-dense block operator | `sbayesc`, `sbayesr`, `sbayesrc` | SBayesRC | diagonal residual policy; `sample_overlap = "not_modeled"` | fixed only for BayesR/RC | native joint chains | Study 07 execution paused pending retained low-rank MT operator |

The `S` prefix denotes a summary-statistics likelihood. It is independent of
the optional MAF-dependent effect-variance parameter `maf_effect_s`.

## Operators and filtering

BED routes decode selected marker and sample subsets directly from packed PLINK
BED files. CSR routes use a disk-backed, float32, upper-triangle sparse
correlation matrix with an implicit unit diagonal and combine it with
trait-specific `ww` cross-product diagonals. Block-eigen routes construct
contiguous, complete, non-overlapping blocks from BED provenance, reconstruct
filtered dense within-block cross-products, store their packed float32 upper
triangles, and set cross-block LD to zero.

The public block filters are `hard_truncate`, `ridge_fixed`, and `ridge_lw`.
Hard truncation has an effective 0.01 eigenvalue floor and projects summary
scores. Fixed ridge uses \(a=\eta/(1+\eta)\); `eta = 0` is the exact no-shrink
path. Ledoit--Wolf ridge estimates a block-specific shrinkage weight and
preserves each block diagonal. The runtime object is reconstructed dense block
storage, not a low-rank eigensystem.

Consequently, spectral filtering changes or regularizes the numerical LD
operator but does not reduce packed dense runtime storage. Rank reduction is
not compression in the current implementation.

## Model and output semantics

BayesC uses a binary null/active prior. BayesR uses one leading null component
and ascending positive component multipliers. BayesRC/SBayesRC use
annotation-dependent probit stick-breaking probabilities; annotations affect
component probabilities, not the LD operator.

Canonical fits record method/operator/input metadata, marker and trait order,
posterior mean effective marker effects (`bm`), posterior binary activity
probabilities (`dm`), variance/covariance summaries, convergence diagnostics,
optional retained traces, and optional compact chains. `dm` is not a component
label. For MT mixture models it is marginalized over joint pattern/component
states. MT BayesR probabilities live on the complete joint
trait-pattern-by-component simplex; marginal pattern and component
probabilities are derived and explicitly labelled.

Posterior summaries pool retained draws across validated logical chains.
Compact chain records and convergence traces have different purposes and must
not be interpreted interchangeably. Covariance fields retain their documented
units: marker-effect covariance, genetic covariance, and residual covariance
are distinct.

## Status and limitations

- Study 05: **interface validated; convergence unsupported under tested pilot
  design; five-replicate benchmark not started**.
- SBayesR: **systematic heritability overestimation observed in the
  five-replicate development benchmark; cause not yet resolved**.
- MT alignment: public CSR and block-eigen routes enforce strict alignment or
  explicit R-side statistic reordering and reject marker intersection. This is
  confirmed implementation behaviour; scientific robustness across more
  complex inputs is not benchmarked here.
- Trait and sample ordering: BED uses explicit selected rows and trait
  metadata. Summary routes preserve trait order but cannot infer unprovided
  cross-trait sample relationships. Treating external metadata as correctly
  ordered remains a user-facing input responsibility.
- Output aggregation and chain summaries: the installed schema documents
  pooled posterior means and separate compact chain records. No sblrbench MT
  benchmark has independently validated every aggregation path; this is an
  unresolved validation gap, not a confirmed bug.
- Joint states and patterns: BayesR/RC MT states are ordered as null followed by
  supplied non-null pattern order and ascending positive components. State
  count is guarded at 4096. Misinterpreting `dm` as a joint state or component
  label is a documented semantic error, not an implementation failure.
- Covariance semantics: individual-level MT BED can use full or diagonal
  residual covariance. Current MT summary routes use a diagonal residual
  policy and `sample_overlap = "not_modeled"`; they do not reconstruct unknown
  cross-trait residual covariance from marginal summaries.
- Sample overlap: unsupported in current MT summary likelihoods. This is a
  confirmed documented limitation.
- Block eigen: exact no-shrink, hard-truncated, fixed-ridge, and LW-ridge
  construction have package contract tests. Study 06 established deterministic
  runtime equivalence and a five-replicate single-trait development benchmark.
- Study 07 scope: two traits and BayesC/SBayesC only, with identical individuals,
  diagonal residual covariance, zero generating residual covariance, and no
  sample-overlap model. The verified two-trait pattern order is `0_0`, `1_0`,
  `0_1`, `1_1`. Global pattern probabilities and their chain histories are
  available, but marker-specific joint-state probabilities are not returned by
  the current BayesC schema; trait-specific `dm` values remain available.
- MT multichain summaries: `bm`, `dm`, covariance means, and probability means
  pool retained draws. Final states and covariance values come from primary
  chain 1. `*_sd`, `*_min`, and `*_max` describe variation among chain posterior
  means and are not posterior uncertainty. Formal convergence must use the
  unpooled chain histories retained by the convergence engine.
- Study 07 pause: MT-BED, MT-CSR, alignment, state, covariance, timing,
  checkpoint, and report infrastructure is preserved. The current
  `mtblr_block_eigen()` path reconstructs dense packed blocks and is explicitly
  guarded from further Study 07 execution. Primary block-eigen validation waits
  for the retained low-rank MT operator in `sblr` after Study 06 v2 review.

## Evidence

Primary evidence is the installed public function signatures and read-only
source at the installed commit, especially:

- `R/stblr-public.R`, `R/stblr-csr-annot.R`, `R/stblr-block-eigen.R`;
- `R/mtblr-bed.R`, `R/mtblr-csr.R`, `R/mtblr-block-eigen.R`;
- `R/mtblr-bayesr.R`, `R/mtblr-bayesrc.R`;
- `src/st_block_eigen.cpp`, `src/blr_block_eigen.h`;
- `docs/dev/blr_output_schema.md`,
  `docs/dev/blr_model_capability_matrix.md`,
  `docs/methods/block_eigen_operator.qmd`,
  `docs/methods/mt_bayesr_sbayesr.qmd`, and
  `docs/methods/mt_bayesrc_sbayesrc.qmd`.

Repository benchmark status comes from frozen Studies 02--06 and the Study 07
development audit. No sibling source was modified or loaded.
