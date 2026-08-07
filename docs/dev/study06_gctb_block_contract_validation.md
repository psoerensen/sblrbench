# Study 06 GCTB-compatible block contract validation

## Purpose

This pinned benchmark-level diagnostic validates the block residual contract
introduced by `sblr` 0.2.0 at source SHA
`0c89234273389e14112ba0e08ef9d11d3e1819dc`. It uses the existing 1,500-marker
Study 06 informative data and the already committed official SBayesRC v0.2.6
single-trajectory artifacts. It is neither a new qualification nor a
multichain convergence comparison.

The immutable small identities are specification
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`
and truth
`169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.
All 1,500 markers, alleles, annotations, GWAS values, train/validation people,
15 blocks, and 100 positive modes per block match the pinned official export.

## Registered comparison

| ID | New `sblr` route | Official reference | Model | Trajectory |
|---|---|---|---|---|
| S0-new | block-eigen SBayesR | D0 | matched four-component | 9,000 total; 3,000 burn; 6,000 retained |
| S1-new | block-eigen SBayesRC | D1 | matched four-component, learned alpha | 9,000 total; 3,000 burn; 6,000 retained |

Both new runs explicitly use `residual_policy = "gctb_block"`,
`block_ve_mode = "allMixVe"`, `resam_thresh = 1.1`, and
`minimum_ve_ratio = 0.7`. They retain all positive modes. The official public
seed contract cannot create demonstrably independent native chains, so these
are descriptive single-trajectory comparisons.

## Residual semantic crosswalk

| Route/evidence | Policy | `fit$ves` meaning | Heritability meaning |
|---|---|---|---|
| historical block-eigen fits | historical `global_projected` contract | global projected residual variance | historical `Vg/(Vg+projected Ve)` |
| new block-eigen fits | `gctb_block`, `allMixVe` | retained mean block residual variance | sum of block Vg divided by phenotype variance |
| BED | global individual-level residual | global individual-level residual variance | `Vg/(Vg+Ve)` |
| CSR | global operator residual | global operator residual variance | operator-contract `Vg/(Vg+Ve)` |

Historical block results retain their original semantics. New block `fit$ves`
is not BED-equivalent residual variance.

## Gates

The SBayesR gate requires PIP Pearson at least 0.99, effect and validation-g
Pearson at least 0.995, mean-block-Ve difference after proved phenotype-scale
conversion at most 0.005, and summary-h2 difference at most 0.05. The SBayesRC
gate differs only in allowing PIP Pearson at least 0.95. Active counts, exact
occupancy, alpha, and `sigmaSqAlpha` are descriptive rather than gates.

Implementation: [gctb-block-contract-validation.R](../../studies/06_annotation_models/gctb-block-contract-validation.R)
and [runner](../../scripts/run_study06_gctb_block_contract_validation.R).
