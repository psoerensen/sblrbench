# Study 06 large information-scale feasibility design

## Status

This is a single-replicate feasibility experiment, not a formal qualification or
a final benchmark. It preserves the existing Study 06 small information-scale
stress test (`n_train = 1,400`, `m = 1,500`) and registers a separate experiment
with all 5,000 canonical individuals and 37,991 chromosome-1 markers.

The ordinary sampler is the object of study: one allocation sweep and one alpha
hierarchy update per cycle. H20 schedules, repeated allocation sweeps,
tempering, partial exchange, official SBayesRC, and sampler changes are excluded.

The original execution phase was technically blocked (`LARGE-F6`). The frozen
continuation uses the validated package source SHA
`0c89234273389e14112ba0e08ef9d11d3e1819dc` and the explicitly registered
GCTB-compatible block residual contract. It does not change any scientific
input, seed, prior, chain length, or model setting.

## Scientific question

Does ordinary learned-alpha BayesRC/SBayesRC achieve multichain convergence and
recover the known annotation architecture when marker, allocation, and GWAS
information are substantially increased?

The six registered fits isolate general BayesR sampling, fixed heterogeneous
marker priors, joint alpha/allocation feedback, BED versus block-eigen routes,
and SNP utility versus latent-architecture convergence.

| ID | Route | Model | Alpha |
|---|---|---|---|
| E0 | individual-level BED | BayesR | absent |
| B0 | block eigen | SBayesR | absent |
| E2 | individual-level BED | BayesRC | fixed to truth |
| B2 | block eigen | SBayesRC | fixed to truth |
| E1 | individual-level BED | BayesRC | learned |
| B1 | block eigen | SBayesRC | learned |

Execution order is E0/B0, E2/B2, then E1/B1. Each scientific fit has four
chains, 12,000 iterations, 3,000 burn-in iterations, 9,000 retained draws, and
registered chain seeds `760121/760222/760323/760424`. Scientific retries with a
different seed are prohibited.

## Canonical data and QC

The source is the deterministic QGG `psoerensen/qgdata` simulated human panel at
commit `6cca5819e711d326cfb2614d7e9d9f34942612cd`. The source contains 5,000
individuals and 50,000 chromosome-1 markers. Existing Study 06 QC is applied
before any truth is generated: MAF 0.05, missingness 0.05, HWE `1e-12`, ambiguous
CG/AT removal, and duplicate/INDEL removal. Retained markers are ordered by
physical position and original PLINK index.

Contiguous blocks are constructed before annotations or phenotype using a
maximum of 500 markers. The summary route uses the same 5,000 individuals,
physical marker order, alleles, frequencies, and phenotype as BED. It requests
the retained low-rank representation, cumulative-positive-mass policy, and
`eigen_prop = 1 - .Machine$double.eps`; execution is permitted only if the audit
shows that every numerically positive mode is retained.

## Annotation and mixture truth

The annotation matrix contains an exact intercept, a deterministic 20% enriched
binary column, a standardized continuous signal correlated about 0.20 with the
binary column, and an independently generated residualized null column. The
annotation seed is `760101`.

The four components use `gamma = c(0, 0.01, 0.1, 1)` and target marginal
probabilities `0.970/0.015/0.010/0.005`. Non-intercept probit-stick truth is:

| Stick | Enriched | Continuous | Null |
|---|---:|---:|---:|
| 1 | 1.00 | 0.40 | 0 |
| 2 | 0.60 | 0.25 | 0 |
| 3 | 0.40 | 0.15 | 0 |

Intercepts are solved over the realized annotation matrix to maximum marginal
error at most `1e-6`. One truth seed, `760202`, samples components, effects, and
residuals once. A common effect multiplier preserves component variance ratios
while setting realized `h2 = 0.50` exactly. Failed truth gates stop the experiment;
they never trigger seed search or resampling.

## Priors and initialization

The learned models use the current small-Study-06 teaching hierarchy: proper
package-default intercept prior, no flat intercept, `sigmaSqAlpha_a = 2`,
`sigmaSqAlpha_b = 2`, initial annotation variances of one, probability floor
`1e-12`, and one hierarchy update per cycle. Learned non-intercept alpha starts at
zero; intercepts start from global target probabilities. Fixed-alpha models use
the exact calibrated truth and `updateAlpha = FALSE`. BayesR/SBayesR use the
matched target mixture initialization.

## Retention and decision contract

All alpha and `sigmaSqAlpha`, variance, component-probability, and annotation
prior summaries are required. A deterministic 300-marker panel stratifies truth
component, effect magnitude, causal status, enrichment, continuous signal, and
block; only this panel retains full effect, active, and component histories.
Genome-wide posterior means, SDs where supported, PIPs, and component
probabilities are retained as accumulators.

The continuation phase requests compact aggregate component-state histories
through `aggregate_component_states = TRUE`. These histories keep explicit
draw and chain dimensions and do not change the frozen specification, truth,
model, seed, or posterior target. Continuation checkpoints are prefixed
`continuation_` so they cannot overwrite the original non-inferential smoke
evidence.

The strict diagnostic thresholds are R-hat 1.01, bulk and tail ESS 400, and
relative MCSE 0.05. Full feasibility requires latent aggregate convergence; stable
PIPs alone cannot establish success. Decision classes are `LARGE-F1` through
`LARGE-F6` as specified in the registered task. The implementation is
[large-feasibility.R](../../studies/06_annotation_models/large-feasibility.R)
and the runner is
[run_study06_large_feasibility.R](../../scripts/run_study06_large_feasibility.R).

## Continuation block semantics

The block fits explicitly set `residual_policy = "gctb_block"`,
`block_ve_mode = "allMixVe"`, `resam_thresh = 1.1`, and
`minimum_ve_ratio = 0.7`. Under this contract, `fit$ves` is the retained mean
block residual variance and `fit$heritability_summary` is total block Vg divided
by phenotype variance. Neither is silently interpreted as the BED global
residual contract. Historical block fits retain the semantics of the package
SHA and residual policy that generated them.

The continuation stores compact component-count, realized-active-count, and
stick eligible/continue/stop histories plus the unchanged 300-marker panel.
Checkpoint identity includes package SHA, installed-tree hash, residual policy,
block Ve mode, frozen input hashes, seeds, and MCMC controls.
