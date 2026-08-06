# Study 06 official SBayesRC parity design

## Purpose and immutable identity

This external-reference diagnostic asks whether official
`zhilizheng/SBayesRC` behaves like `sblr` on the exact Study 06 v2 informative
truth. It is neither a qualification nor permission to change `sblr` or run
the final benchmark.

The frozen specification hash is
`241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56`;
the truth hash is
`169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb`.
The truth has 2,000 people (1,400 training/600 validation), 1,500 ordered
markers, 15 blocks of 100, component counts 1,329/84/50/37, h2 0.50, and
`gamma=c(0,.01,.1,1)`.

## Official dependency

The pinned official release is v0.2.6 at commit
`b95d3fcbad8ff358290922a58fff893439296138`. Source is cloned only under the
ignored parity root and installed into
`results/local/06_annotation_models/gctb_parity/rlib`. It is never a CI or
normal test dependency. The source declares C++11, while Boost shipped by
current BH requires C++14; the unchanged source was therefore built with an
ignored local Makevars mapping `CXX11STD=-std=gnu++14`, BH 1.87.0-1, GCC 13.2,
and `-fopenmp`. No ordinary R library was modified.

The installed DESCRIPTION has version 0.2.6 and the official GitHub URL but no
`RemoteSha`; source-checkout HEAD, tag, clean status, install command, and the
manifest pin jointly establish provenance.

## Source-audited official contract

| Item | Official v0.2.6 behavior |
|---|---|
| defaults | `gamma=c(0,.001,.01,.1,1)`, `startPi=c(.990,.005,.003,.001,.001)`, 3,000 iterations, 1,000 burn-in |
| summary scaling | With `twopq="nbsq"`, sets phenotype variance to 1, uses `ord_std=sqrt(N*se^2+b^2)`, and analyzes `b/ord_std` |
| intercept | Constructed internally and sampled under an explicitly flat prior |
| non-intercept variance | Per stick, `(sum(alpha[-1]^2)+4)/chi-square((numAnno-1)+4)`; scale argument is not used in that source expression |
| residual variance | `allMixVe`: block-specific eligibility based on SNP-effect/LD variance ratio; inverse-chi-square sample accepted only above 0.7 of phenotype variance, otherwise reset to phenotype variance |
| tuning | Pseudo train/validation summary split; evaluates `.995/.99/.95/.9`; usually retains max threshold unless relative validation gain exceeds 1.25 or its correlation is negative |
| PIP | Output is exactly `1 -` posterior null-component assignment frequency |
| alpha | One text trace per stick, columns in annotation-header order |
| annotation probabilities | Conditional and joint component-probability traces are text outputs |
| component histories | `n_comp_hist`, `vg_comp_hist`, `pi_hist` returned in RDS |
| h2 | `hsq_hist`; total genetic variance is block-factor quadratic form divided by unit phenotype variance |
| sigma-alpha trace | Only final `sigma_anno` is returned; a retained per-iteration trace is not exposed |
| seed | Public `seed` formal is declared but is not used by the wrapper or passed to native code; native static/thread-local C++ engines have no seed setter in v0.2.6 |

The official annotation reader requires an explicit second `Intercept` column:
it allocates `ncol(input)-1` annotation columns, overwrites its first internal
column with ones, and reads input columns 3 onward. The export therefore uses
`SNP, Intercept, enriched_binary, continuous_signal, null_annotation`; this is
not a duplicated fitted intercept.

## Deterministic export

`studies/06_annotation_models/gctb-parity.R` reconstructs the exact informative
replicate and writes:

- COJO-style summary statistics with SNP, alleles, training allele frequency,
  marginal beta, standard error, p-value, and `N=1400`;
- the official annotation file described above;
- `ldm.info`, `snp.info`, and 15 official binary eigen blocks computed from the
  exact training-standardized genotypes;
- a local truth RDS and SHA-256 manifest.

The GWAS cross-products are checked against the existing Study 06 sufficient
statistics. LD blocks must be symmetric, have 100 markers, reconstruct through
the official reader within `1e-5`, and preserve marker/allele order. Every
positive mode is stored with file threshold 1.0. G0/G1 request threshold 1.0;
G2 may tune the reader threshold. This is exact within the official float
binary representation, not an external or approximate LD reference.

## Registered diagnostic

Every condition has four intended chains derived deterministically from Study
06 chain seeds 701121/701222/701323/701424:

| Condition | Model | Gamma/start | Tuning | Iterations |
|---|---|---|---|---:|
| G0 | official SBayesR, no annotation | `0,.01,.1,1`; `.88,.06,.036,.024` | off; threshold 1 | 3,000/1,000 burn |
| G1 | matched four-component official SBayesRC | same | off; threshold 1 | 3,000/1,000 burn |
| G2 | native five-component official SBayesRC | `0,.001,.01,.1,1`; `.990,.005,.003,.001,.001` | on, 150/100, grid `.995/.99/.95/.9` | 3,000/1,000 burn |

All runs use `starth2=.5`, `bTunePrior=FALSE`, `sSamVe="allMixVe"`,
`twopq="nbsq"`, detail output, beta histories, and output frequency 1. These
settings expose 2,000 retained states if the registry is executable.

Parity is approximate, not model identity: official G1 differs from `sblr` in
flat versus proper intercept prior, annotation-variance prior, residual update,
summary likelihood/scaling, update order, numerical safeguards, and native RNG.

## Smoke and stop rule

Before the registry, the runner performs 40-iteration/20-burn smokes for G0,
G1, and tuned G2. Two G0 smokes run in fresh R processes with different
preregistered requested seeds. If their scientific histories are identical,
the public official build cannot supply independent chains and the registry
stops as GCTB-P5. It is forbidden to patch the official source, vary settings,
or reinterpret duplicate streams as four chains.

Only if the seed audit passes may the 12 coordinates run. Short-run diagnostics
would be descriptive, using Study 06 thresholds without treating them as formal
qualification. SNP alignment, PIP/effect/rank correlations, top-k overlap,
Bayesian-FDR overlap, alpha directions, occupancy, h2, variance, ranking, and
prediction would then be compared with the committed `sblr` references.

## Decision rule

GCTB-P1 through P4/P6 require four genuinely independent chains. Package,
export, parser, numerical, runtime, or independent-chain incompatibility is
GCTB-P5. A blocked result preserves the validated export and smoke evidence but
must not report convergence or scientific parity from a single duplicated RNG
stream.
