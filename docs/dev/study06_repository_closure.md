# Study 06 repository closure

## Starting state

Closure began on 2026-08-07 from a clean `sblrbench` working tree on branch
`master` at `21ef774ec084d11932e132186ccc085054b9aa56`. The read-only sibling
`sblr` was clean at `7810222f622fb3381208a5aad33fd845ddd37d62`; its source and
installed package versions were both 0.2.0. R was 4.4.1 UCRT. Relevant package
versions included `posterior` 1.6.1, `coda` 0.19-4.1, `ggplot2` 3.5.2,
`testthat` 3.3.2, `jsonlite` 2.0.0, and `digest` 0.6.36.

The immutable identities were verified as:

| Identity | SHA-256 |
|---|---|
| Specification | `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56` |
| Informative truth | `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb` |
| Frozen truth export | `eb4efd63195874c21e273c11114b6144069afa5e841289af9e240a919c0e5013` |
| Learned BED checkpoint | `dea658dc4ad7736128319277d921a92f417c14ae02d4ce1dd5c6c5a05d82210f` |
| Learned block checkpoint | `954100e9a17d96c3ce235c53c75bfaaaee4587d8247bd83fefd3d8be51890495` |

Official marker-order, GWAS, annotation, LD-block, and D1 alpha-trace hashes,
plus compact PMA-R3 hashes, are recorded in
[`final_decision.json`](../../results/reference/06_annotation_models/final_decision.json).
No frozen artifact was changed.

## Final decisions

```text
Primary scientific decision: EST-R2
Official qualifier: EST-R5
Same-posterior sampler endpoint: PMA-R3 — exact but computationally impractical
Study 06 status: CLOSED
Additional final benchmark / scientific fit: not required
```

The official qualifier means formal independently seeded multichain replication
is unavailable under the pinned v0.2.6 native RNG contract. The official D1
trajectory remains useful descriptive evidence; it is not described as a failed
replication.

## Authoritative files audited

- [`studies/06_annotation_models/README.md`](../../studies/06_annotation_models/README.md)
- [`studies/06_annotation_models/report.qmd`](../../studies/06_annotation_models/report.qmd)
- [`docs/studies/study06_final_conclusion.md`](../studies/study06_final_conclusion.md)
- [`docs/studies/study06_estimability_and_contrasts.md`](../studies/study06_estimability_and_contrasts.md)
- [`docs/dev/study06_annotation_inference_evidence_synthesis.md`](study06_annotation_inference_evidence_synthesis.md)
- repository [`README.md`](../../README.md),
  [`studies/index.qmd`](../../studies/index.qmd),
  [`workflows.qmd`](../../workflows.qmd), and [`_quarto.yml`](../../_quarto.yml)
- Study 06 result manifests, registries, local final summaries, and the read-only
  sibling PMA-R3 endpoint documents.

The current machine-readable and compact tracked evidence is:

- [`final_decision.json`](../../results/reference/06_annotation_models/final_decision.json);
- [`final_cross_implementation_comparison.csv`](../../results/reference/06_annotation_models/final_cross_implementation_comparison.csv);
- [`final_hierarchy_of_evidence.csv`](../../results/reference/06_annotation_models/final_hierarchy_of_evidence.csv).

## Stale current status corrected

Current navigation previously described Study 06 as in development, the v2
final benchmark as pending/not authorized, and a focused sampler audit as the
next task. Those statements were correct at intermediate gates but stale as a
current roadmap. The root README, benchmark catalogue, workflow catalogue,
Study README, main report, and current synthesis/roadmap now identify EST-R2 as
authoritative, EST-R5 as the official qualifier, PMA-R3 as the sampler endpoint,
and additional Study 06 fits as unnecessary.

Historical phrases such as “final benchmark not authorized” remain where they
describe the v1/v2 qualification decisions, archived result documents, or
chronological evidence. They are explicitly historical rather than current.
Material retained hits include `study06_annotation_qualification_result.md`,
`study06_v2_design.md`, `study06_gctb_single_trajectory_design.md`, and
`study06_gctb_single_trajectory_result.md`. The older repository snapshots
`final_cleanup_plan.md`, `reorganization_dependencies.md`, and
`reorganization_baseline_capsules.csv` likewise describe their own historical
starting states. The new closure document's mentions of stale phrases are audit
records, not roadmap instructions.

## Historical evidence preserved

The closure did not rewrite the original decisions. The retained sequence is:

```text
initial qualification / stop decisions
→ annotation-identifiability qualification
→ paired power and isolation analyses
→ official SBayesRC descriptive comparison
→ GCTB-compatible block residual contract validation
→ large information-scale feasibility and LARGE-G2
→ BLOCK-MIX-R4
→ AA-R5 / coordinated-transition result
→ PMA-R3 sampler-development endpoint
→ estimability and annotation-contrast analysis
→ EST-R2 scientific closure
```

## Final scientific summary

### Raw first-stick alpha

| Coefficient | Truth | Official D1 | `sblr` BED | `sblr` block |
|---|---:|---:|---:|---:|
| intercept | -1.653 | -2.406 | -2.934 | -2.201 |
| enriched | 1.600 | 1.507 | 1.897 | 1.642 |
| continuous | 0.300 | 0.586 | 0.638 | 0.622 |
| null | 0 | -0.156 | -0.042 | -0.371 |

The official value is single-trajectory descriptive; `sblr` values are selected-
window pooled learned summaries under EST-R2, not converged truth estimates.
Informative directions agree qualitatively, but quantitative equality is not
established. Stick 1 is scientifically central because
`P(c_i>0 | A_i) = Phi(A_i alpha_1)`.

The frozen eligible counts are 1,500/171/87 across the three sticks. The
enriched annotation is present in 225/1,500 SNPs and in 107/171 and 62/87 of the
later eligible sets. Later-stick raw values remain weak because eligible
continuation outcomes shrink and are latent, with strong alpha/allocation
dependence. Official stick-2/stick-3 intercepts near 10.56/8.95 versus truths
near -0.235/-0.431 and stick-3 intercept ESS near 1.4 are descriptive evidence,
not an official multichain result.

### Annotation variance

| Stick | Official D1 final `sigma_anno` | BED `sigmaSqAlpha` mean | Block `sigmaSqAlpha` mean | Fixed ablation |
|---|---:|---:|---:|---:|
| 1 | 1.976 | 2.137 | 2.044 | 1 |
| 2 | 1.043 | 1.897 | 1.540 | 1 |
| 3 | 4.803 | 3.136 | 1.788 | 1 |

This is descriptive scale consistency only. Official final `sigma_anno` and
retained `sblr` `sigmaSqAlpha` summaries are not established as identical
estimands; official variance histories are unavailable; fixed values are
ablations rather than posterior estimates.

### Active annotation contrasts

| Annotation | Truth | `sblr` BED | `sblr` block | Official D1 |
|---|---:|---:|---:|---:|
| enriched | .423 | .192 [.081, .397] | .299 [.119, .642] | .204 [.117, .302] |
| continuous | .0886 | .0555 [.0239, .104] | .112 [.0383, .236] | .077 [.042, .118] |
| null | 0 | -.007 [-.068, .046] | -.068 [-.268, .060] | -.024 [-.096, .024] |

All four evidence sources support enriched much greater than continuous, with
null near zero. Direction and ranking are stable; exact magnitude is uncertain
and route-sensitive.

### Probability/ranking and SNP stability

- Active-prior Spearman is at least .986 BED and .950 block;
  largest-component-prior Spearman is at least .999 and .995.
- Official half-specific `A alpha` vectors can correlate at -.975, yet induced
  nonconstant probability vectors correlate at least .999.
- PIP Pearson is at least .994/.969, signed-beta Pearson .9999/.9991, and
  genetic-value chain correlation .99994/.99925 for BED/block.
- BED chain ranges are AUPRC .591–.597, AUROC .849–.856, genetic-value truth
  correlation .917–.918, and prediction .654–.655. Block ranges are
  .514–.560, .818–.840, .896–.897, and .631–.633.

## Why SBayesRC can work despite difficult alpha parameters

### Post-closure ALPHA-R1 evidence

An isolated standard-SBayesRC probit reference was added after repository
closure. It leaves EST-R2, EST-R5, and PMA-R3 unchanged. Fixed-`z` scalar Gibbs
matched the exact Gaussian alpha posterior (maximum mean/covariance errors
0.00428/0.00084); blocked and scalar known-outcome fits agreed (maximum mean
difference 0.00936); and 20 repeated fits had coverage 0.90--1.00 and maximum
R-hat 1.00199. Enriched-slope posterior SD rose .118, .198, .279 as mean eligible
counts fell 1,500, 170.7, 84.2. Enriched eligible counts were 149.1, 77.2, 42.3,
so selective continuation raised the later eligible fractions to about 45% and
50%.

```text
known compatible hierarchy -> alpha inference works
unknown hierarchy -> alpha changes SNP probabilities
SNP allocations/effects -> define the next alpha regression
coupled joint state -> local Gibbs exploration can be slow
```

PMA-R3 remains the endpoint: coordinated global alpha/allocation movement is
exactly possible but computationally impractical. The Jian Zeng 2024 R-reference
total-`m` intercept detail produced only modest isolated differences from the
eligible-`n` calculation and is retained as a source-contract note, not called
an error or attributed to production GCTB. The tracked script is
`studies/06_annotation_models/alpha_recovery_reference.R`; compact evidence is
under `results/reference/06_annotation_models/alpha_recovery_*`. No production
sampler or frozen Study 06 artifact is involved.

For each stick,

\[
\eta_{ik}=A_i\alpha_k,\qquad q_{ik}=\Phi(\eta_{ik}),
\]

and

\[
\pi_{i0}=1-q_{i1},\quad
\pi_{i1}=q_{i1}(1-q_{i2}),\quad
\pi_{i2}=q_{i1}q_{i2}(1-q_{i3}),\quad
\pi_{i3}=q_{i1}q_{i2}q_{i3}.
\]

```text
alpha → A alpha → q → pi → components/effects → PIP/beta → prediction
```

The design is full rank (rank 4, condition number about 2.23), so Study 06 does
not claim structural non-identifiability. Instead, some latent coefficient
directions are weakly determined and strongly correlated. Compensating alpha
configurations can behave similarly over realized annotations, and the probit
map can compress large latent differences (`Phi(4)` is about .99997 while
`Phi(7)` is effectively 1). Official D1's poor half-to-half `A alpha` agreement
but near-perfect probability-vector agreement illustrates this compression.

MCMC then averages over hierarchy draws:

\[
P(c_i>0\mid y)=E_{\alpha,\theta}[P(c_i>0\mid y,\alpha,\theta)].
\]

The SNP result is not conditioned on one plug-in alpha estimate. This explains
how difficult upstream parameters can coexist with reproducible marginal SNP
inference on this fixture. It does not make poor mixing harmless in general.

### Three levels of estimability

| Level | Result |
|---|---|
| latent alpha / hierarchy scale | weak or difficult in some directions, especially later sticks |
| annotation probability function | substantially more stable; strong direction/ranking, uncertain magnitude |
| marginal SNP inference | highly stable |

R-hat near one is chain agreement, not truth recovery. Expected active count
is about 351.8 BED and 807.4 block versus truth 180 despite near-one R-hat;
stable ranking likewise does not establish calibrated probability magnitude.

## Reporting recommendation

Primary SNP reports should include PIPs, posterior effects, and downstream
genetic-value/prediction performance. Primary annotation reports should include
active/component probabilities, counterfactual contrasts with intervals and
directional posterior probability, and SNP-prior rankings. Raw alpha,
`A alpha`, `sigmaSqAlpha`, posterior covariance/correlation, R-hat, ESS, and
MCSE should be retained as secondary latent-model diagnostics. Raw alpha should
not be the sole or primary biological summary.

## Implementation-validation statement

`sblr` implements the core standard SBayesRC model and has been lightly but
systematically validated against the official implementation and independent
mathematical reference experiments on controlled simulation benchmarks. The
support is the validated residual contract,
baseline/fixed-alpha controls, direct small official comparisons, first-stick
direction and derived-contrast agreement, high SNP-level concordance, exact
alpha-conditional/known-outcome calibration, and transition/oracle audits. The
claim is controlled-benchmark concordance, not
universal GCTB equivalence.

## Reproduction

```powershell
Rscript scripts/run_study06_estimability.R
```

Retained frozen chains must already exist locally. This performs offline
posterior analysis only, invokes no sampler, regenerates no truth, and writes to
the ignored `results/local/06_annotation_models/estimability_and_contrasts/`
directory.

## Next study

The next methodology is a separate Bayesian annotation-selection /
annotation-PIP study. It is motivated by the estimand—annotation relevance and
stable functional summaries—not merely by computational difficulty. This
closure does not implement that model.
