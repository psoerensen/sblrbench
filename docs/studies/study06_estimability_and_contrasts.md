# Study 06 estimability and annotation contrasts

## Decision

The primary classification is **EST-R2 — probability/ranking functions are
stable but annotation contrasts remain uncertain**. The induced SNP-prior
rankings, SNP PIPs, posterior effects, genetic values, and predictions are much
more reproducible than the unrestricted continuous annotation decomposition.
Counterfactual annotation contrasts mix more cleanly than many raw directions
and recover the expected directions, but their magnitudes remain biased,
route-sensitive, and Monte Carlo-uncertain. This does not support strong
quantitative annotation-effect claims.

The pinned official arm is secondarily qualified as **EST-R5 — formal official
multichain replication remains unavailable**. Its one retained trajectory shows the same
qualitative compression from raw alpha movement to induced probabilities, but
official v0.2.6 cannot produce independently seeded native trajectories.

The authoritative closure is recorded in
[`final_decision.json`](../../results/reference/06_annotation_models/final_decision.json)
and interpreted in the [final conclusion](study06_final_conclusion.md). This
document remains the detailed final evidence source.

## Provenance and immutable inputs

The analysis began from clean `sblrbench` `master` at
`fbe80603ff6fa09e0a611a56d09130cb4b2cbc8c` and clean read-only sibling `sblr`
at `3c7b97f3de76a6e19a7e82bf73b2b2b10bf83d34`. The sibling contains the PMA-R3
decision, “exact but computationally impractical,” and the standard SBayesRC
sampler-development endpoint. Source and installed `sblr` are both 0.2.0.

R was 4.4.1 UCRT. Relevant package versions were `posterior` 1.6.1, `coda`
0.19-4.1, `data.table` 1.18.4, `ggplot2` 3.5.2, `dplyr` 1.1.4, `arrow` 22.0.0,
and `digest` 0.6.36.

The immutable Study 06 identities are:

| Identity | SHA-256 |
|---|---|
| Specification | `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56` |
| Informative truth | `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb` |
| Frozen truth export file | `eb4efd63195874c21e273c11114b6144069afa5e841289af9e240a919c0e5013` |
| Learned BED checkpoint | `dea658dc4ad7736128319277d921a92f417c14ae02d4ce1dd5c6c5a05d82210f` |
| Learned block checkpoint | `954100e9a17d96c3ce235c53c75bfaaaee4587d8247bd83fefd3d8be51890495` |

Each primary checkpoint has four chains and 9,000 retained post-burn draws.
The exact paths, sizes, official alpha-file hashes, and frozen identities are
in `analysis_input_manifest.csv` under the ignored local result directory.
No scientific fit was rerun and no truth was regenerated.

## Source-derived evidence

The retained likelihood is already aligned to the official block residual
variance contract. Baseline BayesR/SBayesR and fixed-true-alpha BayesRC/SBayesRC
controls converge. The learned-alpha joint hierarchy fails the established
qualification contract. Audited ordinary, repeated, tempered, pair,
coordinated, PX, particle-Gibbs, and particle-marginal transitions are not
reopened here.

The production mapping was used exactly:

```text
eta = A alpha
q_k = Phi(eta_k)
pi_0 = 1-q_1
pi_1 = q_1(1-q_2)
pi_2 = q_1 q_2(1-q_3)
pi_3 = q_1 q_2 q_3
```

Across every checked draw, probabilities were finite and in `[0,1]`; the
maximum row-sum error was `2.22e-16`.

## New analysis

### Raw alpha and posterior directions

The retained coefficient marginals are not uniformly catastrophic by R-hat,
but they have low effective sample sizes, material MCSE, and poor truth
recovery. This is consistent with poor practical exploration and weakly
determined coefficient decompositions; it is not a proof of formal
non-identifiability.

| Route | Worst R-hat | Median R-hat | Min bulk ESS | Max relative MCSE | Max absolute truth error |
|---|---:|---:|---:|---:|---:|
| BED | 1.025 | 1.006 | 81.1 | 0.132 | 1.222 |
| Block | 1.040 | 1.011 | 70.3 | 0.127 | 0.893 |

The annotation design itself is full rank (rank 4), has condition number 2.23,
and singular values 39.23, 39.08, 38.35, and 13.65. Thus the observed behavior
is not explained by rank deficiency. Pooled alpha correlations show
compensation, especially intercept/enriched tradeoffs (BED stick 1 correlation
−0.667) and continuous/enriched or continuous/intercept tradeoffs. The largest
posterior covariance directions are dominated by enriched coefficients, with
intercept contributions on earlier sticks. Largest eigenvalues rise to 4.49
for BED stick 3 and 1.60 for block stick 3.

### A alpha and prior probabilities

Representative `A alpha` scalars improve to worst R-hat 1.014 (BED) and 1.011
(block). Whole-vector chain agreement is uneven: minimum Pearson correlation is
0.691 for BED and 0.862 for block, with later sticks least stable. Truth RMSE
remains large (0.72–1.31), so improved chain diagnostics do not imply accurate
latent linear predictors.

The nonlinear probability map compresses much of this disagreement:

| Route/quantity | Minimum Pearson | Minimum Spearman | Maximum RMSE | Minimum top-100 overlap |
|---|---:|---:|---:|---:|
| BED active prior | 0.988 | 0.986 | 0.014 | 0.88 |
| Block active prior | 0.954 | 0.950 | 0.054 | 0.83 |
| BED largest-component prior | approximately 1.000 | 0.999 | <0.01 | 0.97 |
| Block largest-component prior | 0.99+ | 0.995 | <0.02 | 0.93 |

Component-probability truth RMSE is 0.021–0.119 for BED and 0.018–0.077 for
block. Magnitudes nevertheless remain unstable: expected active count has mean
351.8, SD 816.1 for BED and mean 807.4, SD 1,263.3 for block, versus truth
180.0. Its R-hat is near one because all chains explore the same heavy-tailed
functional; this is an important example of convergence not implying recovery.

### Counterfactual annotation contrasts

The binary contrast sets enriched status to 1 versus 0 for every SNP while
holding other annotations fixed. Continuous and null contrasts use the frozen
standardized `+1 SD` versus `−1 SD` comparison.

| Route | Annotation | Mean active contrast | 95% interval | P(Delta > 0) | R-hat | Bulk ESS | Truth |
|---|---|---:|---:|---:|---:|---:|---:|
| BED | enriched | 0.192 | [0.081, 0.397] | 1.000 | 1.014 | 307 | 0.423 |
| BED | continuous | 0.0555 | [0.0239, 0.104] | 0.999 | 1.008 | 297 | 0.0886 |
| BED | null | −0.0070 | [−0.0676, 0.0461] | 0.393 | 1.012 | 204 | 0 |
| Block | enriched | 0.299 | [0.119, 0.642] | 1.000 | 1.029 | 146 | 0.423 |
| Block | continuous | 0.112 | [0.0383, 0.236] | 0.998 | 1.009 | 160 | 0.0886 |
| Block | null | −0.0680 | [−0.268, 0.0596] | 0.154 | 1.002 | 84 | 0 |

The enriched direction is decisive in both routes, the continuous direction is
strong, and the null interval covers zero. However, the BED enriched interval
misses truth, maximum contrast bias is 0.231 (BED) and 0.124 (block), relative
MCSE reaches 0.069 and 0.108, and route-level magnitudes differ. This is why the
result is EST-R2 rather than EST-R1.

The secondary observed-group active differences are 0.189 [0.081, 0.391] for
BED and 0.292 [0.117, 0.635] for block. These descriptive contrasts are not
counterfactual effects because correlated annotation profiles can contribute.

### Ranking and SNP-level stability

| Route/quantity | Minimum Spearman | Minimum top-50 overlap | Minimum top-100 overlap |
|---|---:|---:|---:|
| BED prior active | 0.986 | 0.92 | 0.88 |
| Block prior active | 0.950 | 0.74 | 0.83 |
| BED largest-component prior | 0.999 | 0.98 | 0.97 |
| Block largest-component prior | 0.995 | 0.96 | 0.93 |
| BED PIP | 0.952 | 0.90 | 0.91 |
| Block PIP | 0.936 | 0.86 | 0.82 |
| BED absolute beta | 0.842 | 0.98 | 0.96 |
| Block absolute beta | 0.951 | 0.88 | 0.89 |

PIP Pearson correlations are at least 0.994 (BED) and 0.969 (block); signed
effect correlations are at least 0.9999 and 0.9991. Validation genetic-value
correlations between chains are at least 0.99994 and 0.99925. Chain ranges are:

| Route | AUPRC | AUROC | Genetic-value correlation to truth | Phenotype prediction correlation |
|---|---:|---:|---:|---:|
| BED | 0.591–0.597 | 0.849–0.856 | 0.917–0.918 | 0.654–0.655 |
| Block | 0.514–0.560 | 0.818–0.840 | 0.896–0.897 | 0.631–0.633 |

### Official SBayesRC arm

The reused official implementation is SBayesRC 0.2.6 at
`b95d3fcbad8ff358290922a58fff893439296138`. D1 has 6,000 post-burn alpha draws.
Raw alpha single-chain ESS spans 1.4–86.2 and the largest half-to-half mean shift
is 1.25 posterior SD. Half-specific posterior-mean `A alpha` vectors can disagree
severely (minimum Pearson −0.975), while induced component-probability vectors
have minimum nonconstant Pearson correlation 0.999. Active counterfactual
contrasts are 0.204 [0.117, 0.302] for enriched, 0.077 [0.042, 0.118] for
continuous, and −0.024 [−0.096, 0.024] for null; their maximum half shift is
0.12 SD.

This is qualitatively consistent with stronger estimability after the
probability map. It is not independent replication: the public `seed=` does not
control all native static/thread-local RNG state, and distinct fresh processes
were already shown to repeat the same scientific trajectory. No R-hat is
reported and no native RNG was modified. Current `sblr` also uses a proper
default intercept prior and `sigmaSqAlpha` shape/rate 2/2 in these fits; the
official hierarchy has different prior and empty-stick conventions and does
not expose retained `sigmaSqAlpha` history. Exact posterior equality is neither
expected nor claimed.

### PMA-R3 secondary diagnostic

The compact PMA-R3 alpha, convergence, and summary exports were copied with
hash equality into `evidence/pma_r3/`; `provenance.csv` records source paths,
the current sibling head, and the generating commit. The short reference had
acceptance 0.717–0.777, alpha R-hat 1.27–3.34, active-count R-hat 1.14, and
expected-active R-hat 1.22. It shows that a valid global move can coordinate
alpha and occupancy, but 200 retained draws per chain are not qualification
evidence and do not change EST-R2.

## Interpretation

Study 06 separates three claims. The exact continuous alpha decomposition is
weakly determined in practice and costly to explore. Induced SNP-prior rankings
are reproducible. Directional biological contrasts are interpretable and
promising, but their magnitudes are not yet reliable enough for strong
annotation-level claims. SNP prioritization, effects, genetic values, and
prediction remain scientifically useful despite hierarchy-level uncertainty.

## Limitations

- R-hat diagnoses chain agreement, not formal identifiability or truth recovery.
- The frozen fixture is one simulated replicate with four annotations.
- Expected active count is heavy-tailed and poorly recovered despite excellent R-hat.
- Official evidence is one deterministic native trajectory, not a multichain replication.
- BED and block learned-alpha routes have distinct likelihood representations and differ quantitatively.
- PMA-R3 is short diagnostic evidence only.

All detailed tables and figures are under
`results/local/06_annotation_models/estimability_and_contrasts/`.
