# Study 06 large information-scale feasibility result

## Decision

**Primary: LARGE-G2 — learned allocation feedback remains the bottleneck.**

The aggregate latent and variance quantities pass the registered convergence
thresholds in E0/B0 and E2/B2, while both ordinary learned-alpha fits fail
alpha, occupancy, and active/stick-count convergence. B1 also fails genetic and
base-variance convergence. The result supports two secondary flags:

- **LARGE-G3:** SNP effects, ranking, and genetic-value recovery remain useful
  despite learned latent nonconvergence.
- **LARGE-G4:** material BED/block route differences persist in reported
  heritability, Vg, occupancy, and learned marker priors.

This is a single-replicate feasibility result, not a qualification. The formal
small Study 06 qualification remains failed and the final benchmark remains
unauthorized.

## Historical blocker and continuation

The first registered phase ended as `LARGE-F6`: compact allocation histories
were unavailable and B0 rejected a materially negative global-projected
residual scale at iteration zero. That historical result remains intact.

Package-side work subsequently added exact, RNG-neutral compact component and
stick histories and a GCTB-compatible block-local residual policy. Before this
continuation, a pinned small official comparison passed as `SMALL-G1`. The new
large continuation used `sblr` 0.2.0 from clean source SHA
`0c89234273389e14112ba0e08ef9d11d3e1819dc`, installed-tree SHA-256
`e723528e7d5d570a31b5b1d1c90551896ac48f86ab05261c181c8109af971fd0`.
Historical `global_projected` fits were not reinterpreted.

## Frozen identity

No simulation object was regenerated. The continuation reverified all frozen
identities before smoke or scientific execution.

| Identity | SHA-256 |
|---|---|
| specification | `b001bc36a5531e5e6b342286a253fc1fd34dad4265359d89d2feaa026d4533df` |
| truth | `e94a511540f600e61ef47b52947836f19a15388f5e8ce795c179929956817507` |
| marker order | `e394f4324f89ca7ad88691284a6216a24e26c691d7cb59d822d5d7f006b096f2` |
| sample order | `f3c14e98845fc565424a973a5a5850768c587f9d5a5306cc58ff36c5329a256b` |
| alleles | `996fb56147c14801f798a7bb21692a962aba18829c8a3a7abfa4dc768b06b082` |
| annotations | `5cc2a5dc64140b5cb2c3e77044d8a0b7bd698e2b33f1f8cb0f5b39ac278f7bd3` |
| blocks | `3b119c38fcababc5b70892f8fc29dcd773fb73e4e2a2b1cd6daf94556e16294d` |
| GWAS | `f00c01326d20d32c1b389f239ebe520327f4b8c859acf7320b10d531243cacfd` |
| alpha truth | `4766d00b77653825e9130a32ebcde1b16754ee99103f2ab4c4d3f1d715fbbf82` |
| selected panel | `0ae8cc37d0418d54cf52e4cf5271c5859d01506759a6d798f6c512a66b08438f` |

The data contain 5,000 people, 37,991 markers, and 76 blocks (75 by 500 and one
by 491). Truth occupancy is 36,791/618/392/190, hence 1,200 active markers.
Stick eligible/continuation counts are 37,991/1,200, 1,200/582, and 582/190.
Realized Vg and Ve are both one and h2 is exactly .50. All registered
truth-identifiability, marker, allele, GWAS, annotation-rank, and full-positive-
mode block gates passed.

## Execution contract and smokes

All six noninferential 12-draw, four-chain smokes passed (`SMOKE-G1`) with the
registered seeds and compact/selected traces. E2/B2 held alpha exactly at truth.
B0/B1/B2 explicitly used `residual_policy = "gctb_block"`,
`block_ve_mode = "allMixVe"`, `resam_thresh = 1.1`, and
`minimum_ve_ratio = 0.7`; no global projected SSE was evaluated.

Every scientific fit used four chains with seeds
760121/760222/760323/760424, 12,000 total iterations, 3,000 burn-in, 9,000
retained iterations, thinning one, the proper alpha intercept prior,
`sigmaSqAlpha_a = sigmaSqAlpha_b = 2`, probability floor `1e-12`, and the
ordinary one-allocation/one-hierarchy schedule. No changed-seed retry,
extension, substitution, or additional replicate occurred.

| Fit | Model | Runtime | Checkpoint |
|---|---|---:|---:|
| E0 | BED BayesR | 6,475 s (107.9 min) | 497.9 MB |
| B0 | block SBayesR | 1,701 s (28.3 min) | 499.2 MB |
| E2 | fixed-alpha BED BayesRC | 6,070 s (101.2 min) | 500.0 MB |
| B2 | fixed-alpha block SBayesRC | 991 s (16.5 min) | 505.5 MB |
| E1 | learned BED BayesRC | 5,639 s (94.0 min) | 500.0 MB |
| B1 | learned block SBayesRC | 1,486 s (24.8 min) | 505.5 MB |

The compact aggregate trace is about 2 MB per fit, versus roughly 14.6 GB for
a full all-marker component history. The frozen 300-marker panel was retained.

## Convergence

Thresholds were R-hat <= 1.01, bulk and tail ESS >= 400, and relative MCSE <=
.05. “Core failures” below exclude selected-marker diagnostic limitations but
include all registered aggregate, variance, and annotation quantities.

| Fit | Core failures | Max R-hat | Min bulk ESS | Min tail ESS | Max rel. MCSE | Interpretation |
|---|---:|---:|---:|---:|---:|---|
| E0 | 0 | 1.0022 | 2,936 | 3,042 | .0185 | core contract passed |
| B0 | 0 | 1.0013 | 2,304 | 4,528 | .0208 | core contract passed |
| E2 | 0 | 1.0008 | 3,854 | 8,145 | .0161 | core contract passed; alpha fixed |
| B2 | 0 | 1.0016 | 2,573 | 5,235 | .0198 | core contract passed; alpha fixed |
| E1 | 23 | 1.3635 | 9.95 | 24.0 | .339 | failed |
| B1 | 31 | 1.9629 | 5.51 | 11.4 | .454 | failed |

E1 failures comprise 14/15 annotation quantities, three component counts,
realized active count, one stick-eligible, one stick-continuation, and all three
stick-stopping counts. B1 fails all 15 annotation quantities, all four component
counts, realized active count, two stick-eligible, all stick-continuation and
stick-stopping counts, and base/genetic/component-genetic variance quantities.

Some selected-marker binary/component histories in every fit have one or more
constant chains. Their rank/tail diagnostics are explicitly unavailable; they
are not replaced with final states. The complete 300-marker failure table is in
the local compact analysis. This limitation does not affect the exact aggregate
trace results above.

The package did not retain a draw-wise genome-wide annotation-prior expected-
active trace. Posterior expected occupancy is obtained exactly by summing the
posterior marker component probabilities and agrees with the compact realized-
active mean. Annotation-prior expected activity is reported only as a plug-in
at mapped pooled alpha means and is not used as a convergence gate.

## Alpha and annotation variance recovery

Neither learned route meets the convergence precondition for declaring
parameter recovery. Descriptively, all six informative non-intercept means have
the correct sign and all six 95% pooled intervals contain truth on each route;
all three null-annotation intervals contain zero.

| Route | Non-intercept median/max abs. error | Informative sign | Informative coverage | Null coverage |
|---|---:|---:|---:|---:|
| E1 BED | .278/.529 | 6/6 | 6/6 | 3/3 |
| B1 block | .102/.339 | 6/6 | 6/6 | 3/3 |

The pooled means cannot be interpreted as converged recovery. E1 first-stick
means (truth in parentheses) are intercept -3.552 (-2.472), enriched 1.200
(1.0), continuous .293 (.4), and null -.116 (0). B1 gives -2.658, 1.179,
.162, and -.198. Later-stick intercepts depart strongly and chain locations
differ.

`sigmaSqAlpha` is broad and mostly nonconverged. E1 means are
1.195/2.163/1.957; only stick 1 passes. B1 means are 1.208/1.340/1.400 and none
passes. Only three non-intercept coefficients inform each variance, so width by
itself is not failure; here the R-hat/ESS/MCSE failures establish poor mixing.
The value one is the registered initial/reference value, not a randomly drawn
simulation truth for these prespecified alpha coefficients.

## Occupancy, variance, and block residual behavior

| Fit | Posterior component counts 0/1/2/3 | Realized active | Expected active from posterior component probabilities | Vg | h2 |
|---|---|---:|---:|---:|---:|
| E0 | 36919/565/347/159 | 1,071.8 | 1,071.8 | .713 | .349 |
| B0 | 36848/570/380/194 | 1,142.9 | 1,142.9 | 1.054 | .515 |
| E2 | 36903/567/352/169 | 1,087.9 | 1,087.9 | .761 | .372 |
| B2 | 36849/570/379/193 | 1,142.3 | 1,142.3 | 1.049 | .513 |
| E1 | 37847/25/9/110 | 144.4 | 144.4 | .719 | .352 |
| B1 | 36898/655/172/266 | 1,092.7 | 1,092.7 | 1.048 | .512 |

Truth is 36,791/618/392/190 and 1,200 active. Controls broadly recover total
activity, though E0/E2 undercount and their large-component means are 16% and
11% below truth. E1 collapses to a very sparse, large-component-heavy regime;
B1 is closer in total activity but has broad, chain-separated small-component
occupancy. These learned summaries are descriptive under nonconvergence.

BED `ves` is global individual-level residual variance: 1.332 (E0), 1.283
(E2), and 1.324 (E1). Block `ves` is the retained mean block residual variance,
not a BED-equivalent quantity. Under `allMixVe`, no block met the registered
resampling trigger in any large fit; every retained/final block Ve remained at
the initialized phenotype variance 2.045, with zero resamples and zero
minimum-ratio resets. Complete block-Ve histories were not requested.

Reported block heritability is total block Vg divided by phenotype variance.
The BED-minus-block h2 differences are -.166 (baseline), -.141 (fixed alpha),
and -.161 (learned). The offset therefore persists at large scale and remains a
route-contract issue separate from annotation learning. It must not be inferred
from a direct BED-versus-block `ves` comparison.

## SNP-level utility

| Fit | AUPRC | AUROC | Top-50 causal | Top-100 causal | Effect-truth cor. | Genetic-g cor. | Phenotype cor. |
|---|---:|---:|---:|---:|---:|---:|---:|
| E0 | .151 | .618 | 49 | 84 | .783 | .835 | .711 |
| B0 | .134 | .601 | 48 | 78 | .769 | .838 | .748 |
| E2 | .304 | .872 | 50 | 90 | .829 | .866 | .716 |
| B2 | .305 | .872 | 50 | 88 | .814 | .859 | .736 |
| E1 | .312 | .863 | 50 | 93 | .809 | .834 | .644 |
| B1 | .254 | .840 | 50 | 87 | .807 | .855 | .732 |

At Bayesian FDR 5%, selected/true/false counts are E0 42/42/0, B0 55/53/2,
E2 60/59/1, B2 78/74/4, E1 55/54/1, and B1 77/75/2. At 10% they are
50/49/1, 67/62/5, 73/72/1, 96/86/10, 65/64/1, and 95/84/11.

Annotations improve AUPRC over matched no-annotation baselines on both routes
(.312 vs .151 BED; .254 vs .134 block), and fixed-true-alpha results are also
substantially better than baselines. These learned comparisons are useful but
descriptive under nonconvergence; better ranking does not compensate for the
failed joint posterior contract. Genetic-value recovery is in-sample against
known simulated g because no independent registered sample exists.

## BED/block comparisons

| Pair | PIP Pearson/Spearman | Effect Pearson | Genetic-g agreement | Top 50/100 overlap | BED-block h2 |
|---|---:|---:|---:|---:|---:|
| E0/B0 | .894/.581 | .940 | .959 | 43/77 | -.166 |
| E2/B2 | .961/.993 | .949 | .963 | 41/78 | -.141 |
| E1/B1 | .738/.840 | .922 | .933 | 43/77 | -.161 |

Fixed marker priors are identical across routes (marker-prior correlation 1)
and yield the strongest PIP/rank agreement. Learned E1/B1 plug-in marker priors
remain correlated (.916) but differ materially in scale: expected prior-active
counts from pooled alpha means are about 135 and 794. Posterior expected
occupancy is 144 versus 1,093. These plug-in prior summaries are descriptive
under nonconvergence and are not draw-wise prior-probability traces. The route
difference supports `LARGE-G4` in addition to the primary mixing decision.

## Scientific conclusion and next task

> At 5,000 individuals and 37,991 markers, general BayesR mixture sampling and
> heterogeneous marker-prior sampling with alpha fixed to truth have stable
> aggregate latent and variance diagnostics. Ordinary joint learned-alpha
> BayesRC/SBayesRC still fails multichain exploration, so increased information
> scale did not remove learned alpha-allocation feedback as the main
> computational bottleneck. SNP prioritization remains useful, but latent
> architecture is unresolved and the BED/block variance offset persists.

Replicated validation is not recommended because the one registered replicate
did not achieve full feasibility. The next task should be a focused,
posterior-preserving sampler audit of coordinated alpha/allocation movement,
while separately auditing the BED/block Vg and heritability contracts. It must
not reinterpret the current nonconverged learned posterior or launch the final
benchmark.
