# Study 06 final conclusion

## Study status

**CLOSED**

The authoritative machine-readable record is
[`final_decision.json`](../../results/reference/06_annotation_models/final_decision.json).
Historical qualification, audit, feasibility, and sampler-development decisions
retain their original labels; this document states the final scientific endpoint.

## Primary decision

**EST-R2 — probability/ranking functions are stable but annotation-contrast
magnitudes remain uncertain.**

The continuous-alpha decomposition is not a reliable primary biological
endpoint. `A alpha` is somewhat more stable but remains imperfect in some
directions. Annotation-induced SNP-prior rankings, SNP PIPs, posterior effects,
genetic values, and predictions are highly reproducible. Enriched and
continuous annotation directions are reproducibly positive, and the null
annotation is compatible with zero, but precise effect magnitudes are biased,
route-sensitive, and incompletely covered.

## Official qualifier

**EST-R5 — formal official multichain replication remains unavailable.**

The pinned official SBayesRC v0.2.6 trajectory at commit
`b95d3fcbad8ff358290922a58fff893439296138` is valuable descriptive evidence.
Its public seed argument does not control all native static/thread-local RNG
state, so fresh R processes cannot be shown to produce independent native
trajectories. The pinned trajectory shows qualitatively similar behavior, but
formal independently seeded multichain replication remains unavailable because
of this upstream native RNG contract. This is not a failed replication, and no
official multichain convergence claim is made.

## What was validated

- BayesR/SBayesR baseline controls and fixed-alpha BayesRC/SBayesRC controls.
- The retained block likelihood and GCTB-compatible residual-variance contract.
- Strong small-fixture official agreement for block residual variance, PIPs,
  effects, genetic values, prediction, and annotation benefit.
- Exactness of the audited hierarchy transitions. The final global reference is
  **PMA-R3 — exact but computationally impractical**; standard unrestricted
  continuous-alpha same-posterior sampler development is closed.
- Stable annotation-informed SNP prioritization and SNP-level inference on the
  frozen Study 06 fixture.

These results support the restrained statement that `sblr` implements the core
SBayesRC model and has been lightly but systematically validated against the
official implementation on this controlled simulation benchmark. They do not
establish universal GCTB replacement, posterior identity, or equivalence under
all architectures and scales.

## Raw first-stick annotation comparison

Stick-1 values use the calibrated truth directly registered by the final truth
mapping, official D1's single-trajectory mean, and the prespecified selected
qualification window for the learned `sblr` routes. The learned pooled means
are summaries under EST-R2, not converged truth estimates.

| Stick-1 coefficient | Truth | Official D1 | `sblr` BED | `sblr` block | Fixed-alpha oracle |
|---|---:|---:|---:|---:|---:|
| intercept | -1.653 | -2.406 | -2.934 | -2.201 | -1.653 |
| enriched | 1.600 | 1.507 | 1.897 | 1.642 | 1.600 |
| continuous | 0.300 | 0.586 | 0.638 | 0.622 | 0.300 |
| null | 0 | -0.156 | -0.042 | -0.371 | 0 |

First-stick enriched and continuous directions agree across truth, official
SBayesRC, and both learned `sblr` routes. The null coefficient is comparatively
small relative to the informative signals. Official and block enriched values
are close to truth; BED enriched is larger but directionally consistent; every
learned route overestimates the continuous raw coefficient. Quantitative
equality or exact calibration of raw alpha is not established, and raw magnitude
alone is not the biological endpoint.

## Why stick 1 matters

For a four-component stick-breaking model, with component 0 the null component,

\[
P(c_i>0\mid A_i)=q_{i1}=\Phi(A_i\alpha_1).
\]

Stick 1 therefore determines whether a SNP is active at all. Later sticks
redistribute SNPs that have continued past stick 1 among nonzero effect-size
components. Difficulty in later-stick alpha does not automatically imply an
unstable SNP PIP because the primary active-versus-null annotation signal is in
stick 1—and this is the stick showing consistent informative directions across
truth, official SBayesRC, BED, and block.

## Later-stick information loss

The frozen truth has component counts 1,329/84/50/37. The exact eligible counts
are 1,500 for stick 1, 171 for stick 2, and 87 for stick 3. The enriched binary
annotation occurs in 225/1,500 SNPs (15%); among the realized eligible sets it
occurs in 107/171 and 62/87 SNPs. Thus the dominant information loss here is the
sharp reduction in eligible observations and latent continuation outcomes—not
a claim that the later eligible sets contain only a handful of enriched SNPs.

Later-stick raw alpha values are not quantitatively interpretable with the same
confidence. Eligible sets shrink, component membership is latent, posterior
parameter covariance is strong, and alpha/allocation feedback weakens practical
exploration. In official D1, stick-2 and stick-3 intercept means are about 10.56
and 8.95 versus truths about -0.235 and -0.431; stick-3 intercept ESS is about
1.4. These are descriptive single-trajectory observations, not official
multichain diagnostics. The corresponding `sblr` later-stick directions also
show poor ESS, truth recovery, and joint exploration.

The annotation design itself is full rank (rank 4), with condition number about
2.23. There is no established annotation-design null space. The supported
language is weak practical determination, strong posterior dependence, and poor
exploration of some latent directions—not formal mathematical
non-identifiability. Posterior covariance is consistent with compensating
combinations such as a lower intercept with a larger enriched slope over the
realized annotation profiles.

## Annotation-variance comparison

This is a **descriptive scale comparison only**.

| Stick | Official D1 final `sigma_anno` | `sblr` BED mean `sigmaSqAlpha` | `sblr` block mean `sigmaSqAlpha` | Fixed-variance ablation |
|---|---:|---:|---:|---:|
| 1 | 1.976 | 2.137 | 2.044 | 1 |
| 2 | 1.043 | 1.897 | 1.540 | 1 |
| 3 | 4.803 | 3.136 | 1.788 | 1 |

Official `sigma_anno` final values and retained `sblr` `sigmaSqAlpha` posterior
summaries must not be assumed to be identical estimands without a dedicated
source-contract audit. Official retained variance histories are unavailable, so
the final values provide no posterior mean, interval, drift, ESS, or chain
parity. Fixed-variance values are ablations, not posterior estimates. D2's
secondary official final values are 3.480/0.384/0.413/0.605 under its native
five-component contract and are not part of the matched D1 comparison.

## Derived active annotation contrasts

The binary contrast is enriched 1 versus 0 with all other annotations held at
their observed values. Continuous and null contrasts are the frozen standardized
+1 SD versus -1 SD comparison. Brackets are 95% intervals.

| Annotation | Truth | `sblr` BED | `sblr` block | Official SBayesRC D1 |
|---|---:|---:|---:|---:|
| enriched binary | 0.423 | 0.192 [0.081, 0.397] | 0.299 [0.119, 0.642] | 0.204 [0.117, 0.302] |
| continuous signal | 0.0886 | 0.0555 [0.0239, 0.104] | 0.112 [0.0383, 0.236] | 0.077 [0.042, 0.118] |
| null annotation | 0 | -0.007 [-0.068, 0.046] | -0.068 [-0.268, 0.060] | -0.024 [-0.096, 0.024] |

Truth, official SBayesRC, BED, and block share the qualitative ordering

\[
\text{enriched} \gg \text{continuous} > \text{null}\approx0.
\]

Directions and relative ranking are convincing. Exact magnitudes remain
quantitatively uncertain and route-sensitive. In particular, BED's enriched
interval misses truth, while block and official intervals are broad or biased.
This is why the decision is EST-R2 rather than EST-R1.

The continuous annotation is the clearest demonstration. Its raw first-stick
coefficient is 0.300 in truth but 0.586 official, 0.638 BED, and 0.622 block.
After transformation, its active contrast is 0.0886 in truth versus 0.077,
0.0555, and 0.112. The probability-level functional is closer to truth and
qualitatively consistent across implementations, though not perfectly
calibrated. Likewise, the raw null values range from -0.371 to -0.042 while the
induced null contrasts are smaller and every interval includes zero. Isolated
raw coefficients should not be overinterpreted without the induced functional
and its uncertainty.

## Why SBayesRC can work despite difficult alpha parameters

For stick \(k\),

\[
\eta_{ik}=A_i\alpha_k,\qquad q_{ik}=\Phi(\eta_{ik}).
\]

The four component probabilities are

\[
\pi_{i0}=1-q_{i1},\quad
\pi_{i1}=q_{i1}(1-q_{i2}),\quad
\pi_{i2}=q_{i1}q_{i2}(1-q_{i3}),\quad
\pi_{i3}=q_{i1}q_{i2}q_{i3}.
\]

The inferential hierarchy is therefore

```text
alpha
→ A alpha
→ continuation probabilities q
→ SNP-specific mixture probabilities pi
→ latent component/effect states
→ SNP PIP / posterior beta
→ genetic values / prediction
```

Raw alpha values are upstream latent parameters. Different correlated
coefficient configurations can produce similar functions over the realized
annotation profiles. The probit transformation can further compress differences
when latent scores are in saturated regions: `Phi(4)` is about 0.99997, whereas
`Phi(7)` is effectively 1. This does not eliminate all alpha uncertainty; it
means weakly determined latent directions can have little effect on the
SNP-specific probabilities that drive inference.

Official D1 illustrates this compression descriptively. Half-specific
posterior-mean `A alpha` vectors have a minimum Pearson correlation near -0.975,
while the minimum nonconstant induced component-probability vector correlation
is at least 0.999. Raw-alpha ESS spans about 1.4–86.2 and the largest half-mean
shift is about 1.25 posterior SD. These halves belong to one native trajectory,
so they are descriptive stability evidence rather than independent chains.

Across MCMC samples,

\[
\alpha^{(s)}\rightarrow\pi_i^{(s)}\rightarrow c_i^{(s)},\beta_i^{(s)}.
\]

Marginal activity conceptually integrates over the annotation hierarchy,

\[
P(c_i>0\mid y)=E_{\alpha,\theta}\left[P(c_i>0\mid y,\alpha,\theta)\right],
\]

where \(\theta\) denotes the remaining parameters. Final SNP inference does not
normally condition on one plug-in alpha estimate. Posterior averaging propagates
hierarchy uncertainty, so noisy hyperparameters can coexist with stable marginal
SNP outputs. Poor MCMC exploration is not harmless in general; the empirical
result here is specifically that downstream quantities are much more stable
than the latent decomposition on this frozen fixture.

There is therefore no contradiction between difficult alpha inference and
useful SBayesRC SNP inference. Some later-stick combinations are weakly
determined and strongly correlated; the stick/probit map compresses many of
their functional consequences; posterior averaging then yields highly
reproducible SNP PIPs and effects even though every alpha coefficient is not
precisely recovered.

## Three levels of estimability

| Level | Examples | Study 06 result |
|---|---|---|
| Latent annotation decomposition | alpha, `sigmaSqAlpha`, later-stick combinations | Weak or difficult in some posterior directions, especially later sticks. |
| Annotation probability functions | `A alpha`, q, pi, active probability, counterfactual contrasts, SNP-prior ranking | Substantially more stable; strong directional and ranking information, with residual magnitude uncertainty. |
| Marginal SNP inference | SNP PIP, posterior beta, genetic value, prediction | Very stable across retained chains. |

In compact form:

```text
latent alpha decomposition        weak/correlated in some directions
        ↓
probability-level function        substantially more stable
        ↓
marginal SNP inference            highly stable
```

R-hat near one measures chain agreement, not truth recovery. The expected
active count makes this distinction explicit: BED averages about 351.8 and
block about 807.4 versus truth 180, with very large posterior variability,
despite R-hat near one. Thus `convergence != calibration`, and stable ranking
does not imply exact probability-magnitude recovery.

## SNP-level stability

| Quantity | BED minimum/range | Block minimum/range |
|---|---:|---:|
| PIP Pearson across chains | 0.994 | 0.969 |
| signed-beta Pearson | 0.9999 | 0.9991 |
| genetic-value chain correlation | 0.99994 | 0.99925 |
| active-prior Spearman | 0.986 | 0.950 |
| largest-component-prior Spearman | 0.999 | 0.995 |
| AUPRC | 0.591–0.597 | 0.514–0.560 |
| AUROC | 0.849–0.856 | 0.818–0.840 |
| genetic-value correlation to truth | 0.917–0.918 | 0.896–0.897 |
| phenotype prediction correlation | 0.654–0.655 | 0.631–0.633 |

Hierarchy-level uncertainty does not materially destabilize the primary
SNP-level outputs on the frozen Study 06 fixture.

## Scientific hierarchy of evidence

| Quantity | Role | Behavior | Recommended reporting |
|---|---|---|---|
| raw alpha | latent annotation coefficient | weak in some directions, especially later sticks | secondary plus diagnostics |
| `sigmaSqAlpha` / `sigma_anno` | hierarchy scale | limited cross-implementation comparability | diagnostic |
| `A alpha` | latent SNP score | improved but uneven | secondary |
| q / pi | SNP-specific prior architecture | high ranking stability | primary annotation-derived output |
| active annotation contrast | effect on non-null probability | stable direction, uncertain magnitude | primary annotation summary with uncertainty |
| SNP PIP | marginal SNP activity | highly stable | primary SNP output |
| posterior beta | marginal SNP effect | highly stable | primary SNP output |
| genetic value / prediction | downstream inference | highly stable | primary downstream output |

Within this controlled benchmark, saying SBayesRC “works” means it provides
annotation-informed prior differentiation, the correct qualitative annotation
ordering, useful probability functions and contrasts, improved SNP
prioritization, and stable SNP PIPs, effects, genetic values, and predictions.
It does not mean every alpha coefficient mixes rapidly, every annotation-effect
magnitude is calibrated, every variance parameter is well determined, or every
GCTB use case has been reproduced.

## Official comparison

On the controlled fixture, `sblr` reproduces the main behavior of official
SBayesRC: similar first-stick informative directions, broadly comparable active
annotation contrasts, highly concordant SNP inference, and the same qualitative
distinction between difficult latent alpha and much more stable induced
probabilities. This makes the hierarchy-level difficulty not obviously
`sblr`-specific.

The evidence remains delimited. Official D1 is a single trajectory; the
official and `sblr` prior/empty-stick contracts differ; exact posterior identity
is neither expected nor claimed; and `sigma_anno` is not proven equivalent to
`sigmaSqAlpha`. The emphasis on derived annotation probabilities and
enrichment-like quantities is consistent with official SBayesRC/GCTB output
practice, without implying a one-to-one audited output contract.

## Reporting recommendation

Future SBayesRC analyses should report as primary SNP outputs:

- SNP PIPs;
- posterior SNP effects;
- genetic-value and prediction performance where applicable.

Primary annotation outputs should be:

- annotation-induced active/non-null probabilities;
- component probabilities when scientifically relevant;
- counterfactual annotation contrasts with intervals and directional posterior
  probabilities;
- annotation-induced SNP rankings.

Raw alpha, `A alpha`, `sigmaSqAlpha`, posterior covariance/correlation, R-hat,
ESS, MCSE, and credible intervals remain useful secondary latent-model and
diagnostic outputs. Raw alpha can be informative, particularly for first-stick
direction, but should not be the sole or primary evidence for annotation
relevance.

## Closure and next methodology

Study 06 is complete. No further unrestricted continuous-alpha sampler
engineering or additional Study 06 scientific fit is recommended. The next
methodology is a separate Bayesian annotation-selection / annotation-PIP model
targeting `P(delta_j = 1 | y)`. It is motivated by the scientific estimand—not
merely computation—and must not be presented as a sampler fix for standard
SBayesRC. Its conceptual links to Bayesian MAGMA and Bayesian PoPS belong to
that future study; no such model is implemented here.

## Reproduction

Run `Rscript scripts/run_study06_estimability.R` from the repository root. The
retained frozen chains must already exist locally. This is offline posterior
analysis only: it invokes no sampler, regenerates no truth, and writes local
tables and figures under the ignored
`results/local/06_annotation_models/estimability_and_contrasts/` directory.
