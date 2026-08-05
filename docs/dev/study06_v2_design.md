# Study 06 v2 identifiable qualification design

> **Execution status (2026-08-05):** the frozen design below was run for its
> registered one-replicate qualification. All four fits completed, but the
> convergence, informative late-stick, and heritability route-agreement gates
> did not all pass. The final benchmark remains unauthorized. See
> [the v2 qualification result](study06_v2_qualification_result.md).

## Status and version boundary

Study 06 now has three explicit identities:

- `v1_sparse_qualification_failed` is the historical qualification attempt.
  It failed and is preserved as development evidence, not a completed
  benchmark.
- `v1_sparse_stress` preserves the same approximately 37,991-marker,
  50-active-marker architecture for sparse late-stick, prediction, PIP,
  active-count, variance, and numerical-integrity stress testing. It is not the
  primary qualification of precise late-stick annotation recovery.
- `v2_identifiable_qualification` is the new primary qualification design. It
  is designed and scaffolded but has not been run.

The final multi-replicate benchmark is not authorized. The tracked
`results/reference/06_annotation_models/current-stop/` tree remains byte-for-byte
unchanged. V2 outputs must contain the `v2_identifiable_qualification` path
component and cannot be written into or below `current-stop`.

## Why v2 is necessary

V1 combined about 37,991 markers with only 50 expected non-null markers and
active weights 0.60/0.30/0.10. The three active components therefore contained
only about 30, 15, and 5 markers in expectation. The later probit sticks had
very few eligible markers and binary outcomes, so precise annotation recovery
was weakly identified even if an implementation were correct. The completed
BED histories mixed poorly, while approximate sparse CSR failed with an
invalid projected residual scale. Those are valid stop results, but they do
not by themselves prove that ordinary BayesRC/SBayesRC is incorrect under an
identifiable generating design.

Package-side work established the proper probit-stick intercept prior and
correct prior-only annotation updates. It also identified retained low-rank
block eigen as the canonical scalable scalar summary-statistics route. Sparse
CSR failure was tied to approximate LD fidelity and residual scale, so sparse
CSR is not a v2 scientific gate. A mathematically exact two-marker allocation
move did not help because the v1 panel contained essentially no useful high-LD
pair graph; that production move was rejected. A larger collapsed block move
has not been designed and is deferred rather than presumed necessary.

V2 asks the narrower next question: does ordinary BayesRC/SBayesRC behave
scientifically when component membership and annotation effects carry enough
information to be identified?

## Prespecified v2 data design

The design retains the deterministic QGG chromosome-1 source and QC rules:

| Quantity | Value |
|---|---:|
| Total samples | 2,000 |
| Training/test samples | 1,400 / 600 |
| Split seed | 3,101 |
| Retained markers | 1,500 |
| Blocks | 15 |
| Markers per block | 100 |
| Target heritability | 0.50 |
| Target non-null markers | 180 |
| Mixture multipliers | 0, 0.01, 0.1, 1 |
| Active weights | 0.50, 0.30, 0.20 |
| Minimum realized active counts | 60, 30, 20 |
| Retained eigen mass | 0.995 |

After QC, markers are ordered by physical position with marker ID as the
deterministic tie-breaker. The ordered QC panel is divided into 15 equal
genomic strata. The centered 100 consecutive QC markers in each stratum form a
block. Selection uses no causal status, annotation, phenotype, GWAS statistic,
or outcome. The resulting block starts are 1, 101, ..., 1401; BED and block
eigen receive the identical 1,500 marker IDs in identical order.

The deterministic audit reports per-block count, endpoints, physical span,
training-genotype rank, minimum and maximum positive correlation eigenvalues,
and the retained count/mass at 0.995. It also reports complete cross-block
maximum absolute training correlation and method-input marker equality. Exact
zero cross-block LD is not required: the summary route deliberately uses the
prespecified block factorization.

On the cached pinned QGG panel, all 15 blocks contain 100 markers and have rank
100. Physical spans range from 225,241 to 4,469,592 bp; minimum positive
eigenvalues range from 0.5412 to 0.5724 and maxima from 1.5489 to 1.5968. At
0.995 each block retains all 100 positive modes for this panel. The complete
maximum cross-block absolute training correlation is 0.1334. Method marker
identity and order equality passes.

## Genetic and annotation truth

The marginal target is
`c(0.88, 0.06, 0.036, 0.024)`, giving expected active counts 90/54/36.
Effects are sampled with the package's BayesR/BayesRC component semantics and
scaled through the shared simulation framework. Training genetic variance is
scaled to `h2/(1-h2)` with unit residual variance. Target and realized raw
genetic, genetic, residual, phenotype, and heritability quantities are stored.

The ordered annotation matrix is:

1. exact-one intercept;
2. 15% `enriched_binary` membership selected independently of causal truth;
3. sample-standardized `continuous_signal`;
4. sample-standardized `null_annotation`.

For the informative scenario, the prespecified non-intercept coefficient rows
across the three sticks are:

| Annotation | Stick 1 | Stick 2 | Stick 3 |
|---|---:|---:|---:|
| enriched binary | 1.60 | 0.30 | 0.20 |
| continuous signal | 0.30 | 0.15 | 0.10 |
| null annotation | 0 | 0 | 0 |

Only the three intercepts are numerically calibrated, jointly, to the marginal
target. For the deterministic panel they are approximately -1.6528, -0.2349,
and -0.4313. The informative design places 59.79% of expected non-null
probability in the enriched 15%; mean expected non-null probabilities are
0.4783 enriched and 0.05677 unenriched. The continuous signal is directional
but moderate, and the null truth is exactly zero. Enrichment is not oracle:
causal membership is sampled afterward, and validation requires both causal
unenriched and non-causal enriched markers.

The uninformative scenario uses the same marker panel and the same annotation
columns, sets every non-intercept coefficient to exactly zero, and reverse-
calibrates its intercepts to the informative marginal probabilities. The two
scenario marginals match to numerical tolerance; they differ in annotation
information, not sparsity or intended genetic variance. The maximum absolute
correlation among non-intercept annotations is 0.0202 and the design rank is 4.

## Truth-identifiability gates

Simulation fails clearly rather than altering samples or fitted values when:

- realized total non-null count is outside 160--200;
- active component counts are below 60/30/20;
- stick eligible counts are below 500/120/50;
- either stick outcome count is below 100/25/20;
- any required enriched-binary by stick-outcome cell is empty;
- the four-column annotation matrix loses rank on an eligible subset; or
- realized heritability is non-finite or differs from 0.50 by more than 0.02.

With the registered replicate-1 seeds, informative truth realizes
`1329/84/50/37` component counts (171 active) and uninformative truth realizes
`1334/81/51/34` (166 active). Stick eligible/outcome-one counts are
1500/171, 171/87, 87/37 and 1500/166, 166/85, 85/34, respectively. The smallest
required binary-annotation cell is 12 informative and 5 uninformative. Both
realized component-variance heritabilities are 0.50.

## Public method registry

The exact public sibling `sblr` routes resolved at
`f2e3647920ed7e8b1ea9d47a6571b3753285682a` are:

1. `sblr::stblr_bed(..., method = "bayesr")`;
2. `sblr::stblr_bed(..., method = "bayesrc", annotation = ...)`;
3. `sblr::stblr_block_eigen(..., method = "sbayesr",
   representation = "low_rank", eigen_policy =
   "cumulative_positive_mass", eigen_prop = 0.995)`;
4. the same block-eigen call with `method = "sbayesrc"` and `annotation`.

The package's current proper intercept prior default is used; v2 does not pass
`intercept_flat`. No v2 route enables sparse approximate CSR, pair allocation,
or an unimplemented collapsed block-allocation move. Complete CSR is outside
the qualification gate and may only be a separate small engineering fixture.

The optional `fixed_true_annotation_coefficients` diagnostic is disabled and,
when explicitly requested after a failed qualification, initializes at truth
and uses the public `updateAlpha = FALSE` control. Fixed true marker-specific
component probabilities are registered as unavailable because neither public
BED nor block-eigen API accepts them; no workaround is added.

## Estimands and unavailable quantities

The shared extraction supports genetic, residual, and base-effect variance;
heritability; component probabilities; expected active count; chain component
occupancies; posterior SNP effects; held-out genetic values and predictions;
causal ranking; PIP AUROC/AUPRC; annotation coefficient and variance traces;
stick counts; draw-wise marker priors and expected non-null count; enriched,
continuous, and null prior contrasts; and component-specific occupancy.

True retained traces remain mandatory for trace estimands. A final state or a
posterior-mean coefficient is never substituted; unavailable quantities retain
the explicit `unavailable` status.

## Qualification and decision rule

The first v2 qualification is one replicate in each annotation scenario for
BED BayesRC and retained block-eigen SBayesRC, with four registered chains.
Maximum history and candidate windows remain 9,000 iterations, burn-ins
1,000/2,000/3,000, and retained windows 2,000/4,000/6,000. Every required
quantity must meet R-hat <= 1.01, bulk ESS >= 400, tail ESS >= 400, and relative
MCSE <= 0.05. A failure leaves the final benchmark blocked.

Scientific gates are interval- and direction-based, not exact truth recovery.
Directional annotation and prior-contrast claims require posterior direction
probability at least 0.90; null compatibility uses a 95% interval containing
zero. Relative to BayesR/SBayesR, a ranking or prediction-correlation decrease
greater than 0.02 is material. Following Study 05 readiness philosophy, BED
versus retained block eigen allows prediction and heritability differences up
to 0.05 and requires marker-effect correlation at least 0.95. These are
qualification tolerances, not claims that the likelihood representations are
identical.

The exact future launch command, after the isolated R library contains the
pinned sibling build, is:

```powershell
Rscript scripts/run_benchmark.R --study 06_annotation_models --profile benchmark --output-dir results/local/06_annotation_models/v2_identifiable_qualification --resume true --validate-only false --mode qualification
```

Qualification-length execution is never triggered by unit tests or the
default analysis mode. Final mode additionally requires an identity-matched,
passing v2 decision artifact.

## What would justify sampler redesign

A v2 failure after valid truth gates, exact marker alignment, valid operator
scale, retained traces, and the full registered convergence history would
justify the two disabled diagnostic-isolation checks. Failure with fixed true
coefficients would implicate allocation/effect-likelihood or route behavior;
success there but failure under learned coefficients would implicate
annotation-update feedback. Only evidence that remains after this isolation
would justify revisiting larger collapsed allocation moves. V1 sparsity alone,
or failure of approximate sparse CSR, does not justify that conclusion.
