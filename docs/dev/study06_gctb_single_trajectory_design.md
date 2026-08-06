# Study 06 official SBayesRC single-trajectory design

## Scope

This diagnostic asks whether one official SBayesRC v0.2.6 trajectory recovers
approximately the same SNP-level signal as the existing Study 06 `sblr`
results. It is a descriptive implementation comparison, not a convergence
comparison, formal posterior validation, qualification rerun, or final
benchmark.

The prior parity decision is immutable:

```text
GCTB-P5: multichain parity is blocked because official SBayesRC v0.2.6
does not connect its public seed argument to the native RNG streams.
```

Requested seeds are recorded but do not seed native streams. Every condition
runs in a fresh child process with one OpenMP and one BLAS thread and therefore
uses the native default stream. D0/D1 are reproducible under that stream; D2's
wrapper tuning additionally calls unseeded R `rnorm()`, so its pseudo-summary
tuning input is not fixed by the ignored public seed. Duplicate processes must
never be described as independent chains. R-hat is not calculated.

## Immutable inputs

The diagnostic reuses, without regeneration, the validated export under
`results/local/06_annotation_models/gctb_parity/export/`.

| Identity | SHA-256 |
|---|---|
| Study specification | `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56` |
| informative truth | `169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb` |
| marker order | `135c3604e0c4395349475b8126e0957db265c4b348a2e339daa3e7ddf2316a29` |
| GWAS | `1bd0abb220c4e7f9ca58ed11f2c2913e6e142868ba17820b988670cd01cd4610` |
| annotation | `b2442b5c074b1cd6fa0eb047fcedae6d3296a15d2b5f985c2ac288b534a3b156` |
| LD blocks | `369df6bbed513da7913f960b9f79e35d96cdeac321cf3105b303c932d8e2c0a4` |

All 1,500 markers have N=1,400, exact marker and effect-allele order, no
unresolved flips, and 15 blocks of 100 markers. Matched conditions retain all
100 positive modes per block. The annotation file contains the official
wrapper's intercept slot followed by enriched binary, continuous informative,
and null columns.

Official SBayesRC is version 0.2.6 from clean source
`b95d3fcbad8ff358290922a58fff893439296138`, loaded only from the existing
ignored isolated library. The installed DESCRIPTION does not expose
`RemoteSha`; the clean checkout and the existing validated export/smoke
records supply the pin. Neither official source nor sibling `sblr` may change.

## Registry

All conditions use 9,000 iterations, 3,000 burn-in iterations, 6,000 retained
draws, `outFreq=1`, detailed output, sparse beta histories, `allMixVe`, `nbsq`,
`starth2=0.5`, and one thread.

| ID | Model | Gamma | Start pi | Annotations | LD/tuning | nominal seed |
|---|---|---|---|---|---|---:|
| D0 | matched official SBayesR | 0/.01/.1/1 | .88/.06/.036/.024 | none | all positive modes; no tuning | 711121 |
| D1 | matched official SBayesRC | 0/.01/.1/1 | .88/.06/.036/.024 | informative three-column design plus official intercept | all positive modes; no tuning | 721121 |
| D2 | native official SBayesRC | 0/.001/.01/.1/1 | .990/.005/.003/.001/.001 | same annotations | tune 150/100 over .995/.99/.95/.9; prior tuning off | 731121 |

## Sequential stop rules

D0 runs first. Its initial native histories must exactly equal the overlapping
states in the prior fresh-process D0 smoke. D1 starts only after D0 produces
finite, aligned output and passes that equality check. D2 starts only after D1
passes the same structural checks. A numerical failure, malformed output,
marker or allele mismatch, non-finite SNP output, changed source checkout, or
non-reproducible default stream stops later conditions.

## Extraction and mapping

Native RDS arrays, not the wrapper's recycled four-component labels, define
component order. PIP is `1 - null probability`. D1 has sticks p1/p2/p3 and D2
has p1/p2/p3/p4; every alpha file is ordered Intercept, enriched binary,
continuous informative, null. Official p1 is the probability of continuing
out of the null component. Retained `sigmaSqAlpha` history is unavailable and
must not be inferred from final values.

## Descriptive analysis

For every scalar and low-dimensional history, first/middle/final thirds,
halves, running/cumulative means, autocorrelation, single-chain ESS, trend,
and standardized drift are summarized. Drift is labelled little below 0.25
posterior SD, moderate from 0.25 to below 0.75, and substantial at 0.75 or
more. These labels cannot establish convergence or exclude multimodality.

SNP comparisons require exact identifier/order/effect-allele agreement and
use PIP correlation/calibration, ranking overlap at 10/25/50/100/200,
Bayesian FDR, effect agreement, causal truth, and validation prediction. C0-C5
are the comparisons registered in the task. Learned non-converged `sblr`
results are reported both pooled and separately for each chain; pooled values
are not called converged posterior estimates.

Raw official outputs and compact analysis tables stay ignored under
`results/local/06_annotation_models/gctb_single_trajectory/`. Only this design,
the result note, and the compact decision JSON are tracked.
