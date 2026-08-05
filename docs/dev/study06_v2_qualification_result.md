# Study 06 v2 identifiable qualification result

## Decision

**Failed — convergence or mixing, with additional scientific and route-gate
failures. The final benchmark is not authorized and was not launched.**

The versioned `v2_identifiable_qualification` was run on 2026-08-05. All four
registered fits completed at the maximum 9,000-iteration history, but no entry
met the selected-window convergence contract. The informative scenario also
failed the prespecified late-stick direction checks on both routes, and the
BED/block-eigen heritability differences exceeded 0.05 in both scenarios.
This is not a posterior-effect-correlation-only failure and therefore is not a
`criterion review required` result.

```text
v1 sparse qualification: failed and preserved
v2 identifiable qualification: failed with convergence/mixing, scientific-recovery, and route-agreement blockers
final benchmark: not authorized
```

This result does not establish method superiority or final Study 06
validation.

## Frozen identity

| Item | Value |
|---|---|
| `sblrbench` source HEAD | `f8fceaf337087f589d12c69fae113d5211959a57` |
| `sblr` source and installed `RemoteSha` | `f2e3647920ed7e8b1ea9d47a6571b3753285682a` |
| installed `sblr` version | `0.2.0` |
| qgdata commit | `6cca5819e711d326cfb2614d7e9d9f34942612cd` |
| specification hash | `241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56` |
| study version | `v2_identifiable_qualification` |
| selected window for every entry | burn-in 3,000; retained 6,000 |
| chains | four registered chains per entry |

The isolated package library was
`results/local/current_benchmark_refresh/rlib`. A fresh R process resolved
`sblr` there, not the system `0.1.2` installation and not the older isolated
`02e8c74...` build.

The exact qualification command was:

```powershell
Rscript scripts/run_benchmark.R --study 06_annotation_models --profile benchmark --output-dir results/local/06_annotation_models/v2_identifiable_qualification --resume true --validate-only false --mode qualification
```

`SBLR_BENCH_GLIST` pointed to the repository's cached canonical QGG Glist. It
changes no scientific identity; the qgdata file identities remain those in the
frozen specification.

## Deterministic preflight

The preflight reproduced 1,500 markers in 15 separated blocks of 100, a
disjoint 1,400/600 train/test split, exact intercept ones, and identical marker
identity/order for every method input. The complete cross-block correlation
audit gave maximum absolute correlation 0.1333644. Every block had rank 100;
the minimum positive eigenvalue ranged from 0.54118 to 0.57237, the maximum
from 1.54890 to 1.59676, and `eigen_prop = 0.995` retained 100/100 modes with
positive-mass fraction 1 in every block.

Informative and uninformative marginal component probabilities matched to
`1.11e-16`: 0.8800001/0.05999996/0.03599999/0.02399998, or approximately
180 expected active markers with expected active counts 90/54/36. Realized
counts were 1329/84/50/37 and 1334/81/51/34, respectively; realized
heritability was exactly 0.50 under the registered calculation. All truth
identifiability gates passed.

The block-eigen result tests the retained-factor execution route, block
construction/order, summary-statistics model semantics, blockwise likelihood
factorization, and SBayesRC annotation behavior. Because every block retained
all 100 modes, it does **not** test substantial eigenspace rank truncation.

## Execution record and minimal runner fixes

The first launch completed the informative BED history, then exposed a
benchmark extraction defect: `.annotation_parameter_estimates()` shadowed its
simulation bundle with a convergence-trace object. The fix only renamed the
local trace variable and added a focused regression test.

The identity-matched resume reused that BED history, completed the
uninformative BED history, and exposed a second public-API wiring defect. The
pinned public block-eigen SBayesRC route supports component multipliers as
`gamma`, while the benchmark forwarded its semantic `mixture_var` control.
The runner now translates the unchanged vector `c(0, 0.01, 0.1, 1)` to
`gamma` only for public block-eigen SBayesRC. A focused regression test verifies
that BED retains `mixture_var` and block eigen receives the identical values as
`gamma`.

The final identity-matched resume reused both BED histories and completed both
block-eigen histories. No source in the sibling `sblr` repository was changed.
No seed, iteration count, prior, threshold, annotation, marker, block, or
operator setting changed. There were no sampler retries with changed identity.

| Scenario | Route | Reused in final continuation | Runtime (s) | Checkpoint hash |
|---|---|---:|---:|---|
| informative | BED BayesRC | yes | 129.24 | `76f64e9d960e2e3574a05ad0e17f05721154770cdaa561a24b341b5c614bd29a` |
| informative | block-eigen SBayesRC | no | 54.74 | `99225899e89b1e82a61114fa7b49ae24fe5bfec8a5ddc24eafb94c1475fac828` |
| uninformative | BED BayesRC | yes | 109.25 | `bf989a85f7ee67b8a3159d605b3f82671dd3fa69acfc569f8fdeea08456e7afd` |
| uninformative | block-eigen SBayesRC | no | 60.41 | `f17f73c14d42a97ef4930536fd912f06df322fb696c46ede4498a19cf307337e` |

The package reported analytical memory upper bounds of 0.0320 GiB for BED and
0.0197 GiB for block eigen; it did not measure peak RSS. Read-only process
snapshots observed approximately 1.0 GB resident memory near the largest point,
so the analytical figures must not be presented as measured peaks.

## Convergence result

All entries exhausted the maximum supported window. Thresholds remained R-hat
at most 1.01, bulk and tail ESS at least 400, and relative MCSE at most 0.05.

| Scenario | Route | Failed/monitored | Max R-hat | Min bulk ESS | Min tail ESS | Max relative MCSE |
|---|---|---:|---:|---:|---:|---:|
| informative | BED BayesRC | 19/23 | 1.0844 | 34.65 | 32.98 | 0.1770 |
| informative | block-eigen SBayesRC | 22/23 | 1.1607 | 20.86 | 46.66 | 0.2221 |
| uninformative | BED BayesRC | 14/23 | 1.0200 | 121.77 | 117.64 | 0.0965 |
| uninformative | block-eigen SBayesRC | 19/23 | 1.1035 | 33.11 | 36.94 | 0.1771 |

Across entries, 74 of 92 selected quantity checks failed. Failures were
concentrated in alpha, `sigmaSqAlpha`, prior expected-active/contrast summaries,
and effect variance. BED genetic variance, residual variance, and heritability
were chain-consistent; block-eigen residual variance and heritability failed in
both scenarios. This is not merely broad but chain-consistent uncertainty:
many alpha and prior summaries have chain-specific locations, low ESS, high
R-hat, and high MCSE. All retained numerical values were finite; the final
record is not a numerical-failure classification.

Selected-marker component occupancy traces were not retained by the public
fits. The record therefore preserves them as explicitly unavailable and does
not substitute final states for convergence traces. For diagnosis only, local
tables report final per-chain component counts and posterior expected counts.
Uninformative block-eigen posterior expected active counts by chain were
approximately 155, 254, 137, and 140, which is direct evidence of chain-location
disagreement. Informative BED prior expected-active means were approximately
59, 59, 57, and 54; informative block-eigen means were 159, 123, 113, and 141.

## Annotation and prior recovery

In the informative scenario, both routes recovered the first-stick direction
but failed the prespecified 0.90 direction probability for the later sticks.

| Route | Annotation | Stick 0 | Stick 1 | Stick 2 |
|---|---|---:|---:|---:|
| BED | enriched binary | 1.000 | 0.715 | 0.632 |
| BED | continuous signal | 1.000 | 0.725 | 0.665 |
| block eigen | enriched binary | 1.000 | 0.731 | 0.610 |
| block eigen | continuous signal | 0.999 | 0.730 | 0.656 |

All informative null-annotation 95% intervals contained zero. Both enriched
prior contrasts had positive-direction probability 1.000; continuous-prior
contrast probabilities were 1.000 (BED) and 0.999 (block eigen). Thus the
aggregate prior contrast moved in the intended direction even while the late
sticks failed direction and convergence gates.

For the uninformative scenario, all non-intercept coefficient intervals and
both prior-contrast intervals contained zero on both routes. The committed
zero-compatibility scientific checks all passed. However, those intervals are
not treated as reliable recovery claims because the corresponding chains did
not satisfy the convergence contract.

`sigmaSqAlpha` selected-window means for sticks 0/1/2 were approximately
2.14/1.90/3.14 (informative BED), 2.04/1.54/1.79 (informative block eigen),
0.73/2.13/2.95 (uninformative BED), and 0.78/1.33/3.64 (uninformative block
eigen), with broad, heavy-tailed uncertainty and multiple diagnostic failures.

## Prediction, ranking, and route checks

| Scenario | Route | h2 | Genetic-value correlation | Phenotype prediction correlation | PIP AUROC | PIP AUPRC |
|---|---|---:|---:|---:|---:|---:|
| informative | BED | 0.4162 | 0.9178 | 0.6545 | 0.8528 | 0.5948 |
| informative | block eigen | 0.5074 | 0.8962 | 0.6324 | 0.8329 | 0.5503 |
| uninformative | BED | 0.3954 | 0.8845 | 0.6688 | 0.6208 | 0.3413 |
| uninformative | block eigen | 0.4825 | 0.9012 | 0.6863 | 0.6072 | 0.2949 |

| Scenario | Absolute h2 difference | Absolute genetic-correlation difference | Raw effect correlation |
|---|---:|---:|---:|
| informative | 0.0912 (fail) | 0.0216 (pass) | 0.9579 (pass) |
| uninformative | 0.0871 (fail) | 0.0167 (pass) | 0.9681 (pass) |

Direct correlation between BED and block-eigen validation genetic values was
0.9618 informative and 0.9631 uninformative. Median blockwise effect
correlations were 0.9630 and 0.9382, but minima were 0.3982 and 0.6535.
Effect discrepancies were not concentrated at markers with higher within-block
LD: correlation between absolute route discrepancy and maximum within-block
absolute LD was only 0.013 and 0.022, and the top discrepancy decile had nearly
the same mean maximum LD as the remaining markers.

Raw effect correlation is not the sole or decisive blocker. Predictions agree
strongly despite some block-specific redistribution of SNP effects, but the h2
route gates, convergence gates, and informative late-stick scientific gates
independently fail.

The qualification did not include non-annotation BayesR/SBayesR baselines, so
relative annotation advantage over those methods is unavailable rather than
inferred. The absolute ranking and prediction values above are diagnostic
context, not a substitute comparison and not a method-superiority claim.

## Evidence location and next action

The committed framework declares the decision beneath the versioned local
output root but no tracked v2 reference-capsule destination. Consequently, raw
histories and the exact decision remain under:

```text
results/local/06_annotation_models/v2_identifiable_qualification/
```

The compact local record includes the decision JSON, fit status, all candidate
and selected convergence tables, scientific and route checks, block and truth
audits, runtimes, checkpoint SHA-256 inventory, chain summaries, warnings, and
supplemental route context. This document is the smallest tracked result note;
the historical v1 `current-stop` capsule remains byte-for-byte unchanged.

Post-run validation parsed every changed R file and passed all 92 focused
Study 06 expectations. The repository-wide lightweight suite recorded 607
passes, one failure, and six errors; all seven problems are older studies or a
shared helper expecting the previous isolated `sblr` SHA `02e8c74...`, whereas
this qualification correctly requires `f2e364...`. No Study 06 test failed.
`git diff --check` passed.

The Study 06 report and studies index both completed document execution and
standalone HTML rendering. Project-level resource copying also attempted to
walk an unrelated stale broken path under ignored
`results/local/study06-v2-doc-validation/` and failed after HTML conversion;
that local tree was not modified. The standalone rendered outputs were
validation artifacts and were not retained.

The decision artifact SHA-256 is
`80e35cf8a03fc6d42659727cd5104e1788ddc47cc0503a38f241f240e4c56832`.
The final package fit-warning table is empty. The runner emitted four benign
data-frame row-name warnings during table binding; R also noted that the local
`testthat` binary was built under R 4.4.3 while execution used R 4.4.1.

Before any final benchmark can be reconsidered, the failed convergence/mixing
and late-stick recovery must be reviewed in a separate task. The prespecified
disabled diagnostic-isolation profile with fixed true annotation coefficients
is available for that future investigation. No collapsed block move, pair
allocation, sparse CSR gate, altered threshold, or longer chain is authorized
by this result.
