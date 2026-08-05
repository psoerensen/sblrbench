# Study 06 annotation qualification result

> **Historical v1 record — unchanged scientific interpretation.** This page
> records `v1_sparse_qualification_failed`. The failed evidence remains valid
> development evidence and is not a completed benchmark. The separately
> versioned `v2_identifiable_qualification` has now run and failed independently;
> see `study06_v2_qualification_result.md`. The final benchmark remains
> unauthorized.

## Decision

**Qualification failed. The final 40-coordinate benchmark was not launched.**

The prespecified four-entry qualification was run on 2026-08-04 from
`02088ea8da42a5b9532d4ab843bb38456474e105` with `sblr` 0.2.0 at
`02e8c74baa906e83c4a08d42a9cc6339b4e81072` and qgdata at
`6cca5819e711d326cfb2614d7e9d9f34942612cd`. The decision artifact records
specification hash
`061ea41eb89da9a7bdee5a0dc221b8bb9321961e8bf0b9c3d9718827856357fb` and
`overall_decision = "failed"`.

No threshold, annotation, seed, prior, control, method, history length, or
candidate window was changed. Final mode was invoked only to audit the gate
and rejected the failed decision before any sampler call.

## Gate 0

- Package load succeeded.
- The deterministic repository suite passed: 589 tests, 0 failures, 0
  warnings, and 0 skips.
- Validation-only resolved exactly four qualification entries and 40 final
  coordinates, with zero sampler calls.
- The pinned package and data SHAs matched.
- True selected alpha and `sigmaSqAlpha` traces were available in the pinned
  BED and CSR interfaces.

The first completed fit exposed a framework descriptor mismatch: pinned
`sblr` uses `component_0_stick` through `component_2_stick`, whereas the
fixture expected `step_1` through `step_3`. The shared extractor and fixture
were corrected without changing the scientific estimand. The saved fit was
then reused by exact semantic hash. Draw-wise probit-stick reconstruction was
batched and tested against `sblr::sbayesrc_marker_pi()` to tolerance `1e-12`;
qualification omits only unused marker-component tables while retaining the
same registered draw-level non-null summaries.

## Fit and trace integrity

| Scenario | Method | Fit status | Recorded fit runtime (seconds) | Trace status |
|---|---|---:|---:|---|
| informative | BED BayesRC | completed | 2770.86 | 9,000 draws in each of four identifiable chains |
| informative | CSR SBayesRC | failed | unavailable | no valid retained history |
| uninformative | BED BayesRC | completed | 1938.42 | 9,000 draws in each of four identifiable chains |
| uninformative | CSR SBayesRC | failed | unavailable | no valid retained history |

Both CSR entries failed Gate 1 inside the pinned native backend:

- informative, chain 1: `sampleE_ST_operator: invalid projected residual scale`;
- uninformative, chain 0: `sampleE_ST_operator: invalid projected residual scale`.

The failed CSR calls do not expose a completed-fit runtime or reusable true
trace history. No final state, posterior mean, pooled vector, or compact
summary was substituted.

## Convergence gate

Neither completed BED history had a passing supported candidate. The
deterministic least-failing choice was the prespecified maximum supported view,
burn-in 3,000 and 6,000 retained draws.

| Scenario | Selected burn-in | Retained | Failed quantities | Maximum R-hat | Minimum bulk ESS | Minimum tail ESS | Maximum relative MCSE |
|---|---:|---:|---:|---:|---:|---:|---:|
| informative | 3000 | 6000 | 17 of 23 | 1.9601 | 5.5082 | 12.4688 | 0.4347 |
| uninformative | 3000 | 6000 | 21 of 23 | 1.6438 | 6.5572 | 11.0372 | 0.3968 |

The informative BED failures were:

- all four component-0 alpha coefficients;
- all four component-1 alpha coefficients;
- component-2 intercept and continuous-signal alpha coefficients;
- component-1 `sigmaSqAlpha`;
- effect variance and heritability;
- expected active count, enriched and unannotated mean non-null priors, and
  their contrast.

The uninformative BED failures were:

- all 12 alpha coefficients;
- component-1 and component-2 `sigmaSqAlpha`;
- effect variance, genetic variance, and heritability;
- expected active count, enriched and unannotated mean non-null priors, and
  their contrast.

The complete per-window and per-quantity values remain in the ignored local
evidence under `results/local/06_annotation_models/qualification/`:
`candidate_windows.csv`, `convergence.csv`, `fit_status.csv`, `runtime.csv`,
and `qualification_decision.json`. The two completed BED histories remain as
semantic checkpoints and must not be rerun when their identities match.

## Interpretation and next action

Failure is common across the annotation-aware routes but differs in form:
BED BayesRC completes and mixes far too poorly, while CSR SBayesRC fails a
native residual-scale integrity check in both scenarios. This is not evidence
about annotation-model performance. It is evidence that the prespecified
qualification gate correctly blocks the benchmark.

A focused `sblr` investigation is required before Study 06 can proceed:

1. reproduce and diagnose the CSR `sampleE_ST_operator` residual-scale failure
   for both semantic qualification identities;
2. investigate BayesRC alpha and group-scale mixing, especially intercept and
   early-stick coefficients, without removing the continuous annotation or
   weakening the gate;
3. validate any package-side correction with fixed-alpha reduction,
   directional-enrichment, seed, threading, and true-trace tests;
4. rerun the unchanged four-entry qualification at the same maximum history.

Study 06 remains **In development**. The `current-stop` capsule remains partial
evidence only; no `current` capsule was created, and no report or website
completion claim is authorized.
