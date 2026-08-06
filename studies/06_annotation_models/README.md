# Study 06 — Annotation-informed models

**Status: In development — formal qualifications failed; final benchmark not
authorized.**

## Current status

```text
v1 sparse qualification: failed and preserved
v2 identifiable qualification: failed
paired power isolation: completed
package-side hierarchy and transition audits: completed
official SBayesRC multichain parity: blocked by the v0.2.6 seed contract
official SBayesRC single-trajectory descriptive comparison: completed
final benchmark: not authorized
large n=5000, m=37,991 feasibility experiment: technically blocked (LARGE-F6); no scientific fits run
```

Study 06 asks whether correctly specified annotations improve causal-marker
prioritization and signal recovery in BayesRC and SBayesRC, and whether the
annotation hierarchy, mixture allocation, and individual-level versus
summary-statistic likelihood routes behave reliably. The immutable v2 identity
is:

```text
Specification:      241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56
Informative truth:  169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb
```

The [main report](report.qmd) is the primary readable account. The
[evidence synthesis](../../docs/dev/study06_annotation_inference_evidence_synthesis.md)
is the authoritative chronological technical history. Experiment-specific
result documents and decision JSON files remain authoritative for their own
formal or descriptive decisions.

## Current conclusions

### Established

- BayesR/SBayesR baselines and fixed-true-alpha BayesRC/SBayesRC converge.
- The independently audited probit-stick and `sigmaSqAlpha` conditionals are
  correct, and the hierarchy converges when allocations are frozen.
- Fixing `sigmaSqAlpha` or using the production-equivalent variance prior does
  not resolve dynamic learned-alpha mixing.
- Repeating hierarchy updates improves BED alpha marginals, but not the full
  joint posterior; repeating allocation sweeps is unfavorable.
- The tested complete-state coupling ladder is algebraically correct on a tiny
  state space but has effectively zero overlap on Study 06. Existing retained
  state cannot evaluate exact partial exchange.
- The BED/block-eigen heritability offset persists without annotations, with
  alpha fixed, and with all block modes retained.
- Informative annotations improve causal-marker ranking. One matched official
  SBayesRC trajectory independently reproduces that benefit and strongly agrees
  with `sblr` for SNP effects, PIPs, ranks, and validation genetic values.

### Descriptive but useful

Learned-alpha AUPRC and prediction summaries are descriptive under
non-convergence. Official alpha, occupancy, architecture, native five-component
results, runtime, and drift are single-trajectory descriptive. Pooled learned
BED/block summaries and per-chain official-versus-`sblr` comparisons are not
formal posterior validation.

### Unresolved

Official multichain convergence, quantitative alpha and annotation-variance
parity, latent architecture parity, active/component-count semantics,
residual/effect-scale contracts, the BED/block variance calibration, and
whether larger information scale resolves joint mixing remain unresolved.

## Authoritative documents

- [Main Study 06 report](report.qmd) — concise scientific overview.
- [Technical evidence synthesis](../../docs/dev/study06_annotation_inference_evidence_synthesis.md)
  — authoritative chronological ledger and evidence map.
- [Documentation inventory](../../docs/dev/study06_documentation_inventory.md)
  — tracked benchmark and package-side evidence map.
- [Documentation cleanup result](../../docs/dev/study06_documentation_cleanup_result.md)
  and [decision](../../docs/dev/study06_documentation_cleanup_decision.json).

Formal decisions remain separate:

- [v1 failed qualification](../../docs/dev/study06_annotation_qualification_result.md);
- [v2 failed qualification](../../docs/dev/study06_v2_qualification_result.md);
- [paired power isolation](../../docs/dev/study06_v2_power_isolation_result.md)
  and [decision](../../docs/dev/study06_v2_power_isolation_decision.json);
- [official multichain blocker](../../docs/dev/study06_gctb_parity_result.md)
  and [GCTB-P5 decision](../../docs/dev/study06_gctb_parity_decision.json);
- [official single-trajectory result](../../docs/dev/study06_gctb_single_trajectory_result.md)
  and [GCTB-D decision](../../docs/dev/study06_gctb_single_trajectory_decision.json).
- [large information-scale design](../../docs/dev/study06_large_feasibility_design.md),
  [blocked result](../../docs/dev/study06_large_feasibility_result.md), and
  [LARGE-F6 decision](../../docs/dev/study06_large_feasibility_decision.json).

## Recommended reading paths

### Quick scientific overview

1. [Main report](report.qmd).
2. [Evidence synthesis](../../docs/dev/study06_annotation_inference_evidence_synthesis.md).
3. [Official single-trajectory result](../../docs/dev/study06_gctb_single_trajectory_result.md).

### Formal qualification history

1. [v1 qualification result](../../docs/dev/study06_annotation_qualification_result.md).
2. [v2 design](../../docs/dev/study06_v2_design.md).
3. [v2 qualification result](../../docs/dev/study06_v2_qualification_result.md).
4. [paired power-isolation result](../../docs/dev/study06_v2_power_isolation_result.md).

### Implementation and sampler-diagnostic history

1. [`sblr` alpha-hierarchy audit](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_alpha_hierarchy_joint_sampling_audit.md).
2. Component-trace correction in that audit.
3. [`sblr` kernel-composition audit](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_allocation_hierarchy_kernel_composition.md).
4. [`sblr` coupling-tempering screen](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_bed_coupling_tempering_screen.md).
5. [`sblr` partial-exchange feasibility audit](https://github.com/psoerensen/sblr/blob/master/docs/dev/study06_partial_exchange_feasibility.md).

### Official SBayesRC comparison

1. [Parity design](../../docs/dev/study06_gctb_parity_design.md).
2. [GCTB-P5 seed-contract blocker](../../docs/dev/study06_gctb_parity_result.md).
3. [Single-trajectory design](../../docs/dev/study06_gctb_single_trajectory_design.md).
4. [Single-trajectory result](../../docs/dev/study06_gctb_single_trajectory_result.md).

### Reproducibility and evidence

- Immutable [v1 current-stop capsule](../../results/reference/06_annotation_models/current-stop/README.md),
  [manifest](../../results/reference/06_annotation_models/current-stop/benchmark_manifest.json),
  and [checksums](../../results/reference/06_annotation_models/current-stop/checksums.csv).
- Versioned [specification](spec.R), [analysis](analysis.R), and
  [annotation design](annotation-design.R).
- Diagnostic helpers: [paired isolation](power-isolation.R),
  [official export/parity](gctb-parity.R), and
  [single trajectory](gctb-single-trajectory.R), plus the frozen
  [large-feasibility profile](large-feasibility.R).
- Entry points: [benchmark runner](../../scripts/run_benchmark.R),
  [paired isolation](../../scripts/run_study06_power_isolation.R),
  [official parity](../../scripts/run_study06_gctb_parity.R), and
  [single trajectory](../../scripts/run_study06_gctb_single_trajectory.R), plus
  the [large-feasibility runner](../../scripts/run_study06_large_feasibility.R).

## Historical evidence and simulation designs

The 37,991-marker, approximately 50-active-marker v1 design is preserved as
failed qualification evidence and as the `v1_sparse_stress` concept. It is a
sparse late-stick/numerical stress scenario, not proof that annotation sampling
is incorrect and not a completed benchmark.

The v2 identifiable design uses 2,000 people (1,400 training, 600 validation),
1,500 markers in 15 separated blocks of 100, 171 realized active markers in the
informative truth (component counts 1,329/84/50/37), heritability 0.50, and
`gamma = c(0, 0.01, 0.1, 1)`. BED and block-eigen routes share marker order,
truth, phenotype, annotations, and priors. Every block retains 100/100 positive
modes, so the block route tests factorization and summary-model semantics, not
substantial eigen truncation.

## Next phases

1. **Documentation cleanup** — this consolidation.
2. **Large information-scale feasibility experiment** — designed and frozen,
   but [technically blocked](../../docs/dev/study06_large_feasibility_result.md)
   before scientific fitting (`LARGE-F6`). The current public API lacks compact
   genome-wide occupancy-count traces, and the required full-positive-mode B0
   smoke rejected its iteration-0 residual scale.
3. **Package-side contract audit** — add/test compact occupancy aggregates and
   diagnose the frozen B0 residual-scale rejection without changing this truth.
4. **Exact feasibility execution and replicated validation** — only after those
   contracts are resolved; replicated validation still requires the exact
   single replicate to succeed first.

An optional fixed-state official-versus-`sblr` latent-contract audit may examine
residual variance, effect/`nbsq` scaling, p1 orientation, pi updates, active and
component-count definitions, and annotation-variance output semantics. Official
multichain parity still requires a seedable official release or a documented
independent-chain mechanism. The final benchmark remains unauthorized.
