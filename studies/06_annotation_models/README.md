# Study 06 — Annotation-informed models

## Status

**CLOSED — EST-R2**

Functional annotations produce reproducible SNP-prior rankings and stable
SNP-level inference, while precise recovery of the unrestricted continuous-
alpha decomposition and exact annotation-effect magnitudes is substantially
less reliable.

## Current status

```text
Study 06 status: CLOSED
Primary final decision: EST-R2
Official qualifier: EST-R5
Standard unrestricted continuous-alpha sampler development: closed at PMA-R3
Additional final scientific benchmark: not required
Additional Study 06 scientific fits: not required
Next methodology: separate Bayesian annotation-selection / annotation-PIP study
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

The [final conclusion](../../docs/studies/study06_final_conclusion.md) is the
authoritative scientific endpoint, and the machine-readable
[final decision](../../results/reference/06_annotation_models/final_decision.json)
records its frozen identity. The [main report](report.qmd) provides the readable
study account. The
[evidence synthesis](../../docs/dev/study06_annotation_inference_evidence_synthesis.md)
preserves the chronological technical history; experiment-specific reports and
decision files remain authoritative for their historical gates.

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
  SBayesRC trajectory shows the same benefit descriptively and strongly agrees
  with `sblr` for SNP effects, PIPs, ranks, and validation genetic values.
- At 5,000 people and 37,991 markers, baseline and fixed-true-alpha controls
  pass their aggregate convergence contracts, but both ordinary learned-alpha
  routes still fail joint alpha/allocation convergence.
- The final estimability analysis finds highly reproducible induced prior
  rankings and SNP outputs, directionally stable informative-annotation
  contrasts, and quantitatively uncertain contrast magnitudes: **EST-R2**.

### Descriptive but useful

Learned-alpha AUPRC and prediction summaries are descriptive under
non-convergence. Official alpha, occupancy, architecture, native five-component
results, runtime, and drift are single-trajectory descriptive. Pooled learned
BED/block summaries and per-chain official-versus-`sblr` comparisons are not
formal posterior validation.

### Difficult but closed for Study 06

Unrestricted continuous-alpha joint learning, later-stick decomposition,
global alpha/allocation coupling, quantitative annotation-variance parity, and
exact annotation-effect magnitudes remain difficult. Official independently
seeded multichain inference is unavailable under the v0.2.6 native RNG contract
(**EST-R5**). These limitations do not reopen Study 06: same-posterior sampler
development ended at **PMA-R3**, exact but computationally impractical.

## Authoritative documents

- [Final Study 06 conclusion](../../docs/studies/study06_final_conclusion.md) —
  authoritative current scientific decision.
- [Estimability and annotation contrasts](../../docs/studies/study06_estimability_and_contrasts.md)
  — detailed final evidence.
- [Repository closure](../../docs/dev/study06_repository_closure.md) — current
  provenance, status reconciliation, and preserved-history map.
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
- [GCTB-compatible block contract validation](../../docs/dev/study06_gctb_block_contract_validation_result.md)
  and [SMALL-G1 decision](../../docs/dev/study06_gctb_block_contract_validation_decision.json);
- [large information-scale design](../../docs/dev/study06_large_feasibility_design.md),
  [completed result](../../docs/dev/study06_large_feasibility_result.md), and
  [LARGE-G2 decision](../../docs/dev/study06_large_feasibility_decision.json).

## Recommended reading paths

### Quick scientific overview

1. [Final Study 06 conclusion](../../docs/studies/study06_final_conclusion.md).
2. [Estimability and annotation contrasts](../../docs/studies/study06_estimability_and_contrasts.md).
3. [Official SBayesRC comparison](../../docs/dev/study06_gctb_single_trajectory_result.md).
4. [Large information-scale result](../../docs/dev/study06_large_feasibility_result.md).
5. [Sampler-development endpoint](https://github.com/psoerensen/sblr/blob/master/docs/dev/sbayesrc_sampler_development_endpoint.md).
6. [Historical evidence synthesis](../../docs/dev/study06_annotation_inference_evidence_synthesis.md).

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
- Final [decision JSON](../../results/reference/06_annotation_models/final_decision.json),
  [cross-implementation comparison](../../results/reference/06_annotation_models/final_cross_implementation_comparison.csv),
  and [hierarchy-of-evidence summary](../../results/reference/06_annotation_models/final_hierarchy_of_evidence.csv).
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

## Closure and next methodology

1. Study 06 is closed at **EST-R2**.
2. Standard same-posterior continuous-alpha sampler development is closed at
   **PMA-R3**.
3. No additional Study 06 scientific fit or final benchmark is required.
4. The next methodology is a separate Bayesian annotation-selection /
   annotation-PIP study, not a sampler fix for standard SBayesRC.

Offline reproduction of the final derived analysis uses
`Rscript scripts/run_study06_estimability.R`. Retained frozen chains must
already exist locally. The command invokes no sampler, regenerates no truth,
and writes tables and figures to the ignored
`results/local/06_annotation_models/estimability_and_contrasts/` directory.
