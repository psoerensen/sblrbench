# Study 06 shared-framework migration

Study 06 remains **Study 06 — LD-operator validation**. Study 05 remains the
annotation-model study and Study 07 remains reserved for future MTBLR
validation. The main retained-low-rank benchmark and supplemental SBayesR
LD-robustness diagnostic retain separate frozen capsules.

## Responsibility map

| Historical responsibility | Authoritative location after migration |
|---|---|
| Design values spread across `config.R` and `v2/config.R` | `studies/06_ld_operator/spec.R` |
| Block and retained-eigen scientific policy | `studies/06_ld_operator/operator-design.R` |
| Reusable matrix/action error metrics | `R/metrics-operator.R` |
| Coordinates, seeds, profiles, validation | `R/benchmark-spec.R` |
| Task dispatch and validation-only output | `R/benchmark-execution.R` |
| Design tables and named plots | `R/benchmark-reporting.R` |
| Exact main workflow | `studies/06_ld_operator/analysis.R` |
| Main capsule-only report | `studies/06_ld_operator/report.qmd` |
| Supplemental deterministic audits | `studies/06_ld_operator/sbayesr_ld_robustness/` |
| Reusable example | `inst/templates/operator-analysis.R` |

The per-study targets graph, launchers, promotion code, duplicated simulation,
method and diagnostic helpers, contract smoke test, and v2 orchestration were
removed. `operators.R` and `operator_validation.R` remain because paused Study
07 still calls their exact operator functions. They contain real
operator logic, not active Study 06 orchestration; consolidating that dependency
is deferred to Study 07 to avoid changing its scientific contract here.

## Preserved scientific contract

- Two architectures, five paired benchmark replicates, and six configurations
  produce 60 benchmark coordinates (12 for structural workshop validation).
- The 70/30 split uses seed 3101 and preserves the pinned qgdata sample and
  chromosome-1 marker order.
- The simulation mixture remains 0.60/0.30/0.10 with multipliers 0.01/0.1/1.
  The distinct BayesR prior remains 0.99 inactive with active multipliers
  0.01/0.1/1.
- Blocks remain contiguous at 1,000 markers. `block_low_rank_v1` retains
  near-full, 0.999, and 0.995 positive spectral mass and rebuilds residuals
  every 100 updates.
- Seed arithmetic, four-chain policy, convergence thresholds, equivalence
  tolerances, projected SSE identity, recovery metrics, and conclusions are
  unchanged.
- The supplemental 1,500-marker window, seed 17002, four chain seeds,
  exact/sparse/block definitions, 0.995 policy, and rank 1490 are unchanged.

## Checkpoints, capsules, and validation

New work uses the shared `sblrbench-semantic-v2` identity with operator type,
ordered inputs, blocks, eigen policy, retained mass/rank, controls, seeds,
package SHA and data provenance. Source and documentation paths are excluded;
historical caches remain retired and were not translated.

The main numerical fits retain `sblr` SHA
`bd8e2c8148a0d9540dc20716455706beeb16fa86`; compatibility was validated at
`02e8c74baa906e83c4a08d42a9cc6339b4e81072`. The supplemental capsule records
the latter SHA. No capsule was regenerated or promoted.

Validation comprises capsule checksums and recursive SHA-256 comparison to
HEAD, deterministic fixtures, validation-only CLI runs, package tests, parsing,
and capsule-only rendering. No sampler, benchmark/targets pipeline, package
build, or package check is run. Study 05 and Study 07 remain deferred.
