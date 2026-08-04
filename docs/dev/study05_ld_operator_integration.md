# Study 05 LD-operator integration

## Outcome

Study 05 is the completed LD-operator validation. It contains one ordinary-list
specification, one exact analysis, one operator-design implementation, one
capsule-only report, and one compact reference capsule. Study 06 is the
annotation-informed-model study and remains in development. Study 07 remains
multitrait validation and remains in development.

## Integrated components

| Component | Preserved evidence |
|---|---|
| `operator_validation` | Two architectures, five paired replicates, BED/full CSR/block CSR/near-full/0.999/0.995 operators, deterministic identities, convergence, recovery, runtime and compatibility evidence |
| `sbayesr_ld_sensitivity` | Scheduler comparison, full-sweep BED, exact CSR, hard-sparse CSR, full-rank block eigen, retained-0.995 block eigen, corrected-score, quadratic/residual, spectral and recovery audits |

All settings, coordinates, seeds, priors, controls, blocks, retained ranks,
tolerances, metrics, and scientific conclusions are unchanged.

## Source consolidation

The former main and nested analyses are integrated into
`studies/05_ld_operator/analysis.R`. Reusable and Study-05-specific operator
construction, block, spectral, and deterministic-audit functions are
consolidated in `operator-design.R`. The nested SBayesR source directory and
separate report were removed. Paused Study 07 now sources only the authoritative
operator-design file.

## Capsule integration

The former main and supplemental capsules were merged into
`results/reference/05_ld_operator/current/`. Main result names are retained;
former supplemental tables use an `sbayesr_` prefix. One integrated manifest,
README, source inventory, and SHA-256 checksum index describe both components.
Obsolete capsule reproduction wrappers were omitted. Numerical and scientific
values were not recalculated; detailed file-level hashes and intentional
metadata changes are recorded in `study05_capsule_integration.csv`.

## Execution and provenance

Future checkpoints use semantic Study 05 identities. Historical local caches
under old numbered paths remain ignored and retired; fit objects were neither
copied nor translated. The validated package remains `sblr` 0.2.0 at SHA
`02e8c74baa906e83c4a08d42a9cc6339b4e81072`. No sampler, benchmark pipeline,
targets pipeline, package build, or package check was run for the integration.

