# Study 03 diagnostic checkpoint retirement

Historical diagnostic checkpoints identified reuse through hashes of Study 03
implementation files. Those ignored local caches are execution artifacts, not
authoritative scientific records. The frozen reference capsules, manifests,
checksums, committed specifications, reports, provenance, and Git history remain
the scientific record.

| caller | old dependency | purpose | replacement | checkpoint impact | action |
|---|---|---|---|---|---|
| `studies/03_parameter_estimation/diagnostics/sbayesr-gctb-comparison.R` | `config.R`, `simulation.R`, `methods.R`, `pilot.R` and their hashes | exact Study 03 coordinate and future diagnostic-fit identity | `spec.R`, shared data/simulation/extraction/checkpoint helpers | old coordinate and fit caches are not reusable | semantic identity; explicit legacy rejection |
| `studies/06_ld_operator/sbayesr_ld_robustness/scripts/scheduler-diagnostic.R` | Study 03 config/simulation/pilot sources and prior coordinate | scheduler comparison coordinate and fits | Study 03 `spec.R`, shared simulation/statistic helpers, semantic coordinate and fit identities | old scheduler caches and the old source-hashed Study 03 coordinate are not reusable | semantic identity; explicit legacy rejection |
| `studies/06_ld_operator/sbayesr_ld_robustness/scripts/exact-sparse-diagnostic.R` | Study 03 config/simulation/pilot sources and prior coordinate | reduced exact/sparse operator inputs, LD, and fits | Study 03 `spec.R`, shared checkpoint helpers, explicit scientific payloads | old subset, simulation, LD, and fit caches are not reusable | semantic identity; explicit legacy rejection |

## Decision

The active diagnostic schema is `sblrbench-semantic-v2`. Identities contain a
diagnostic ID and only relevant scientific/computational inputs: coordinates,
ordered sample and marker identities or hashes, simulation/phenotype/truth
hashes where available, LD settings and file hashes, methods, priors and
controls, simulation/fit/chain seeds, the installed `sblr` SHA, qgdata
provenance, and the schema version. Script, helper, report, documentation,
working-directory, and timestamp fields are rejected.

No legacy loader or object translation is provided. Encountering a checkpoint
without the semantic schema produces a clear retirement error. A future
explicit diagnostic run may create new semantic checkpoints; no such fit was
run during this cleanup. Existing ignored checkpoint directories were left
untouched.

Retired local cache roots observed during the cleanup were
`results/local/sbayesr_gctb_diagnostic/`,
`results/local/bed_vs_csr_bayesr_scheduler/`,
`results/local/bed_vs_csr_bayesr_exact_ld/`, and
`results/local/06_ld_operator/sbayesr_ld_robustness/`.
