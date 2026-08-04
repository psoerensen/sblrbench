# Study 01 fine-mapping migration

Study 01 now uses `spec.R`, the common runner and CLI, semantic checkpoints,
shared qgdata/alignment and method helpers, shared extraction, and
`R/metrics-finemapping.R`. The report remains frozen-capsule-only.

The migration preserves chromosome 1; 5,000 samples; 37,991 ordered markers;
the separated ten-causal-marker design; target heritability 0.2; ten benchmark
replicates; causal, simulation, method, and chain seeds; the four BED/CSR
BayesC/BayesR methods; 500 retained and 250 burn-in iterations in one chain;
and the native `sblr::make_credible_sets()` policy. Separated causal selection
and locus boundaries remain in `locus-design.R`.

| Retired source | Authoritative replacement |
|---|---|
| `config.R` | `spec.R` |
| `targets.R` | `run_benchmark()` and the common CLI |
| `setup_example_data.R` | `R/benchmark-data.R` |
| `pilot.R` | shared mechanics, `locus-design.R`, and `metrics-finemapping.R` |
| `promotion.R` | shared capsule mechanics; no promotion occurred |
| `fine-mapping.qmd` | `report.qmd` |

Studies 03, 05, 06, and paused Study 07 were updated mechanically to use
shared data helpers; their scientific designs were not migrated. Historical
targets caches are retired internal artifacts. New local fits use semantic-v2
identities based on scientific inputs, not source paths. Frozen capsules remain
authoritative and were not regenerated.
