# Study 04 convergence migration

Study 04 is migrated to the common ordinary-R framework. Its scientific record
remains the frozen selection and validation capsules; neither was regenerated
or promoted.

## Responsibility map

| Previous responsibility | Authoritative replacement |
|---|---|
| `config.R` and diagnostic registry | `studies/04_convergence/spec.R` |
| Study 03 scenario/method copies | matched Study 03 spec plus `benchmark_convergence_coordinates()` |
| chain extraction | `extract_convergence_traces()`; true native traces only |
| windowing and scalar diagnostics | `R/benchmark-convergence.R` |
| burn-in stability and recommendations | `benchmark_burnin_stability()` and `benchmark_convergence_recommendations()` |
| selection and validation targets graphs | convergence task in `run_benchmark()` |
| source-path caches | `sblrbench-semantic-v2` scientific identities |
| launch scripts and worked example | common CLI, exact `analysis.R`, and convergence template |
| promotion wrappers | frozen capsules plus shared checksum validation; no promotion during migration |
| `convergence.qmd` | `report.qmd`, still frozen-capsule-only |

## Scientific preservation

The matched Study 03 grid remains homogeneous BayesC/SBayesC and mixture
BayesR/SBayesR. Selection remains four 3,000-draw, zero-warm-up, four-chain
fits. Validation remains twenty fixed-setting fits over five deterministic
replicates. Simulation seeds, Study 04 fit/chain seeds, priors, controls,
candidate burn-ins (250, 500, 1,000), retained draws (250, 500, 1,000, 2,000),
four required scalar quantities, R-hat/ESS/MCSE definitions, thresholds, and
earliest-stable recommendation rule are unchanged.

The benchmark profile represents the complete historical design: four
selection plus twenty validation coordinates. Workshop represents the four
selection coordinates only and is unsuitable for performance claims.

## Checkpoints and callers

New local fits use semantic identities containing stage, matched coordinate,
ordered samples/markers, simulation state, controls, seeds, diagnostic design,
package SHA, and qgdata provenance. Source paths and timestamps are excluded.
Historical targets stores and ignored caches remain untouched; no legacy loader
was retained. Study 06 now reads the authoritative Study 03 data/scenario spec
directly instead of sourcing Study 04 config. The current-refresh driver uses
the common runner. Root targets dispatch gives a clear message for migrated
Studies 02--04 and remains available to unmigrated studies.

## Validation record

Focused deterministic tests compare the grid and seed map with the historical
formulas, reproduce the frozen recommendation table from frozen diagnostics,
exercise true-trace extraction and convergence fixtures, verify zero-dispatch
validation-only execution, and require table-only plotting/report contracts.
All authoritative capsules are checksum-validated and recursively compared to
HEAD during handoff. No sampler, benchmark pipeline, targets pipeline, capsule
promotion, package build, or package check is run.

Study 01 fine-mapping remains the next migration candidate. Annotation,
operator, and paused MTBLR interfaces remain deferred.
