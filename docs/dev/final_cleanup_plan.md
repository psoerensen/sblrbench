# Final repository cleanup plan

This pass starts from clean commit `71f627d` and preserves the frozen science
of completed Studies 01--05. Study 06 and Study 07 remain development work.

## Decisions

1. Remove the `00_contract_smoke` pseudo-study and standalone smoke scripts.
   Preserve useful deterministic assertions in `tests/testthat/`.
2. Remove report and five-replicate compatibility loaders after updating every
   active caller to package/shared functions.
3. Limit root targets dispatch to the genuinely active Study 06 and Study 07
   development graphs. Completed studies use `scripts/run_benchmark.R`.
4. Remove migration-era refresh, per-study launch, promotion, targets, and
   smoke infrastructure that is neither scientific source nor authoritative
   evidence.
5. Preserve Study 06 annotation design/simulation/method/diagnostic prototypes
   and Study 07 MT state/simulation/alignment/method/diagnostic prototypes.
6. Replace development launch instructions with concise README status pages.
7. Retain all useful shared `R/` responsibility files and documented public
   APIs; no new architecture is introduced.
8. Rewrite the root README and website text around validation and practical
   workflows, then validate all capsules, reports, links, CLI profiles, tests,
   and parses without fitting.

Uncertain scientific material is retained. Git history remains the archive for
deleted internal orchestration and superseded run documentation.
