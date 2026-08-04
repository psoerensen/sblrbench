# Study 05/06 renumbering audit

## Decision

The committed LD-operator validation becomes **Study 05 — LD-operator
validation**. The annotation-model work becomes **Study 06 —
Annotation-informed models** and remains explicitly in development. Study 07
remains multitrait validation and remains in development.

The previous internal identifiers `06_ld_operator` and
`05_annotation_models` are retired. Repository callers are updated directly;
no compatibility aliases are retained.

## Pre-rename dependency inventory

| Area | Old dependency | Required action |
|---|---|---|
| Study source | `studies/06_ld_operator/` | Move to `studies/05_ld_operator/`; integrate the nested SBayesR analysis, helpers, and report into the top-level workflow |
| Annotation source | `studies/05_annotation_models/` | Move mechanically to `studies/06_annotation_models/`; retain development code and honest status |
| Main LD capsule | `results/reference/06_ld_operator/current/` | Merge into the single `results/reference/05_ld_operator/current/` capsule |
| Supplemental LD capsule | `results/reference/06_ld_operator/sbayesr_ld_robustness/current/` | Merge compact tables into the same Study 05 capsule with collision-safe names |
| Annotation stop capsule | `results/reference/05_annotation_models/current-stop/` | Move mechanically to `results/reference/06_annotation_models/current-stop/`; update path/study metadata and checksums only |
| Shared runner/spec | `R/benchmark-execution.R`, `R/benchmark-spec.R` | Make `05_ld_operator` authoritative and reject retired IDs clearly |
| Root dispatch | `_targets.R` | Replace the LD ID only; do not run or redesign targets |
| Website | `_quarto.yml`, `index.qmd`, `studies/index.qmd`, `reproducibility.qmd` | Renumber studies, remove the Additional validation entry, publish one Study 05 report |
| Annotation launch/refresh | `scripts/run_study05_annotation_models.*`, `scripts/run_current_benchmark_refresh.R` | Rename mechanically to Study 06 paths and labels; do not execute |
| Tests | Study 05 annotation and Study 06 LD tests | Renumber and add integrated-report/capsule assertions |
| Developer documentation | baseline, inventories, dependency map, Study 06 migration and LD robustness notes | Record new numbering and retired historical paths |
| Paused Study 07 | sources/copies Study 06 operator internals | Point mechanically to Study 05 without changing Study 07 science |

## Integration policy

The main and SBayesR evidence streams retain their existing scientific tables.
The final report and analysis present them as components of one Study 05. The
integrated capsule contains one README, one manifest, one checksum index, all
main compact files, and supplemental compact tables prefixed with `sbayesr_`
where necessary. Native fits, matrices, eigenvectors, checkpoints, and logs are
excluded.

Historical ignored local caches under old numbered paths are retired execution
caches. They are not copied or translated.

