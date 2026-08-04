# Final reorganization summary

## Architecture

Completed Studies 01--05 use ordinary-list specifications, shared functions
under `R/`, readable exact analyses, the common `run_benchmark()` entry point,
semantic checkpoints, task-specific metrics, and frozen-capsule-only reports.
Study 01 retains locus design locally; Study 05 retains operator design locally.

The final completed studies are fine-mapping, prediction, parameter estimation,
convergence, and integrated LD-operator validation. Study 06 annotation-informed
models and Study 07 multitrait validation remain explicitly in development.

## Cleanup decisions

- Removed the `00_contract_smoke` pseudo-study and standalone smoke scripts;
  useful sampler-free assertions remain in focused `testthat` tests.
- Removed report and five-replicate compatibility loaders after callers moved
  to authoritative shared implementations.
- Removed superseded refresh, overnight, worked-example, Study 06, and Study 07
  launcher scripts.
- Retained root targets support only for the two development study graphs.
- Preserved Study 06 annotation and Study 07 multitrait scientific prototypes,
  while adding honest development README/status documentation.
- Retained five concise templates and the single completed-study CLI.
- Removed superseded operational migration documents while preserving baseline,
  migration, checkpoint, renumbering, capsule, and scientific decision records.

## CLI and website

`scripts/run_benchmark.R` supports benchmark and workshop profiles for Studies
01--05. Study 06 and Study 07 identifiers fail with an in-development message;
retired identifiers fail with their authoritative replacement. The website
separates completed validation studies from development work and contains no
Additional validation category or separate SBayesR page.

## Scientific record

Authoritative completed-study capsules remain under
`results/reference/01_finemapping/` through `05_ld_operator/`. The Study 06
`current-stop` directory is labelled development evidence, not a completed
benchmark capsule. Study 07 has no authoritative capsule. No scientific result,
coordinate, seed, prior, control, metric, or conclusion was changed or rerun.

## Validation

Final validation covers package loading, the focused test suite, parsing of all
modified R sources, benchmark and workshop validation-only CLI calls, retired
and development ID rejection, all capsule checksums, capsule-only report
dependencies, five direct completed-report renders, a full website render, and
a rendered local-link audit. No sampler, benchmark pipeline, targets pipeline,
package build, or package check is part of this cleanup.

The focused suite passed 580 tests. The validation-only coordinate counts were
4/40, 8/40, 8/40, 4/24, and 12/60 for workshop/benchmark profiles in Studies
01--05, respectively. Seven checksum inventories verified 164 reference files.
The full website produced 12 HTML pages; the local-link audit found zero missing
links and zero absolute Windows paths. Generated `_site/` output and the ten
validation-only result directories were removed after validation, while
`_targets/`, `_freeze/`, and historical ignored scientific caches were kept.
