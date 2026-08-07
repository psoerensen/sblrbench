# AGENTS.md

## Repository role

`sblrbench` is the scientific validation and benchmarking repository for `sblr`.

Statistical models, samplers, inference algorithms, and production implementation belong in the sibling `sblr` repository. `sblrbench` owns simulation designs, benchmark specifications, validation analyses, metrics, provenance, compact reference evidence, and reports.

Do not modify `sblr` from this repository unless the task explicitly requests separate package-development work.

## Read before changing a study

Before modifying an existing study, read:

1. the repository `README.md`;
2. `framework.qmd`;
3. `reproducibility.qmd`;
4. the study's `README.md` or status page;
5. its `spec.R`;
6. any current final decision, closure, or authoritative study document referenced by that README.

Historical documents under `docs/dev/` preserve development history and may describe superseded states. Do not treat an older development decision as the current roadmap when a later closure or final-decision document exists.

## Study design

Scientific choices for a benchmark belong in the study `spec.R`.

For a new benchmark, define the scientific design, comparison models, simulation scenarios, primary metrics, provenance pins, and decision criteria before generating primary benchmark results.

Do not change primary hypotheses, scenarios, metrics, or qualification gates after inspecting primary results without recording the change explicitly as a new design version or addendum.

Reuse shared mechanics under `R/` where their contracts apply. Keep genuinely study-specific scientific logic within the study.

## Implementation versus benchmarking

Use `sblr` to answer:

> Is the statistical implementation correct?

Use `sblrbench` to answer:

> How does the validated method behave scientifically?

A benchmark failure must not be silently repaired by changing the statistical model or sampler from within `sblrbench`.

If benchmark evidence suggests a package defect, document the evidence and stop or route the implementation work to `sblr`.

## Provenance and reproducibility

Every benchmark must record the exact `sblr` source revision and relevant execution provenance.

Use semantic checkpoints and existing repository provenance mechanisms.

Working runs belong under `results/local/` and remain non-authoritative.

Tracked evidence under `results/reference/` must be compact, reviewed, reproducible evidence rather than raw fit objects or arbitrary working output.

Do not promote or replace a reference capsule unless the task explicitly calls for that reviewed action.

## Completed studies

Preserve completed and closed studies.

Do not silently rewrite historical decisions because later work changes their interpretation. Add a clearly identified addendum, closure update, or subsequent study instead.

Current study status is determined by the latest authoritative study README/status page, final decision, and closure documents; historical `docs/dev/` material remains evidence of the development path.

## Repository safety

Before substantial work, record the branch, HEAD, and working-tree status.

Do not overwrite unrelated local changes.

Do not commit, push, tag, publish, or deploy unless explicitly requested.

Keep generated local results, build products, caches, and temporary files out of tracked source.

Run the relevant tests and validation checks before reporting completion.
