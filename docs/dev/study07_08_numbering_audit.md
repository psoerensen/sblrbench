# Study 07/08 numbering audit

## Current numbering

- Study 06: annotation-informed BayesRC/SBayesRC closure and sampler audit.
- Study 07: one-replicate SBayesRC joint-versus-EM demonstration.
- Study 08: paused multitrait validation scaffold.

The multitrait source tree and its focused test were moved mechanically from
`studies/07_mt_validation/` to `studies/08_mt_validation/` and from
`test-study07-mt-validation.R` to `test-study08-mt-validation.R`. Active study
IDs, paths, headings, launch guards, navigation, and test references now use
Study 08. Seeds, model controls, scientific assertions, expected behavior, and
the paused status were not changed. No multitrait analysis was run and no
multitrait result capsule existed to rename.

**STUDY08-MOVE PASS — existing MTBLR Study 07 mechanically renumbered to Study
08.**

## Intentionally retained historical references

Development inventories and reorganization reports under `docs/dev/` retain
references to the former Study 07 multitrait scaffold when they describe the
repository state at that historical commit. In particular, the
`final_cleanup_*`, `reorganization_*`, `study01_migration.md`,
`study05_06_renumbering.md`, `study05_ld_operator_integration.md`, and
`study06_annotation_completion_plan.md` records are historical audit evidence,
not active navigation or executable identifiers. Frozen artifact names would
also be retained if present; this audit found no Study 07 multitrait reference
capsule.

Current navigation, execution guards, the implementation map, and tests identify
the paused multitrait work as Study 08. A repository search after the move found
no ambiguous active Study-07/MTBLR reference.
