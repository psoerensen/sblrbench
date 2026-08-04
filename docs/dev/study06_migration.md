# Historical Study 06 LD-operator migration

This document records the framework migration committed while LD-operator
validation still used the internal identifier `06_ld_operator`. The subsequent
renumbering moved that completed scientific work to **Study 05** and moved the
in-development annotation work to **Study 06**. See
[`study05_ld_operator_integration.md`](study05_ld_operator_integration.md) and
[`study05_06_renumbering.md`](study05_06_renumbering.md) for the authoritative
current structure.

The migration established the ordinary-list specification, shared task
dispatch, reusable operator metrics, complete exact analysis, capsule-only
report, and semantic checkpoint policy now located at:

- `studies/05_ld_operator/spec.R`
- `studies/05_ld_operator/analysis.R`
- `studies/05_ld_operator/operator-design.R`
- `studies/05_ld_operator/report.qmd`
- `R/metrics-operator.R`

The scientific contract remains two architectures, five paired benchmark
replicates, six main configurations, the frozen seed/control/tolerance scheme,
and the deterministic 1,500-marker SBayesR window with retained rank 1,490.
Historical source-number and nested-supplement paths are retired and are not
compatibility interfaces.
