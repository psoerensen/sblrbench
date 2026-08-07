# Studies

Completed validation studies use a small source contract: `spec.R`, a readable
exact `analysis.R`, and a frozen-capsule-only `report.qmd`. Study-specific
scientific logic remains alongside that contract when it does not belong in the
shared framework.

- Studies 01–05 are completed and validated.
- Study 06 annotation-informed models is in development; v2 qualification
  failed and the final benchmark is blocked. The paired isolation and
  package-side audits are complete; official multichain parity is blocked, but
  its single-trajectory descriptive comparison completed. The frozen large
  information-scale feasibility profile completed all six fits: controls pass,
  learned-alpha joint mixing fails (`LARGE-G2`), SNP utility remains useful,
  and route differences persist. Start with the
  [Study 06 navigation page](06_annotation_models/README.md) or
  [report](06_annotation_models/report.qmd).
- Study 07 multitrait validation is in development and paused.

Move mechanics into package-level `R/` only when their shared contract is
demonstrated. This directory is repository research infrastructure and is
excluded from the installed package.
