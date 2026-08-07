# Studies

Completed validation studies use a small source contract: `spec.R`, a readable
exact `analysis.R`, and a frozen-capsule-only `report.qmd`. Study-specific
scientific logic remains alongside that contract when it does not belong in the
shared framework.

- Studies 01–06 are completed and validated within their stated scopes.
- Study 06 annotation-informed BayesRC/SBayesRC is closed at **EST-R2**:
  probability rankings and SNP inference are stable, while exact annotation-
  contrast magnitudes remain uncertain. The official qualifier is **EST-R5**
  because formal independently seeded multichain replication is unavailable;
  its retained trajectory is descriptive evidence. Same-posterior sampler
  development ended at **PMA-R3**. Start with the
  [final conclusion](../docs/studies/study06_final_conclusion.md),
  [Study 06 navigation page](06_annotation_models/README.md), or
  [report](06_annotation_models/report.qmd).
- Study 07 multitrait validation is in development and paused.

Move mechanics into package-level `R/` only when their shared contract is
demonstrated. This directory is repository research infrastructure and is
excluded from the installed package.
