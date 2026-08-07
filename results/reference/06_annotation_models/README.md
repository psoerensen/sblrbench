# Study 06 reference evidence

## Current status

**CLOSED — EST-R2**, with official qualifier **EST-R5** and same-posterior
sampler-development endpoint **PMA-R3**.

Post-closure supporting evidence **ALPHA-R1** is registered in
`alpha_recovery_metadata.json`, `alpha_recovery_summary.csv`, and
`alpha_recovery_single_fixture.csv`. These compact files describe an isolated
probit hierarchy reference only; they contain no production SBayesRC draws.

Current compact closure evidence:

- `final_decision.json` — authoritative machine-readable status, source
  identities, frozen artifact hashes, and reproduction scope;
- `final_cross_implementation_comparison.csv` — truth, official D1, learned
  BED, learned block, and fixed-oracle/ablation summaries;
- `final_hierarchy_of_evidence.csv` — reporting interpretation by inferential
  level.

The `current-stop/` capsule is immutable historical v1 failed-qualification
evidence. Its “stop” label is not the current Study 06 status.

## Reproduction

```powershell
Rscript scripts/run_study06_estimability.R
```

Retained frozen chains must already exist under ignored local results. This is
offline posterior analysis only: it invokes no sampler, regenerates no truth,
and writes no scientific fit artifact.
