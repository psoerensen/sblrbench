# Study 07 — SBayesRC joint versus EM inference

## Status

**COMPLETE — STUDY07-R5 (review required)**

Study 07 reuses exactly one frozen Study 06 simulated dataset to ask whether
the completed MCEM inference lines reach stable annotation-informed solutions
on the benchmark that exposed difficult joint learned-alpha/allocation mixing.
Inference chains and dispersed EM starts are diagnostics on the same dataset;
they are not additional simulated-data replicates.

The five prespecified methods were SBayesR, joint SBayesRC, SBayesRC-EM, joint
SBayesRC-S, and SBayesRC-S-EM. The study distinguishes joint posterior
inference from MCEM/MAP inference, and full posterior `annotation_pip` from the
responsibility-conditioned MCEM-Laplace `annotation_pip_eb`.

All four dispersed starts for both EM lines reached the qualified 50-outer
iteration cap without satisfying convergence tolerances, and both methods had
materially start-dependent annotation/prior endpoints. Genomic inference was
often sensible, but the prespecified stability demonstration failed. See the
[concise report](report.qmd) for the evidence and limits.

The scientific design is frozen in [spec.R](spec.R). Working fits, tables, and
figures are written below `results/local/07_joint_em_sbayesrc/`. They remain
ignored working evidence; no reference capsule was promoted before review.

Run the gates in order with:

```powershell
Rscript scripts/run_study07_joint_em.R --gate baseline
Rscript scripts/run_study07_joint_em.R --gate joint
Rscript scripts/run_study07_joint_em.R --gate em
Rscript scripts/run_study07_joint_em.R --gate selection
Rscript scripts/run_study07_joint_em.R --gate summarize
```

No command simulates data or modifies `sblr`.
