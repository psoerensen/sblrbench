# Study 06 — Annotation-informed models

**Status: In development — qualification failed**

Study 06 now has a final software specification and shared-framework
implementation for deterministic annotations, matched BayesR/BayesRC and
SBayesR/SBayesRC comparisons, comparable draw-wise annotation-prior extraction,
convergence qualification, semantic checkpoints, recovery metrics, and exact
analysis/report contracts. The prespecified four-entry qualification ran on
2026-08-04 and failed, so the final benchmark was not launched.

The tracked `results/reference/06_annotation_models/current-stop/` directory is
development evidence for that stop decision. It is not an authoritative
completed benchmark capsule and must not be used for method-performance claims.
This is not a completed benchmark.

The pinned `sblr` interfaces expose the selected annotation coefficient traces
needed for a consistent draw-wise marker-prior estimand. If a future fit lacks
those true retained traces, extraction reports the quantity as unavailable; it
does not substitute a final state or posterior mean.

Both BED BayesRC histories completed but failed the convergence thresholds;
both CSR SBayesRC entries failed with an invalid projected residual scale. The
next required action is a focused package-side mixing and residual-scale
investigation, followed by the unchanged four-entry qualification. Only a
passing, identity-matched decision can authorize the separate 40-fit,
160-chain final benchmark. The old local target cache is retired.

The concrete audit and completion plan are documented in
[`docs/dev/study06_annotation_audit.md`](../../docs/dev/study06_annotation_audit.md)
and
[`docs/dev/study06_annotation_completion_plan.md`](../../docs/dev/study06_annotation_completion_plan.md).
The study now uses the shared data, simulation, method, semantic-checkpoint,
convergence, extraction, metric, and reporting framework. The implementation
record is in
[`docs/dev/study06_annotation_implementation.md`](../../docs/dev/study06_annotation_implementation.md).
The qualification evidence and exact blockers are recorded in
[`docs/dev/study06_annotation_qualification_result.md`](../../docs/dev/study06_annotation_qualification_result.md).
No final validated benchmark is currently claimed.
