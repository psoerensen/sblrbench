# Study 06 — Annotation-informed models

**Status: In development**

Study 06 now has a final software specification and shared-framework
implementation for deterministic annotations, matched BayesR/BayesRC and
SBayesR/SBayesRC comparisons, comparable draw-wise annotation-prior extraction,
convergence qualification, semantic checkpoints, recovery metrics, and exact
analysis/report contracts. Scientific qualification has not been run.

The tracked `results/reference/06_annotation_models/current-stop/` directory is
development evidence for that stop decision. It is not an authoritative
completed benchmark capsule and must not be used for method-performance claims.
This is not a completed benchmark.

The pinned `sblr` interfaces expose the selected annotation coefficient traces
needed for a consistent draw-wise marker-prior estimand. If a future fit lacks
those true retained traces, extraction reports the quantity as unavailable; it
does not substitute a final state or posterior mean.

The next required scientific action is a four-entry qualification run. Only a
passing, identity-matched qualification decision can authorize the separate
40-fit, 160-chain final benchmark. The old local target cache is retired.

The concrete audit and completion plan are documented in
[`docs/dev/study06_annotation_audit.md`](../../docs/dev/study06_annotation_audit.md)
and
[`docs/dev/study06_annotation_completion_plan.md`](../../docs/dev/study06_annotation_completion_plan.md).
The study now uses the shared data, simulation, method, semantic-checkpoint,
convergence, extraction, metric, and reporting framework. The implementation
record is in
[`docs/dev/study06_annotation_implementation.md`](../../docs/dev/study06_annotation_implementation.md).
No final validated benchmark is currently claimed.
