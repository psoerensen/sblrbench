# Study 06 — Annotation-informed models

**Status: In development**

Study 06 contains deterministic synthetic annotation construction, two
provisional simulation scenarios, matched BayesR/BayesRC and SBayesR/SBayesRC
comparisons, chain extraction, convergence diagnostics, and draft recovery
metrics. The prespecified maximum-history convergence gate stopped both
annotation-aware methods before the five-replicate benchmark was started.

The tracked `results/reference/06_annotation_models/current-stop/` directory is
development evidence for that stop decision. It is not an authoritative
completed benchmark capsule and must not be used for method-performance claims.
This is not a completed benchmark.

The repository audit found that the intended package interfaces exist at the
pinned `sblr` SHA, but annotation-parameter convergence and consistent
draw-based prior-probability extraction remain completion gates. The old local
target cache is not accepted as reusable scientific evidence.

The concrete audit and completion plan are documented in
[`docs/dev/study06_annotation_audit.md`](../../docs/dev/study06_annotation_audit.md)
and
[`docs/dev/study06_annotation_completion_plan.md`](../../docs/dev/study06_annotation_completion_plan.md).
Future work should migrate the study directly to the shared data, simulation,
method, semantic-checkpoint, convergence, capsule, and reporting framework.
No final validated benchmark is currently claimed.
