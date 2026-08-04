# Study 06 — Annotation-informed models

**Status: In development**

Study 06 contains deterministic annotation construction, provisional
simulation scenarios, BayesRC/SBayesRC method experiments, chain extraction,
convergence diagnostics, and draft recovery metrics. The prespecified
maximum-history convergence gate stopped both annotation-aware methods before a
five-replicate benchmark was started.

The tracked `results/reference/06_annotation_models/current-stop/` directory is
development evidence for that stop decision. It is not an authoritative
completed benchmark capsule and must not be used for method-performance claims.
This is not a completed benchmark.

Scientific questions still unresolved include a supported convergence design,
the final annotation scenario grid, and validation of annotation-effect and
mixture-probability recovery. A future migration should use the shared data,
simulation, method, checkpoint, convergence, capsule, and reporting helpers.
The existing target graph remains only as preserved development machinery; it
is not run by the completed-study CLI.
