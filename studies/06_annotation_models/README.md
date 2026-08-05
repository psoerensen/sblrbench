# Study 06 — Annotation-informed models

**Status: In development — versioned qualification work**

- v1 sparse qualification: **failed and preserved**;
- v2 identifiable qualification: **failed with recorded convergence/mixing,
  scientific-recovery, and route-agreement blockers**;
- final benchmark: **not authorized**.

The tracked
`results/reference/06_annotation_models/current-stop/` directory is immutable
failed v1 development evidence. It is not a completed benchmark capsule and
does not support method-performance claims. The original approximately
37,991-marker, 50-active-marker design is also recorded as `v1_sparse_stress`:
useful for sparse late-stick, prediction, PIP, active-count, variance, and
numerical-integrity stress, but not the primary qualification of precise
late-stick annotation recovery.

The current specification is `v2_identifiable_qualification`. It deterministically
selects 1,500 chromosome-1 QC markers in 15 separated 100-marker blocks,
targets 180 non-null markers with expected active counts 90/54/36, and uses
matched informative and uninformative annotation scenarios. Deterministic
truth checks prevent an accidentally unidentified qualification replicate.

The four scientific routes are public `sblr::stblr_bed()` BayesR/BayesRC and
public `sblr::stblr_block_eigen()` SBayesR/SBayesRC with canonical
`representation = "low_rank"` and `eigen_prop = 0.995`. The current proper
annotation-intercept prior default is used. Approximate sparse CSR, pair
allocation, and unimplemented collapsed block allocation are not v2 gates.

V2 outputs and checkpoints must use the
`results/local/06_annotation_models/v2_identifiable_qualification/` namespace;
the runner rejects v1/current-stop collisions. Validation-only is the default
analysis mode and performs no sampler call. The registered four-chain
qualification has now run. All four histories completed, but every entry
failed the convergence gate; informative late-stick directions and both
BED/block-eigen heritability comparisons also failed. The final benchmark
remains blocked.

See the [v2 design](../../docs/dev/study06_v2_design.md) for the full scientific
contract, audits, decision rule, and diagnostic isolation profiles. The
[v2 qualification result](../../docs/dev/study06_v2_qualification_result.md)
records the exact failed decision. The historical
[v1 qualification result](../../docs/dev/study06_annotation_qualification_result.md)
remains the authoritative record of why v1 stopped.
