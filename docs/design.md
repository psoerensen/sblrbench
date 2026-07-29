# Design

Study 01 is a configuration-driven separated-locus vertical path. Replicates branch from shared genotype and sparse-LD targets and cross with four scalar-trait methods. It uses seeded randomized-greedy physical separation, exact-causal simulation verification, canonical-tie marker ranking, public `sblr` credible sets, and exact/reference-LD proxy scoring. Method failures remain structured records. Development MCMC controls are replaceable as one config list and are not scientific defaults.

```text
sblr or external simulation
          ↓
sblrbench_simulation
          ↓
strict alignment and oracle validation
          ↓
sblrbench_method
          ↓
native fit
          ↓
sblrbench_result
          ↓
truth-aware metric rows
          ↓
versioned manifest and summaries
```

Implementation correctness remains in `sblr`; larger benchmarks and robustness studies belong in `sblrbench`. The sibling/source `sblr` is a read-only external dependency and only exported installed-package functions are called. Marker, sample, and trait alignment are strict. Absent quantities remain `NULL`, and method-specific outputs are not falsely equated. Native fits are retained optionally. Future external methods can implement the same fit/extract/predict contract without a registry.
