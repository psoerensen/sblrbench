# Design

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
