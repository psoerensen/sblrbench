# Implementation plan

1. Preserve the existing package metadata and add focused S3 list contracts for simulations, methods, results, and manifests.
2. Implement strict axis alignment, oracle validation, the public-API-only `sblr` adapter/extractor, and four truth-aware long-format metrics.
3. Add fast contract tests, documentation, a small non-MCMC smoke study, and plain-R future-study configurations.
4. Generate roxygen files, run tests/check/smoke from `sblrbench`, then prove the sibling commit and status are unchanged.

The requested structure is used. A small internal helper layer in `R/utils.R` centralizes identifier and matrix validation; no framework or registry is introduced.
