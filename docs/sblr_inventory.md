# `sblr` inventory

## Immutable reference state

- Source repository: `../sblr`
- Commit inspected: `fd76b18828bc77756948aa3138de07ae4dc75513`
- Initial `git status --short`: empty (clean)
- Installed package inspected: `sblr` 0.1.0

The source tree was inspected read-only. The installed package is used for all execution. It is older than the inspected source: it exports `mtsim()`, `make_summary_stats()`, `make_sparse_ld()`, `stblr_csr()`, `stblr_csr_annot()`, `stblr_bed()`, `check_stblr_consistency()`, and credible-set helpers, but does not export the source tree's `stblr_block_eigen()`, `mtblr_csr()`, `mtblr_block_eigen()`, or `mtblr_bed()`. Those adapters can be specified but execution is deferred until a compatible installed `sblr` is available.

## Public contracts

The current source defines seven canonical fit interfaces: `stblr_csr(stats, Glist = NULL, ld_prefix = NULL, method, ...)`, `stblr_csr_annot(stats, Glist = NULL, annotations, annotation_model, method = NULL, ...)`, `stblr_block_eigen(stats, Glist, block_start, method, ...)`, `stblr_bed(y, Glist, method, ...)`, `mtblr_csr(stats, Glist = NULL, ld_prefix = NULL, ..., method)`, `mtblr_block_eigen(stats, Glist, block_start, ..., method)`, and `mtblr_bed(y, Glist, ..., method)`. `make_summary_stats()` and `make_sparse_ld()` are the current preparation names. `make_credible_sets()` and related exported helpers own credible-set semantics. `check_stblr_consistency()` is an ST structural checker, not a convergence diagnostic or universal MT checker.

`mtsim()` returns `y`, standardized or raw `W` according to its controls, marker-by-trait `B`, sample-by-trait `G` and `E`, `causal`, `rsids`, `ids`, target/observed `h2`, `rg`, and `re`, plus `Sigma_e`, shared/specific indices, causal marker IDs, and observed effect-correlation summaries. Rows of `W`, `G`, `E`, and `y` follow `ids`; columns of `W` and rows of `B` follow `rsids`; trait columns follow the simulation trait names. Its causal list has `shared`, a trait-named `specific` list, and `all`.

Canonical fits preserve marker and trait order in `fit$data` metadata. `fit$bm` is the marker-by-trait posterior mean effective effect. `fit$dm` is the marker-by-trait posterior probability of binary activity/non-null status, including marginalized MT pattern/component states; it is not a component label. MT `fit$cov_g_mean` is the full posterior-mean trait covariance matrix. `fit$convergence` is the convergence summary; `fit$diagnostics` holds native execution/numerical diagnostics; `fit$convergence_traces`, optional compact `fit$chains`, and analytical `fit$memory_estimate` have separate meanings.

## Responsibility boundary and useful references

`sblrbench` calls the exported functions above directly. It does not reproduce simulation kernels, summary statistics, sparse LD, block eigen preparation, BED processing, samplers, fit formatting, convergence algorithms, structural diagnostics, or credible-set algorithms. The scripts in `examples/workflows/`, including `workflow_helpers.R`, are examples rather than package APIs and are not sourced. Small patterns adapted here are limited to explicit named-list interfaces and preservation of canonical marker/trait names.

Relevant reference tests include the CSR/BED/block-eigen interface tests, unified public-contract and convergence tests, annotation-interface tests, MT model/operator tests, backend consistency tests, and credible-set tests. They confirm that implementation correctness and model-specific semantics remain in `sblr`.

Public-API limitations relevant here are the installed/source export mismatch noted above, absence of a universal MT structural checker, and model-specific posterior fields that cannot safely be normalized (`component_probabilities`, prior probabilities, `pi_*`, pattern/annotation probabilities, and annotation coefficients). These remain only in optional `native_fit`. External methods and full prediction/fine-mapping studies are deferred.

## Completion verification

- Final source commit: `fd76b18828bc77756948aa3138de07ae4dc75513`
- Final `git status --short`: empty (clean)
- The final commit and status are identical to the recorded initial values.
