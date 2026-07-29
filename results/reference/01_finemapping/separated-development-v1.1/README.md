# Separated-locus fine-mapping development pilot v1.1

This capsule contains the same numerical benchmark results as `separated-development-v1`. Version 1.1 improves public-data provenance, examples, and reproduction guidance; it does not rerun or alter the benchmark.

The benchmark remains a structural development run with limited LD, 500 iterations, 250 burn-in iterations, one chain, and no method-ranking claims.

## Worked user example

`worked_finemapping_example.R` downloads the public simulated PLINK files, creates and QC-filters a Glist, simulates one phenotype, fits a real ST-BED BayesC model, and evaluates causal recovery. Its short settings demonstrate the workflow and are not convergence recommendations.

## Developer contract smoke test

`contract_smoke_test.R` quickly checks simulation, oracle, result-object, and metric contracts without fitting a sampler.

## Complete benchmark reproduction

`run_benchmark.R` runs the configuration-driven 10-replicate targets workflow. Downloaded qgdata files are cached and checksum-validated; `_targets/` should normally be retained so only missing or outdated work runs.

## Public simulated data

The five inputs are publicly accessible from `psoerensen/qgdata` at commit `6cca5819e711d326cfb2614d7e9d9f34942612cd`. `example_data_manifest.csv` records pinned URLs, sizes, MD5 checksums, and roles. The cached benchmark files match that revision exactly.

No explicit data licence was found in the qgdata repository root at the pinned revision. Public accessibility is verified, but reuse terms should be clarified with the repository owner.

Exact numerical reproduction uses the pinned files, checksums, package versions, seeds, platform, compiler, and numerical libraries. Small numerical differences may occur across platforms.

The capsule excludes genotype data, sparse LD, fit objects, posterior samples, and `_targets/`. Verify files against `checksums.csv`; that table intentionally excludes itself.
