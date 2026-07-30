# Single-trait prediction one-replicate development benchmark

## Benchmark status

This is a complete full-size one-replicate development benchmark.

## Scientific question

BayesC versus BayesR and BED versus CSR for held-out single-trait prediction.

## What is complete

Two architectures, one replicate per architecture, four methods, eight successful fits, complete metrics and paired comparisons, and passing leakage and oracle checks.

## What is not claimed

This capsule provides no replicate-to-replicate uncertainty, convergence evidence, final method ranking, general runtime superiority, or comparison beyond single-trait models.

## Reproduction

The public simulated qgdata files are pinned to commit 6cca5819e711d326cfb2614d7e9d9f34942612cd and checksum validated.
The clean analysis source commit is ea978c67d207297dc5bbbc3bb46ad6d9c5c2b5bb .
Run `run_prediction_benchmark.R` from a compatible clone. Valid targets are reused; development settings are short and single-chain. Small numerical differences may occur across platforms and numerical libraries.

No explicit data licence was identified for qgdata; clarify reuse terms before redistribution.
