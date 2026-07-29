# Separated-locus fine-mapping development pilot

This frozen reference capsule records a 10-replicate structural benchmark of four single-trait sblr methods on 37,991 canonical markers with 10 separated causal markers per replicate.

The 500-iteration, 250-burn-in, one-chain settings are for development validation only. The example markers have limited LD; these results do not support scientific method rankings.

## Workflow reproduction

The code and configuration can be used to rerun the same analytical workflow with compatible genotype input.

## Exact numerical reproduction

Exact numerical reproduction requires the same original genotype data, software versions, seeds and computing environment. The genotype data are not included in the public reference snapshot.

Run `source("minimal_example.R")` for the small synthetic contract demonstration. Run `source("run_benchmark.R")` only in a suitable clone with the required non-redistributed genotype/Glist input.

The capsule includes compact metrics, summaries, statuses, the frozen manifest/configuration, reproduction scripts, source inventory, session information and checksums. It excludes genotype data, sparse LD, fits, posterior samples and `_targets/`.

Source snapshot: `b10167876d741c1977b8d91779d66043d7567d15`; fit provenance: `f178b332fb3ac6cee95914d519a84301a305ff24`; sblr `0.1.2` (`92ff3f6e7a0b1228f9f04b693d91a36d86934b0f`); qgg `1.1.6`.

Verify a file with `tools::md5sum()` and compare it with `checksums.csv`. The checksum table intentionally does not checksum itself.
