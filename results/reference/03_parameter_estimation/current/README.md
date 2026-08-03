# Current single-trait parameter-estimation benchmark

Fresh five-replicate evidence for two architectures and four methods: 40
validated four-chain fits using the current Study 04 recommendations and pinned
`sblr` source recorded in `benchmark_manifest.json`.

The capsule retains compact truth, posterior, recovery, paired-comparison,
runtime, seed, source, and checksum tables. Native fits and targets stores are
excluded. Reproduce with:

    powershell -ExecutionPolicy Bypass -File scripts/run_current_benchmark_refresh.ps1 -Phase study03 -Resume
