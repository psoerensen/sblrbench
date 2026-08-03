# Current single-trait prediction benchmark

Fresh five-replicate evidence for two architectures and four methods: 40
validated four-chain fits using the current Study 04 recommendations and pinned
`sblr` source recorded in `benchmark_manifest.json`.

This capsule contains compact tabular outputs and provenance only; native fits,
sparse LD, and targets stores are excluded. Reproduce with:

    powershell -ExecutionPolicy Bypass -File scripts/run_current_benchmark_refresh.ps1 -Phase study02 -Resume
