# Current Study 06 annotation-convergence stop

The BayesRC and SBayesRC interface audit and two four-chain maximum-history fits
completed with the explicit mixture grid `c(0, 0.01, 0.1, 1)`. Neither method
passed every prespecified convergence threshold, so the stop rule was triggered.

No five-replicate scientific benchmark was started. This capsule preserves the
annotation design, seeds, raw scalar chains, diagnostics, limiting quantities,
source provenance, and checksums. Reproduce the convergence decision with:

    powershell -ExecutionPolicy Bypass -File scripts/run_current_benchmark_refresh.ps1 -Phase study06 -Resume
