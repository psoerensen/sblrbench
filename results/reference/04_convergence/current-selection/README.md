# Current multichain convergence selection

Four maximum-history fits and 16 identifiable chains select operational MCMC
settings under the frozen diagnostic thresholds. This stage does not assess
prediction or scientific superiority, and passing these datasets does not prove
universal convergence.

The separate `current-validation` capsule evaluates the selected settings in five
replicates. Native fits are excluded. Reproduce both stages with:

    powershell -ExecutionPolicy Bypass -File scripts/run_current_benchmark_refresh.ps1 -Phase study04 -Resume
