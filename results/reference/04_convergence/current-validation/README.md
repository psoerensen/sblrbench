# Current five-replicate convergence validation

Twenty method fits and 80 identifiable chains evaluate the settings selected in
the `current-selection` capsule. All four method recommendations were supported
in all five tested replicates. This is development evidence, not universal
convergence validation.

Native fits are excluded. Reproduce with:

    powershell -ExecutionPolicy Bypass -File scripts/run_current_benchmark_refresh.ps1 -Phase study04 -Resume
