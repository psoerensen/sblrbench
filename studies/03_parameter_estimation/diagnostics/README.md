# CSR SBayesR / GCTB prior diagnostic

This developer diagnostic fixes one existing Study 03 coordinate and was designed
to isolate BayesR mixture-proportion prior strength, `adjE`, and initial sparsity.
It is not a benchmark refresh and never reads or writes a benchmark target store
or reference capsule.

Run from the repository root:

```powershell
Rscript studies/03_parameter_estimation/diagnostics/sbayesr-gctb-comparison.R
```

Generated inputs, checkpoints, tables, figures, logs, provenance, and session
information belong under the ignored directory:

```text
results/local/sbayesr_gctb_diagnostic/
```

The script requires the isolated `sblr` 0.2.0 installation at source SHA
`02e8c74baa906e83c4a08d42a9cc6339b4e81072`. It reuses the Study 03 data,
simulation, summary-statistic, sparse-LD, method-seed, and chain-seed contracts.

## Isolation gate

Before fitting, the script resolves the baseline BayesR variance priors and
checks whether the public CSR SBayesR API accepts explicit `B`, `E`,
`ssb_prior`, and `sse_prior` arguments. All four must be fixed across variants;
otherwise changing `alpha` also changes the variance prior.

At the pinned package SHA, `sblr::stblr_csr_bayesr()` does not expose these
arguments and rejects unknown `...` arguments. The wrapper
`sblr::stblr_csr()` cannot forward them through its internal implementation
either. The script therefore records a controlled `blocked_before_fitting`
result and exits successfully without calling a sampler. This behavior is
required by the diagnostic protocol; bypassing it with an unexported sampler
would no longer test the supported public interface.

If a future public interface supports all four explicit prior quantities, the
gate must be reviewed together with the checkpoint identity before any fit is
allowed to run.
