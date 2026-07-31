# Study 05 annotation-model development run

## Scope and provenance

Study 05 is a package-specific development benchmark of the installed `sblr`
implementation. It compares ST-BED BayesR, ST-BED BayesRC, ST-CSR SBayesR and
ST-CSR SBayesRC on five paired simulations in each of two annotation scenarios.
It is not a reproduction of another package's implementation, a definitive
method ranking, or universal convergence validation. Block-eigen models are out
of scope.

The initial repository gate was recorded on `master` at
`765346d0d5fd5d63c77f9b537badd415a7a894bc` with an empty
`git status --short`. The installed packages were `sblr` 0.1.2 (source commit
`92ff3f6e7a0b1228f9f04b693d91a36d86934b0f`) and `sblrbench`
0.0.0.9000. The qgdata source commit is
`6cca5819e711d326cfb2614d7e9d9f34942612cd`. The cached BED/BIM/FAM,
phenotype and covariate files matched the pinned sizes and MD5 values in
`config.R`. No network access or package installation is part of this run.

## Verified installed interfaces

| Method | Public entry point | Model string | Annotation argument |
|---|---|---|---|
| ST-BED BayesR | `sblr::stblr_bed()` | `bayesr` | none |
| ST-BED BayesRC | `sblr::stblr_bed()` | `bayesrc` | `annotation`, numeric marker-by-annotation matrix |
| ST-CSR SBayesR | `sblr::stblr_csr()` | `sbayesr` | none |
| ST-CSR SBayesRC | `sblr::stblr_csr_annot()` | `sbayesrc` | `annotations`, numeric marker-by-annotation matrix |

All fits use the explicit component variance grid `c(0, 0.01, 0.1, 1)`.
The first component is the exact null and the remaining components are active.
The annotation matrix has exactly the canonical marker rows, explicit marker
IDs, and columns `Intercept`, `enriched_binary`, `continuous_signal`, and
`null_annotation`, in that order. The intercept is supplied explicitly;
automatic intercept addition and annotation standardization are disabled
because preprocessing is frozen in the benchmark.

For four components, `alpha` is a 4-by-3 matrix: annotation columns by
stick-breaking steps. For step \(k\), `eta = A %*% alpha[, k]` and
`stick = pnorm(eta)`. Component \(k\) receives
`remaining * (1 - stick)` and the remaining mass is multiplied by `stick`;
the fourth component receives the final remainder. Thus step 1 controls
null-versus-non-null mass in the verified component order. The permanent
implementation calls `sblr::sbayesrc_marker_pi()` and tests its result against
this orientation, including bounds and unit row sums.

`sigmaSqAlpha` has one value per stick (length three). The BED public `...`
controls are named `annot_alpha_init`, `annot_sigma_sq_alpha_init`, and
`annot_alpha_update_every`; the CSR annotation entry point uses `alpha_init`,
`sigmaSqAlpha_init`, and `alpha_update_every`. The installed defaults
are shape/rate hyperparameters 2 and 2, initial value 1, a flat intercept
coefficient, and `alpha_update_every = 1`, meaning an annotation update every
MCMC iteration. Study 05 records these controls explicitly.

`nit` is retained post-warm-up draws, `nburn` is additional warm-up,
`nthin = 1` means no thinning, and native multichain execution uses
`nchains = 4`, `ncores = 4`, four distinct `chain_seeds`, and
`keep_chains = TRUE`. Extended convergence retention with the `annotations`
and `probability` groups exposes iteration-by-chain draws for all alpha
coefficients, all `sigmaSqAlpha` values, global mixture summaries, and
marker-averaged prior non-null summaries. This satisfies the chain-level
output gate for the convergence pilot.

Returned objects are kept conceptually separate:

- BED `annotation_prior` contains annotation-implied prior component
  probabilities per marker. CSR does not return that field, so the same
  quantity is reconstructed exactly from its returned posterior-mean `alpha`
  and the frozen annotation matrix with `sbayesrc_marker_pi()`.
- `component_probabilities` are posterior component-allocation probabilities,
  per marker.
- global mixture summaries are posterior global component proportions.
- posterior marker non-null probability is one minus posterior allocation
  probability for the null component.

Marker effects, PIPs/component allocations, variance components, annotation
coefficients and their compact chain summaries are available. A metric is
emitted only from the corresponding returned quantity; unavailable estimands
remain explicit with a machine-readable reason.

## Read-only source audit

The installed interface was checked against the installed-package source commit
using read-only `git show`, plus read-only inspection of the sibling checkout.
No sibling file was changed or loaded as package code. The consulted paths are:

- `R/mtblr-bayesrc.R`
- `R/sbayesrc-helpers.R`
- `R/blr-extended-convergence.R`
- `R/mtblr-bed.R`
- `R/mtblr-csr.R`
- `R/stblr-csr-annot.R`
- `R/stblr-csr-sbayesrc.R`
- `R/sparse_ld_bed_helper.R`
- `tests/testthat/test-stblr-bed-interface.R`
- `tests/testthat/test-stblr-csr-interface.R`
- `tests/testthat/test-blr-extended-parameter-diagnostics.R`
- `tests/testthat/test-blr-extended-probability-diagnostics.R`

The final capsules include this inventory.

## Convergence stop (2026-07-31)

Both maximum-history fits completed and remain cached: ST-BED BayesRC took
19m48s and ST-CSR SBayesRC took 9m07s. Four requested seeds and four distinct
effective package task seeds were retained for each fit. No main-benchmark fit
was started.

No candidate setting passed the unchanged Study 04 thresholds for all required
quantities. This triggered the prespecified stop condition, so neither
convergence evidence nor the five-replicate benchmark was promoted as a frozen
capsule. The local evidence is in:

- `results/local/study05_annotation_models/unsupported_convergence_diagnostics.csv`
- `results/local/study05_annotation_models/unsupported_candidate_settings.csv`
- `results/local/study05_annotation_models/_targets_convergence`

The least-failing BED candidate had burn-in 1,000 and 500 retained draws, with
18 failed quantities; its maximum R-hat was 1.961, minimum bulk ESS 5.65,
minimum tail ESS 10.93, and maximum relative MCSE 0.487. The least-failing CSR
candidate had burn-in 1,000 and 1,000 retained draws, with 23 failed
quantities; its maximum R-hat was 2.960, minimum bulk ESS 4.57, minimum tail
ESS 11.22, and maximum relative MCSE 0.477. Alpha coefficients across all
three sticks, global non-null proportion, annotation-implied enriched and
unenriched non-null summaries, and annotation-variance quantities were among
the recurring limiters. The thresholds were not weakened and the completed
fits were not rerun.

## Frozen simulation contract

One deterministic annotation design is reused by both scenarios and every
replicate. Exactly 10% of markers receive `enriched_binary`; both continuous
columns are independently generated and standardized to sample mean zero and
sample standard deviation one.

The informative non-intercept alpha rows, in stick order 1–3, are:

| annotation | step 1 | step 2 | step 3 |
|---|---:|---:|---:|
| enriched_binary | 0.98 | 0.20 | 0.00 |
| continuous_signal | 0.25 | 0.00 | 0.20 |
| null_annotation | 0.00 | 0.00 | 0.00 |

The first intercept is calibrated once to 50 expected non-null markers. The
remaining intercepts target active component weights 0.60, 0.30 and 0.10. The
uninformative scenario sets every non-intercept coefficient to zero and
reverse-calibrates its intercepts to the informative scenario's marker-average
component probabilities. No replicate-specific calibration or rejection is
allowed. A realized non-null count outside 20–100 is a prespecified failure.

## Execution

Sampler-free validation:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_study05_annotation_models.ps1 --validate-only
```

Detached convergence selection:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_study05_annotation_models.ps1 --phase convergence --resume
```

After the convergence capsule validates, detached benchmark execution:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_study05_annotation_models.ps1 --phase benchmark --resume
```

Run both sequentially in the foreground:

```powershell
Rscript scripts/run_study05_annotation_models.R --phase all --resume
```

Progress and logs:

```powershell
Import-Csv results/local/study05_annotation_models/study05_status.csv | Format-Table
Get-Content results/local/study05_annotation_models/study05_summary.txt -Wait
Get-Content results/local/study05_annotation_models/convergence.log -Wait
Get-Content results/local/study05_annotation_models/benchmark.log -Wait
```

Check the detached process:

```powershell
$pid = Get-Content results/local/study05_annotation_models/process_id.txt
Get-Process -Id $pid
```

Stop it cleanly with `Stop-Process -Id $pid`. Resume with the same phase and
`--resume`; validated target branches remain cached. Capsule promotion occurs
only after its complete grid validates.
