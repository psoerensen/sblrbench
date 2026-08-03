# Removed benchmark artifacts

Git history is the archive for earlier benchmark states. The refresh removed
duplicate active capsules and chronology-based reports only after their current
replacement validated and all active dependencies were redirected.

| Former path | Study | Reason removed | Last containing commit before refresh | Current replacement |
|---|---|---|---|---|
| `results/reference/01_finemapping/separated-development-v1` | 01 | Superseded numerical capsule | `ace87b742fbda202dba7aadaf3ac17631e6ce26c` | `results/reference/01_finemapping/current` |
| `results/reference/01_finemapping/separated-development-v1.1` | 01 | Duplicate provenance revision | `32ecbf19e2641f3126196c046e998e05a0aacc95` | `results/reference/01_finemapping/current` |
| `results/reference/02_prediction/st-bayesc-bayesr-one-replicate-development-v1` | 02 | Preliminary one-replicate evidence | `aebd65d9af08ff6ecc5d2389db7fc45faa222d14` | `results/reference/02_prediction/current` |
| `results/reference/02_prediction/st-bayesc-bayesr-five-replicate-development-v1` | 02 | Superseded package/prior run | `39fb5802507fc867b34dc4a32644f72705202686` | `results/reference/02_prediction/current` |
| `results/reference/03_parameter_estimation/st-parameter-estimation-one-replicate-development-v1` | 03 | Preliminary one-replicate evidence | `6296dce0f4c0507227d9d25017b1f977cf35bde9` | `results/reference/03_parameter_estimation/current` |
| `results/reference/03_parameter_estimation/st-parameter-estimation-five-replicate-development-v1` | 03 | Superseded package/prior run | `39fb5802507fc867b34dc4a32644f72705202686` | `results/reference/03_parameter_estimation/current` |
| `results/reference/04_convergence/st-multichain-convergence-development-v1` | 04 | Superseded recommendation selection | `c40eb3893b6ff8946b64bd0650ac6ebd1cc3df31` | `results/reference/04_convergence/current-selection` |
| `results/reference/04_convergence/st-multichain-convergence-validation-five-replicate-v1` | 04 | Superseded fixed-setting validation | `39fb5802507fc867b34dc4a32644f72705202686` | `results/reference/04_convergence/current-validation` |
| `results/reference/06_ld_operator/st-ld-operator-convergence-development-v1` | 06 | Reconstructed-dense historical stage | `a4bd3fdaedbe6e36400c97318299f1965fba72bf` | Git history; current low-rank evidence below |
| `results/reference/06_ld_operator/st-ld-operator-five-replicate-development-v1` | 06 | Reconstructed-dense historical benchmark | `a4bd3fdaedbe6e36400c97318299f1965fba72bf` | `results/reference/06_ld_operator/current` |
| `results/reference/06_ld_operator/st-low-rank-operator-convergence-development-v2` | 06 | Chronology-specific duplicate stage | `285a33ae9fcaba93af2b5bd001c72507da2949e6` | `results/reference/06_ld_operator/current` |
| `results/reference/06_ld_operator/st-low-rank-operator-five-replicate-development-v2` | 06 | Versioned duplicate of compatible current evidence | `285a33ae9fcaba93af2b5bd001c72507da2949e6` | `results/reference/06_ld_operator/current` |
| `studies/01_finemapping/separated-development-pilot.qmd` | 01 | Chronology-based report | `9007f48f2c69fc13cc74c4ae18bc926dc0cc5e25` | `studies/01_finemapping/fine-mapping.qmd` |
| `studies/02_prediction/prediction-development-pilot.qmd` | 02 | Fallback reader for old capsules | `0fffcb87d0c3e85d7685803caa2299cfcd8409fd` | `studies/02_prediction/report.qmd` |
| `studies/03_parameter_estimation/parameter-estimation-development-pilot.qmd` | 03 | Fallback reader for old capsules | `0fffcb87d0c3e85d7685803caa2299cfcd8409fd` | `studies/03_parameter_estimation/report.qmd` |
| `studies/04_convergence/convergence-development-pilot.qmd` | 04 | Version-specific report | `39fb5802507fc867b34dc4a32644f72705202686` | `studies/04_convergence/convergence.qmd` |
| `studies/05_annotation_models/annotation-models-development-pilot.qmd` | 05 | Pre-refresh pilot presentation | `0d5b7d854e655c88aac69cef59279be513f4b37d` | `studies/05_annotation_models/annotation-convergence.qmd` |
| `studies/06_ld_operator/ld-operator-development-pilot.qmd` | 06 | Reconstructed-dense historical report | `a4bd3fdaedbe6e36400c97318299f1965fba72bf` | Git history |
| `studies/06_ld_operator/retained-low-rank-operator-development-v2.qmd` | 06 | Chronology-based report | `285a33ae9fcaba93af2b5bd001c72507da2949e6` | `studies/06_ld_operator/low-rank-operator.qmd` |

Obsolete one-off promotion and overnight-run controllers were removed with the
same cleanup because the current refresh runner and study-specific promotion
contracts replace them. Their last containing commits are recorded in Git:
`scripts/promote_reference_results.R` (`32ecbf1`),
`scripts/promote_prediction_results.R` (`aebd65d`),
`scripts/run_five_replicate_overnight.R` (`39fb580`),
`studies/five_replicate_promotion.R` (`ca78ed2`), and
`studies/03_parameter_estimation/promote_parameter_estimation_results.R`
(`6296dce`). The obsolete execution note
`docs/dev/five_replicate_overnight_run.md` (`39fb580`) was removed with those
controllers; the current refresh record replaces it. Its dedicated obsolete
orchestrator test, `tests/testthat/test-five-replicate-overnight.R`
(`471dbb5`), was also removed after current study and refresh tests covered the
replacement paths.
