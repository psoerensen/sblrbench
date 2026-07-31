# Five-replicate overnight run

The Windows launcher performs an offline validation pass before detaching the
sequential Study 02, Study 03, and Study 04 run. It sets four OpenMP threads per
native fit, one BLAS thread, and runs only one four-core method fit at a time.
`tarchetypes` is not required: the active pipelines call `targets` directly and
contain no `tarchetypes` functions.

Validate without sampling:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_five_replicate_overnight.ps1 -ValidateOnly
```

Launch in the background after validation succeeds:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/run_five_replicate_overnight.ps1
```

The foreground fallback is:

```powershell
Rscript scripts/run_five_replicate_overnight.R
```

## Inspect progress

Follow the summary:

```powershell
Get-Content results/local/five_replicate_overnight/overnight_summary.txt -Wait
```

Inspect the machine-readable phase table:

```powershell
Import-Csv results/local/five_replicate_overnight/overnight_status.csv |
  Format-Table
```

Check the detached process:

```powershell
$overnightPid = Get-Content results/local/five_replicate_overnight/process_id.txt
Get-Process -Id $overnightPid
```

Study-specific logs are `study02.log`, `study03.log`, and `study04.log` in
`results/local/five_replicate_overnight`. Rendering, tests, and package-check
logs are `quarto.log`, `tests.log`, and `check.log` in the same directory.

## Stop and resume

Request a clean process stop (Windows does not deliver Unix signals to R):

```powershell
$overnightPid = Get-Content results/local/five_replicate_overnight/process_id.txt
Stop-Process -Id $overnightPid
```

Completed targets remain in the three ignored stores under
`results/local/five_replicate_overnight/_targets_02`, `_targets_03`, and
`_targets_04`. Relaunch the same PowerShell script to resume; `targets` skips
completed branches. Do not delete these stores. A deterministically failing
branch is recorded in its study log and status table before later safe phases
continue.

## Inspect promotion and Git state

Promotion succeeded only when the corresponding directory exists and its
validator passes:

```powershell
Get-ChildItem results/reference/02_prediction/st-bayesc-bayesr-five-replicate-development-v1
Get-ChildItem results/reference/03_parameter_estimation/st-parameter-estimation-five-replicate-development-v1
Get-ChildItem results/reference/04_convergence/st-multichain-convergence-validation-five-replicate-v1
```

The orchestrator contains no staging, commit, push, tag, publish, or deployment
operation. Confirm the resulting unstaged workspace with:

```powershell
git status --short
git diff --stat
git diff --cached --stat
git log -1 --oneline
```

The run never needs interactive approval between studies and never downloads
data or accesses `../sblr`.
