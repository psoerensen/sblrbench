<#
Study 06 detached launcher.

Foreground fallback:
  Rscript scripts/run_study06_ld_operator.R --phase all --resume

Validation-only:
  powershell -ExecutionPolicy Bypass -File scripts/run_study06_ld_operator.ps1 --validate-only
#>
param(
  [Parameter(ValueFromRemainingArguments = $true)]
  [string[]] $RunnerArguments
)

$ErrorActionPreference = "Stop"
$scriptPath = $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent (Split-Path -Parent $scriptPath)
Set-Location -LiteralPath $repoRoot

$localDir = Join-Path $repoRoot "results/local/study06_ld_operator"
New-Item -ItemType Directory -Force -Path $localDir | Out-Null
$pidPath = Join-Path $localDir "run.pid"
$statusPath = Join-Path $localDir "status.csv"
$stdoutPath = Join-Path $localDir "launcher.stdout.log"
$stderrPath = Join-Path $localDir "launcher.stderr.log"

if (Test-Path -LiteralPath $pidPath) {
  $priorPid = Get-Content -LiteralPath $pidPath -TotalCount 1
  if ($priorPid -match '^\d+$') {
    $prior = Get-Process -Id ([int]$priorPid) -ErrorAction SilentlyContinue
    if ($null -ne $prior) {
      throw "A Study 06 process is already active with PID $priorPid."
    }
  }
}

$rscript = Get-Command Rscript -ErrorAction SilentlyContinue
if ($null -eq $rscript) {
  throw "Rscript was not found on PATH."
}

$env:OMP_NUM_THREADS = "4"
$env:OMP_THREAD_LIMIT = "4"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:VECLIB_MAXIMUM_THREADS = "1"

$runner = Join-Path $repoRoot "scripts/run_study06_ld_operator.R"
$arguments = @($runner) + @($RunnerArguments)
$validateOnly = $RunnerArguments -contains "--validate-only"

if ($validateOnly) {
  $priorErrorAction = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  & $rscript.Source @arguments 1> $stdoutPath 2> $stderrPath
  $exitCode = $LASTEXITCODE
  $ErrorActionPreference = $priorErrorAction
  if ($exitCode -ne 0) {
    throw "Study 06 validation-only run failed with exit code $exitCode. Inspect $stderrPath and results/local/study06_ld_operator/operator_validation.log."
  }
  Write-Host "Study 06 validation-only run passed."
  exit 0
}

$process = Start-Process -FilePath $rscript.Source `
  -ArgumentList $arguments `
  -WorkingDirectory $repoRoot `
  -RedirectStandardOutput $stdoutPath `
  -RedirectStandardError $stderrPath `
  -WindowStyle Hidden `
  -PassThru

Set-Content -LiteralPath $pidPath -Value $process.Id -Encoding ascii

$deadline = (Get-Date).AddSeconds(30)
do {
  Start-Sleep -Milliseconds 500
  $alive = Get-Process -Id $process.Id -ErrorAction SilentlyContinue
} until (($null -ne $alive -and (Test-Path -LiteralPath $statusPath)) -or
  (Get-Date) -ge $deadline)

if ($null -eq (Get-Process -Id $process.Id -ErrorAction SilentlyContinue)) {
  throw "Study 06 process exited during launch. Inspect $stderrPath."
}
if (-not (Test-Path -LiteralPath $statusPath)) {
  throw "Study 06 process is alive, but status.csv was not created within 30 seconds."
}

Write-Host "Study 06 started with PID $($process.Id)."
Write-Host "Status: Import-Csv results/local/study06_ld_operator/status.csv | Format-Table"
Write-Host "Summary: Get-Content results/local/study06_ld_operator/summary.txt -Wait"
Write-Host "Operator log: Get-Content results/local/study06_ld_operator/operator_validation.log -Wait"
Write-Host "Convergence log: Get-Content results/local/study06_ld_operator/convergence.log -Wait"
Write-Host "Benchmark log: Get-Content results/local/study06_ld_operator/benchmark.log -Wait"
Write-Host "Stop cleanly: Stop-Process -Id (Get-Content results/local/study06_ld_operator/run.pid)"
