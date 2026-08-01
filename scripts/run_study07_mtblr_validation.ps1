# Foreground fallback:
# Rscript scripts/run_study07_mtblr_validation.R --phase all --resume
[CmdletBinding()]
param(
  [switch]$ValidateOnly,
  [switch]$Resume,
  [ValidateSet('audit','contract','runtime','convergence','benchmark','stress','aggregate','verify','all')]
  [string]$Phase = 'all'
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$rscript = (Get-Command Rscript -ErrorAction Stop).Source
$localDir = Join-Path $root 'results/local/study07_mtblr_validation'
New-Item -ItemType Directory -Force -Path $localDir | Out-Null
$pidPath = Join-Path $localDir 'run.pid'
if (Test-Path $pidPath) {
  $oldPid = [int](Get-Content $pidPath -First 1)
  if (Get-Process -Id $oldPid -ErrorAction SilentlyContinue) {
    throw "Another Study 07 process is active: $oldPid"
  }
}
$env:OMP_NUM_THREADS = '4'
$env:OMP_THREAD_LIMIT = '4'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:VECLIB_MAXIMUM_THREADS = '1'
$arguments = @('scripts/run_study07_mtblr_validation.R', '--phase', $Phase)
if ($ValidateOnly) { $arguments += '--validate-only' }
if ($Resume) { $arguments += '--resume' }
if ($ValidateOnly) {
  & $rscript @arguments
  exit $LASTEXITCODE
}
$stdout = Join-Path $localDir 'launcher_stdout.log'
$stderr = Join-Path $localDir 'launcher_stderr.log'
$process = Start-Process -FilePath $rscript -ArgumentList $arguments `
  -WorkingDirectory $root -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr -PassThru -WindowStyle Hidden
$process.Id | Set-Content $pidPath
Start-Sleep -Seconds 2
if (-not (Get-Process -Id $process.Id -ErrorAction SilentlyContinue)) {
  throw "Study 07 process did not remain alive; inspect $stderr"
}
Write-Host "Study 07 started with PID $($process.Id)."
Write-Host "Status: Import-Csv results/local/study07_mtblr_validation/status.csv | Format-Table"
Write-Host "Summary: Get-Content results/local/study07_mtblr_validation/summary.txt -Wait"
Write-Host "Stop: Stop-Process -Id (Get-Content results/local/study07_mtblr_validation/run.pid)"
