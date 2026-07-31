<#
Foreground fallback:
  Rscript scripts/run_study05_annotation_models.R --phase all --resume

Examples:
  powershell -ExecutionPolicy Bypass -File scripts/run_study05_annotation_models.ps1 --validate-only
  powershell -ExecutionPolicy Bypass -File scripts/run_study05_annotation_models.ps1 --phase convergence --resume
#>
[CmdletBinding()]
param(
  [Alias("validate-only")]
  [switch]$ValidateOnly,
  [ValidateSet("convergence", "benchmark", "all")]
  [string]$Phase = "all",
  [switch]$Resume
)

$ErrorActionPreference = "Stop"
$scriptPath = $MyInvocation.MyCommand.Path
$repoRoot = (Resolve-Path (Join-Path (Split-Path $scriptPath -Parent) "..")).Path
if (-not (Test-Path (Join-Path $repoRoot "sblrbench.Rproj"))) {
  throw "Could not resolve the sblrbench repository root."
}
Set-Location $repoRoot

$rscript = (Get-Command Rscript -ErrorAction SilentlyContinue).Source
if (-not $rscript) {
  throw "Rscript was not found on PATH."
}

$logDir = Join-Path $repoRoot "results/local/study05_annotation_models"
New-Item -ItemType Directory -Path $logDir -Force | Out-Null
$runner = Join-Path $repoRoot "scripts/run_study05_annotation_models.R"
$common = @($runner, "--phase", $Phase)
if ($Resume) { $common += "--resume" }

$env:OMP_NUM_THREADS = "4"
$env:OMP_THREAD_LIMIT = "4"
$env:OPENBLAS_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
$env:VECLIB_MAXIMUM_THREADS = "1"

if ($ValidateOnly) {
  & $rscript @common "--validate-only"
  exit $LASTEXITCODE
}

$pidPath = Join-Path $logDir "process_id.txt"
if (Test-Path $pidPath) {
  $oldPid = [int](Get-Content $pidPath -First 1 -ErrorAction SilentlyContinue)
  if ($oldPid -and (Get-Process -Id $oldPid -ErrorAction SilentlyContinue)) {
    throw "A Study 05 process is already active with PID $oldPid."
  }
}

Write-Host "Running the sampler-free validation gate..."
& $rscript @common "--validate-only"
if ($LASTEXITCODE -ne 0) {
  throw "Validation-only preflight failed; no sampler process was launched."
}

$stdout = Join-Path $logDir "launcher.stdout.log"
$stderr = Join-Path $logDir "launcher.stderr.log"
$process = Start-Process -FilePath $rscript -ArgumentList $common `
  -WorkingDirectory $repoRoot -WindowStyle Hidden -PassThru `
  -RedirectStandardOutput $stdout -RedirectStandardError $stderr
$process.Id | Set-Content -Path $pidPath -Encoding ascii

Start-Sleep -Seconds 2
if (-not (Get-Process -Id $process.Id -ErrorAction SilentlyContinue)) {
  throw "The Study 05 process exited immediately. Inspect $stderr."
}

$statusPath = Join-Path $logDir "study05_status.csv"
$deadline = (Get-Date).AddSeconds(60)
while (-not (Test-Path $statusPath) -and (Get-Date) -lt $deadline) {
  Start-Sleep -Milliseconds 500
}
if (-not (Test-Path $statusPath)) {
  throw "The process is alive, but the status file was not created within 60 seconds."
}

Write-Host "Study 05 started with PID $($process.Id)."
Write-Host "Progress:"
Write-Host "  Import-Csv results/local/study05_annotation_models/study05_status.csv | Format-Table"
Write-Host "  Get-Content results/local/study05_annotation_models/study05_summary.txt -Wait"
Write-Host "  Get-Content results/local/study05_annotation_models/convergence.log -Wait"
Write-Host "  Get-Content results/local/study05_annotation_models/benchmark.log -Wait"
