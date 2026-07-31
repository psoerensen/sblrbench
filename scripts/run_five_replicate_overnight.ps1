<#
Foreground fallback:
  Rscript scripts/run_five_replicate_overnight.R
#>
[CmdletBinding()]
param([switch]$ValidateOnly, [switch]$SkipCheck)

$ErrorActionPreference = 'Stop'
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Resolve-Path (Join-Path $scriptDir '..')
$rscript = Get-Command Rscript -ErrorAction SilentlyContinue
if (-not $rscript) { throw 'Rscript was not found on PATH.' }
$logDir = Join-Path $repoRoot 'results/local/five_replicate_overnight'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$runner = Join-Path $repoRoot 'scripts/run_five_replicate_overnight.R'
$arguments = @($runner)
if ($ValidateOnly) { $arguments += '--validate-only' }
if ($SkipCheck) { $arguments += '--skip-check' }

if ($ValidateOnly) {
  Push-Location $repoRoot
  try { & $rscript.Source @arguments; exit $LASTEXITCODE } finally { Pop-Location }
}

$stdout = Join-Path $logDir 'launcher_stdout.log'
$stderr = Join-Path $logDir 'launcher_stderr.log'
$process = Start-Process -FilePath $rscript.Source -ArgumentList $arguments `
  -WorkingDirectory $repoRoot -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr -WindowStyle Hidden -PassThru
$process.Id | Set-Content -Path (Join-Path $logDir 'process_id.txt')

$deadline = (Get-Date).AddSeconds(60)
$statusPath = Join-Path $logDir 'overnight_status.csv'
do {
  Start-Sleep -Seconds 1
  if ($process.HasExited) { throw "Overnight runner exited early with code $($process.ExitCode)." }
  if (Test-Path $statusPath) {
    $preflight = Import-Csv $statusPath | Where-Object { $_.phase -eq 'preflight' } | Select-Object -Last 1
    if ($preflight -and $preflight.status -eq 'failed') { throw "Overnight preflight failed: $($preflight.message)" }
    if ($preflight -and $preflight.status -eq 'success') { break }
  }
} while ((Get-Date) -lt $deadline)
if (-not $preflight -or $preflight.status -ne 'success') { throw 'Timed out waiting for successful preflight status.' }

Write-Host "Five-replicate overnight run started with PID $($process.Id)."
Write-Host "Progress: Import-Csv '$statusPath' | Format-Table"
Write-Host "Summary:  Get-Content '$(Join-Path $logDir 'overnight_summary.txt')' -Wait"
Write-Host "Process:  Get-Process -Id $($process.Id)"
