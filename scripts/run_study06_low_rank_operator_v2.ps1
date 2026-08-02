param(
  [Parameter()][ValidateSet('audit','deterministic','operator-pilot',
    'convergence','benchmark','aggregate','report','verify','all')]
  [string]$Phase = 'all',
  [switch]$Resume,
  [Alias('validate-only')][switch]$ValidateOnly,
  [switch]$Foreground
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location $root
$local = Join-Path $root 'results/local/study06_low_rank_operator_v2'
$logs = Join-Path $local 'logs'
New-Item -ItemType Directory -Path $logs -Force | Out-Null
$rscript = (Get-Command Rscript -ErrorAction Stop).Source
$arguments = @('scripts/run_study06_low_rank_operator_v2.R', '--phase', $Phase)
if ($Resume) { $arguments += '--resume' }
if ($ValidateOnly) { $arguments += '--validate-only' }
$env:OMP_NUM_THREADS = '4'
$env:OMP_THREAD_LIMIT = '4'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:VECLIB_MAXIMUM_THREADS = '1'
$studyLib = Join-Path $local 'rlib'
if (-not (Test-Path (Join-Path $studyLib 'sblr/DESCRIPTION'))) {
  throw "Pinned Study 06 v2 installed sblr library is unavailable."
}
$env:R_LIBS_USER = "$studyLib;$env:LOCALAPPDATA\R\win-library\4.4"

if ($Foreground -or $ValidateOnly -or $Phase -in @('audit','deterministic')) {
  & $rscript @arguments
  exit $LASTEXITCODE
}

$pidPath = Join-Path $local 'run.pid'
if (Test-Path $pidPath) {
  $oldPid = [int](Get-Content $pidPath -First 1)
  if (Get-Process -Id $oldPid -ErrorAction SilentlyContinue) {
    throw "Study 06 v2 process $oldPid is already active."
  }
}
$stdout = Join-Path $logs "$Phase.stdout.log"
$stderr = Join-Path $logs "$Phase.stderr.log"
$process = Start-Process -FilePath $rscript -ArgumentList $arguments `
  -WorkingDirectory $root -RedirectStandardOutput $stdout `
  -RedirectStandardError $stderr -WindowStyle Hidden -PassThru
$process.Id | Set-Content $pidPath
Start-Sleep -Seconds 2
if (-not (Get-Process -Id $process.Id -ErrorAction SilentlyContinue)) {
  throw "Study 06 v2 process exited before launch verification."
}
Write-Host "Study 06 v2 $Phase phase started with PID $($process.Id)."
Write-Host "Status: Import-Csv results/local/study06_low_rank_operator_v2/status.csv | Format-Table"
Write-Host "Log: Get-Content $stdout -Wait"
