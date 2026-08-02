param(
  [ValidateSet('preflight','audit','study04-selection','study04-validation','study03',
    'study02','study01','study05','study06-compatibility','cleanup','verify','all')]
  [string]$Phase = 'audit',
  [switch]$Resume
)

$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
Set-Location $root
$library = Join-Path $root 'results/local/current_benchmark_refresh/rlib'
if (-not (Test-Path $library)) { throw "The isolated refresh library is missing: $library" }
$env:R_LIBS = $library
$env:OMP_NUM_THREADS = '4'
$env:OMP_THREAD_LIMIT = '4'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:VECLIB_MAXIMUM_THREADS = '1'
$arguments = @('scripts/run_current_benchmark_refresh.R', '--phase', $Phase)
if ($Resume) { $arguments += '--resume' }
& Rscript @arguments
if ($LASTEXITCODE -ne 0) { throw "Current benchmark refresh phase failed: $Phase" }
