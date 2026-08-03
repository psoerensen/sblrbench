param(
  [Parameter(Mandatory = $true)][string]$Study,
  [ValidateSet('workshop', 'benchmark')][string]$Profile = 'benchmark',
  [Parameter(Mandatory = $true)][string]$OutputDir,
  [bool]$Resume = $true,
  [bool]$ValidateOnly = $false
)

$scriptPath = Join-Path $PSScriptRoot 'run_benchmark.R'
& Rscript $scriptPath `
  --study $Study `
  --profile $Profile `
  --output-dir $OutputDir `
  --resume $Resume.ToString().ToLowerInvariant() `
  --validate-only $ValidateOnly.ToString().ToLowerInvariant()
exit $LASTEXITCODE
