param(
  [Parameter(ValueFromRemainingArguments = $true)]
  [string[]]$ScriptArgs
)

$scriptPath = Join-Path $PSScriptRoot "normal_chamber_analysis.R"
if (-not (Test-Path -LiteralPath $scriptPath)) {
  throw "Cannot find script: $scriptPath"
}

$rscript = $null
$cmd = Get-Command Rscript -ErrorAction SilentlyContinue
if ($cmd) {
  $rscript = $cmd.Source
}

if (-not $rscript) {
  $candidates = Get-ChildItem "$env:ProgramFiles\R" -Recurse -Filter Rscript.exe -ErrorAction SilentlyContinue |
    Sort-Object FullName -Descending
  if ($candidates) {
    $rscript = $candidates[0].FullName
  }
}

if (-not $rscript) {
  throw "Rscript.exe was not found. Install R or add Rscript to PATH first."
}

& $rscript $scriptPath @ScriptArgs
exit $LASTEXITCODE
