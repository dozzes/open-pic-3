# Verifies a behavior-preserving change is bit-exact against a stored baseline
# case by re-running the case with the freshly built exe and diffing grids
# against the pre-existing diag/ output already sitting in that case folder.
#
# Usage (run manually -- this script does not build anything, per project
# convention that only the user triggers cmake/msbuild):
#   cmake --build build --config Release -j 8
#   powershell -ExecutionPolicy Bypass -File .\tools\Compare-BaselineAfterFix.ps1

param(
    [string]$CaseDir = "sim\tasks\01_BASELINE\MA2_D44_UNIF_COLD",
    [string]$Stamp = (Get-Date -Format "yyyyMMdd_HHmmss")
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$OpenPic = Join-Path $RepoRoot "build\bin\Release\open-pic.exe"
$Compare = Join-Path $RepoRoot "tools\Compare-GridOutputs.ps1"
$CaseRoot = Join-Path $RepoRoot $CaseDir

if (-not (Test-Path -LiteralPath $OpenPic)) {
    throw "Non-MPI build not found: $OpenPic -- build it first (cmake --build build --config Release -j 8)"
}
if (-not (Test-Path -LiteralPath $CaseRoot)) {
    throw "Case folder not found: $CaseRoot"
}

$diagDir = Join-Path $CaseRoot "diag"
if (-not (Test-Path -LiteralPath $diagDir)) {
    throw "No existing diag/ output in $CaseRoot to compare against -- nothing to diff"
}

# Preserve the pre-change output under a stamped name so repeated runs of this
# script never collide with each other or with the unrelated scratch_before/
# scratch_after folders already present under 01_BASELINE.
$backupDir = Join-Path $CaseRoot "diag_before_$Stamp"
Write-Host "Backing up existing output: diag -> diag_before_$Stamp"
Move-Item -LiteralPath $diagDir -Destination $backupDir

Write-Host "Copying freshly built exe into case folder and re-running"
Copy-Item -LiteralPath $OpenPic -Destination (Join-Path $CaseRoot "open-pic.exe") -Force
Push-Location $CaseRoot
try {
    & ".\open-pic.exe" "main.lua"
    $exitCode = $LASTEXITCODE
}
finally {
    Pop-Location
}

if ($exitCode -ne 0) {
    throw "open-pic.exe exited with code $exitCode -- see console output above"
}
if (-not (Test-Path -LiteralPath $diagDir)) {
    throw "Re-run did not produce a new diag/ folder: $diagDir"
}

Write-Host ""
Write-Host "Comparing $backupDir vs $diagDir"
& $Compare $backupDir $diagDir
$compareOk = $LASTEXITCODE -eq 0

Write-Host ""
if ($compareOk) {
    Write-Host "PASS: grid outputs match bit-for-bit"
    exit 0
}
else {
    Write-Host "FAIL: grid outputs differ -- see Compare-GridOutputs.ps1 output above"
    exit 1
}
