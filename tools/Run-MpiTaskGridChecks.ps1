param(
    [string]$Stamp = (Get-Date -Format "yyyyMMdd_HHmmss"),
    [int]$MpiRanks = 2,
    [int]$TimeSteps = 2,
    [double]$BreakTimes = 0.8,
    [int]$CloudPartsOnStep = 3,
    [int]$BackgrPartsOnStep = 1,
    [double]$H_R0 = 4.0,
    [switch]$KeepPassedArtifacts,
    [string[]]$TaskMain
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

. (Join-Path $PSScriptRoot "lib\OpenPicTools.ps1")

$RepoRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "..")).Path
$OpenPic = Join-Path $RepoRoot "build-mpi\bin\Release\open-pic.exe"
$MpiExec = "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"
$Compare = Join-Path $RepoRoot "tools\Compare-GridOutputs.ps1"

if (-not (Test-Path -LiteralPath $OpenPic)) {
    throw "MPI build not found: $OpenPic"
}
if (-not (Test-Path -LiteralPath $MpiExec)) {
    throw "MS-MPI mpiexec not found: $MpiExec"
}

function Get-DefaultTasks {
    # Dedicated MPI regression cases live in sim\tasks\MPI_Tests -- one short
    # case per distinct engine code path (cold FDTD baseline, thermal
    # electrons, PSTD solver, background flow, dipole field). Production
    # cases are NOT enumerated; check one explicitly via -TaskMain.
    # Adding a case = dropping a new <NAME>\main.lua folder into MPI_Tests.
    $testsRoot = Join-Path $RepoRoot "sim\tasks\MPI_Tests"
    if (-not (Test-Path -LiteralPath $testsRoot)) {
        throw "Dedicated MPI test case folder not found: $testsRoot"
    }

    $mains = Get-ChildItem -LiteralPath $testsRoot -Recurse -Filter main.lua |
        Where-Object { $_.FullName -notmatch "MPI_CHECK_" } |  # skip harness run artifacts
        ForEach-Object { $_.FullName }

    if (-not $mains) {
        throw "No main.lua cases found under $testsRoot"
    }

    $mains
}

function Convert-ToShortCase {
    param(
        [string]$SourceMain,
        [string]$TargetMain
    )

    $taskName = Split-Path -Leaf (Split-Path -Parent $SourceMain)
    $content = Get-Content -LiteralPath $SourceMain -Raw

    $content = Set-LuaAssignment -Content $content -Name "break_times" -Value "$BreakTimes" -SourceName $taskName
    $content = Set-LuaAssignment -Content $content -Name "h_R0" -Value "$H_R0" -SourceName $taskName
    $content = Set-LuaAssignment -Content $content -Name "cloud_parts_on_step" -Value "$CloudPartsOnStep" -SourceName $taskName
    $content = Set-LuaAssignment -Content $content -Name "backgr_parts_on_step" -Value "$BackgrPartsOnStep" -SourceName $taskName
    $content = Set-LuaAssignment -Content $content -Name "cloud_jitter_enabled" -Value "false" -SourceName $taskName -Optional
    $content = Set-LuaAssignment -Content $content -Name "backgr_jitter_enabled" -Value "false" -SourceName $taskName -Optional
    $content = Set-LuaAssignment -Content $content -Name "pstd_quick_test_steps" -Value "$TimeSteps" -SourceName $taskName -Optional

    $content += @"

-- MPI regression harness limits. Original task file is unchanged.
pic_parameters.time_steps = math.min(pic_parameters.time_steps, $TimeSteps)
pic_parameters.save_time_steps = 1
"@

    Set-Content -LiteralPath $TargetMain -Value $content -Encoding UTF8
}

function Invoke-CheckedCommand {
    param(
        [string]$FilePath,
        [string[]]$Arguments,
        [string]$WorkingDirectory,
        [string]$LogPath
    )

    $exitCode = Invoke-OpenPicProcess -FilePath $FilePath -Arguments $Arguments -WorkingDirectory $WorkingDirectory -LogPath $LogPath

    return [pscustomobject]@{
        Ok = ($exitCode -eq 0)
        ExitCode = $exitCode
        Log = $LogPath
    }
}

$tasks = if ($TaskMain -and $TaskMain.Count -gt 0) { $TaskMain } else { Get-DefaultTasks }
$results = New-Object System.Collections.Generic.List[object]

foreach ($main in $tasks) {
    $sourceMain = (Resolve-Path -LiteralPath $main).Path
    $sourceDir = Split-Path -Parent $sourceMain
    $parentDir = Split-Path -Parent $sourceDir
    $taskName = Split-Path -Leaf $sourceDir
    $caseRoot = Join-Path $parentDir ("MPI_CHECK_{0}_{1}" -f $taskName, $Stamp)
    $serialDir = Join-Path $parentDir ("MPI_CHECK_{0}_{1}_serial" -f $taskName, $Stamp)
    $mpiDir = Join-Path $parentDir ("MPI_CHECK_{0}_{1}_mpi{2}" -f $taskName, $Stamp, $MpiRanks)

    New-Item -ItemType Directory -Force -Path $caseRoot | Out-Null
    New-Item -ItemType Directory -Force -Path $serialDir | Out-Null
    New-Item -ItemType Directory -Force -Path $mpiDir | Out-Null
    Convert-ToShortCase -SourceMain $sourceMain -TargetMain (Join-Path $serialDir "main.lua")
    Convert-ToShortCase -SourceMain $sourceMain -TargetMain (Join-Path $mpiDir "main.lua")

    Write-Host ""
    Write-Host "== $taskName =="

    $serial = Invoke-CheckedCommand `
        -FilePath $OpenPic `
        -Arguments @("main.lua") `
        -WorkingDirectory $serialDir `
        -LogPath (Join-Path $caseRoot "serial.log")

    if (-not $serial.Ok) {
        Write-Host "serial FAIL exit=$($serial.ExitCode)"
        $results.Add([pscustomobject]@{ Task = $taskName; Status = "serial_fail"; Path = $caseRoot })
        continue
    }

    $mpi = Invoke-CheckedCommand `
        -FilePath $MpiExec `
        -Arguments @("-n", [string]$MpiRanks, $OpenPic, "-mpi", "main.lua") `
        -WorkingDirectory $mpiDir `
        -LogPath (Join-Path $caseRoot "mpi.log")

    if (-not $mpi.Ok) {
        Write-Host "mpi FAIL exit=$($mpi.ExitCode)"
        $results.Add([pscustomobject]@{ Task = $taskName; Status = "mpi_fail"; Path = $caseRoot })
        continue
    }

    # In-process: Compare-GridOutputs always exits explicitly, so $LASTEXITCODE
    # is fresh, and a nested powershell.exe per task is not needed. A thrown
    # error must degrade to compare_fail, not abort the whole task loop.
    $compareOutput = @()
    $compareOk = $false
    try {
        $compareOutput = & $Compare $serialDir $mpiDir -Quiet
        $compareOk = $LASTEXITCODE -eq 0
    } catch {
        $compareOutput = @("compare error: $($_.Exception.Message)")
    }
    if ($compareOk) {
        Write-Host "PASS"
        $shownPath = $caseRoot
        if (-not $KeepPassedArtifacts) {
            Remove-Item -LiteralPath $caseRoot, $serialDir, $mpiDir -Recurse -Force
            $shownPath = "(removed)"
        }
        $results.Add([pscustomobject]@{ Task = $taskName; Status = "pass"; Path = $shownPath })
    }
    else {
        Write-Host "compare FAIL"
        $compareOutput | ForEach-Object { Write-Host $_ }
        $results.Add([pscustomobject]@{ Task = $taskName; Status = "compare_fail"; Path = $caseRoot })
    }
}

Write-Host ""
Write-Host "Summary:"
$results | Format-Table -AutoSize

$failedResults = @($results | Where-Object { $_.Status -ne "pass" })
if ($failedResults.Count -gt 0) {
    exit 1
}

exit 0
