param(
    [Parameter(Mandatory = $true)]
    [string]$RunA,

    [Parameter(Mandatory = $true)]
    [string]$RunB,

    [double]$Atol = 1e-10,
    [double]$Rtol = 1e-8,

    [switch]$Quiet
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

. (Join-Path $PSScriptRoot "lib\OpenPicTools.ps1")

function Get-GridFiles {
    param([string]$Root)

    $rootPath = (Resolve-Path -LiteralPath $Root).Path
    $files = @{}

    Get-ChildItem -LiteralPath $rootPath -Recurse -File -Filter "*grd*.dat" | ForEach-Object {
        $relative = $_.FullName.Substring($rootPath.Length).TrimStart("\", "/").Replace("\", "/")
        $files[$relative] = $_.FullName
    }

    return $files
}

function Get-RelativeError {
    param(
        [double]$Diff,
        [double]$A,
        [double]$B
    )

    $scale = [Math]::Max([Math]::Max([Math]::Abs($A), [Math]::Abs($B)), 1.0)
    return $Diff / $scale
}

function Compare-GridFile {
    param(
        [string]$Name,
        [string]$PathA,
        [string]$PathB,
        [double]$Atol,
        [double]$Rtol
    )

    $tableA = Read-GridTable -Path $PathA
    $tableB = Read-GridTable -Path $PathB

    if (($tableA.Header -join "`t") -ne ($tableB.Header -join "`t")) {
        return [pscustomobject]@{ Ok = $false; Message = "${Name}: header mismatch" }
    }

    if ($tableA.Rows.Count -ne $tableB.Rows.Count) {
        return [pscustomobject]@{
            Ok = $false
            Message = "${Name}: row count mismatch $($tableA.Rows.Count) != $($tableB.Rows.Count)"
        }
    }

    $maxAbs = 0.0
    $maxRel = 0.0
    $maxLocation = ""

    for ($rowIndex = 0; $rowIndex -lt $tableA.Rows.Count; ++$rowIndex) {
        $rowA = $tableA.Rows[$rowIndex]
        $rowB = $tableB.Rows[$rowIndex]

        for ($colIndex = 0; $colIndex -lt $rowA.Count; ++$colIndex) {
            $a = $rowA[$colIndex]
            $b = $rowB[$colIndex]

            if ([double]::IsNaN($a) -or [double]::IsNaN($b) -or
                [double]::IsInfinity($a) -or [double]::IsInfinity($b)) {
                if ($a -ne $b) {
                    $lineNo = $rowIndex + 1
                    $colName = $tableA.Header[$colIndex]
                    return [pscustomobject]@{
                        Ok = $false
                        Message = "${Name}: non-finite mismatch at row ${lineNo}, column ${colName}"
                    }
                }
                continue
            }

            $diff = [Math]::Abs($a - $b)
            $rel = Get-RelativeError -Diff $diff -A $a -B $b

            if ($diff -gt $maxAbs -or $rel -gt $maxRel) {
                $maxAbs = [Math]::Max($maxAbs, $diff)
                $maxRel = [Math]::Max($maxRel, $rel)
                $lineNo = $rowIndex + 1
                $maxLocation = "row ${lineNo}, column $($tableA.Header[$colIndex])"
            }

            if ($diff -gt $Atol -and $rel -gt $Rtol) {
                $lineNo = $rowIndex + 1
                $colName = $tableA.Header[$colIndex]
                return [pscustomobject]@{
                    Ok = $false
                    Message = ("{0}: tolerance exceeded at row {1}, column {2}: a={3:R}, b={4:R}, abs={5:E3}, rel={6:E3}" -f
                        $Name, $lineNo, $colName, $a, $b, $diff, $rel)
                }
            }
        }
    }

    $detail = "max_abs={0:E3}, max_rel={1:E3}" -f $maxAbs, $maxRel
    if ($maxLocation.Length -gt 0) {
        $detail = "${detail} at ${maxLocation}"
    }

    return [pscustomobject]@{
        Ok = $true
        Message = "${Name}: OK (${detail})"
    }
}

$filesA = Get-GridFiles -Root $RunA
$filesB = Get-GridFiles -Root $RunB

$namesA = [System.Collections.Generic.HashSet[string]]::new([string[]]$filesA.Keys)
$namesB = [System.Collections.Generic.HashSet[string]]::new([string[]]$filesB.Keys)

$missingInB = @($filesA.Keys | Where-Object { -not $filesB.ContainsKey($_) } | Sort-Object)
$missingInA = @($filesB.Keys | Where-Object { -not $filesA.ContainsKey($_) } | Sort-Object)
$common = @($filesA.Keys | Where-Object { $filesB.ContainsKey($_) } | Sort-Object)

$ok = $true

foreach ($name in $missingInB) {
    $ok = $false
    Write-Host "Missing in B: $name"
}
foreach ($name in $missingInA) {
    $ok = $false
    Write-Host "Missing in A: $name"
}

$compared = 0
foreach ($name in $common) {
    $result = Compare-GridFile -Name $name -PathA $filesA[$name] -PathB $filesB[$name] -Atol $Atol -Rtol $Rtol
    ++$compared
    if (-not $result.Ok) {
        $ok = $false
    }
    if (-not $Quiet -or -not $result.Ok) {
        Write-Host $result.Message
    }
}

Write-Host ""
Write-Host "Compared $compared common grid files; missing_in_a=$($missingInA.Count), missing_in_b=$($missingInB.Count)."

if ($ok) {
    Write-Host "PASS"
    exit 0
}

Write-Host "FAIL"
exit 1
