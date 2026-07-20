param(
    [string]$Version = "10.1.3",
    [string]$DownloadId = "105289",
    [string]$CacheDir = "$PSScriptRoot\.cache\msmpi",
    [switch]$ForceDownload,
    [switch]$NoInstall
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RuntimeFile = "msmpisetup.exe"
$SdkFile = "msmpisdk.msi"
$DownloadPage = "https://www.microsoft.com/en-us/download/details.aspx?id=$DownloadId"
$ConfirmationPage = "https://www.microsoft.com/en-us/download/confirmation.aspx?id=$DownloadId"
$KnownLinks = @{
    "10.1.3" = @{
        $RuntimeFile = "https://download.microsoft.com/download/7/2/7/72731ebb-b63c-4170-ade7-836966263a8f/msmpisetup.exe"
        $SdkFile = "https://download.microsoft.com/download/7/2/7/72731ebb-b63c-4170-ade7-836966263a8f/msmpisdk.msi"
    }
}

function Test-IsAdmin {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = [Security.Principal.WindowsPrincipal]::new($identity)
    return $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)
}

function Find-File {
    param([string[]]$Candidates)

    foreach ($candidate in $Candidates) {
        if ($candidate -and (Test-Path -LiteralPath $candidate)) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }

    return $null
}

function Get-RegistryValue {
    param(
        [string[]]$Keys,
        [string[]]$Names
    )

    foreach ($key in $Keys) {
        if (-not (Test-Path -LiteralPath $key)) {
            continue
        }

        $item = Get-ItemProperty -LiteralPath $key -ErrorAction SilentlyContinue
        if (-not $item) {
            continue
        }

        foreach ($name in $Names) {
            $property = $item.PSObject.Properties[$name]
            if ($property -and $property.Value) {
                return [string]$property.Value
            }
        }
    }

    return $null
}

function Get-MsMpiRuntime {
    $fromPath = Get-Command mpiexec.exe -ErrorAction SilentlyContinue
    if ($fromPath) {
        return [pscustomobject]@{
            Installed = $true
            Path = $fromPath.Source
            Source = "PATH"
        }
    }

    $runtimeCandidates = @(
        "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe",
        "C:\Program Files (x86)\Microsoft MPI\Bin\mpiexec.exe"
    )
    if ($env:MSMPI_BIN) {
        $runtimeCandidates = @((Join-Path -Path $env:MSMPI_BIN -ChildPath "mpiexec.exe")) + $runtimeCandidates
    }

    $binFromEnv = Find-File $runtimeCandidates
    if ($binFromEnv) {
        return [pscustomobject]@{
            Installed = $true
            Path = $binFromEnv
            Source = "MSMPI_BIN/default path"
        }
    }

    $installRoot = Get-RegistryValue `
        -Keys @(
            "HKLM:\SOFTWARE\Microsoft\MPI",
            "HKLM:\SOFTWARE\WOW6432Node\Microsoft\MPI"
        ) `
        -Names @("InstallRoot", "Path", "InstallDir")

    if ($installRoot) {
        $mpiexec = Find-File @(
            (Join-Path -Path $installRoot -ChildPath "Bin\mpiexec.exe"),
            (Join-Path -Path $installRoot -ChildPath "mpiexec.exe")
        )
        if ($mpiexec) {
            return [pscustomobject]@{
                Installed = $true
                Path = $mpiexec
                Source = "registry"
            }
        }
    }

    return [pscustomobject]@{
        Installed = $false
        Path = $null
        Source = $null
    }
}

function Get-MsMpiSdk {
    $sdkCandidates = @(
        "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h",
        "C:\Program Files\Microsoft SDKs\MPI\Include\mpi.h"
    )
    if ($env:MSMPI_INC) {
        $sdkCandidates = @((Join-Path -Path $env:MSMPI_INC -ChildPath "mpi.h")) + $sdkCandidates
    }

    $header = Find-File $sdkCandidates

    if ($header) {
        return [pscustomobject]@{
            Installed = $true
            Header = $header
            Source = "MSMPI_INC/default path"
        }
    }

    $sdkRoot = Get-RegistryValue `
        -Keys @(
            "HKLM:\SOFTWARE\Microsoft\Microsoft SDKs\MPI",
            "HKLM:\SOFTWARE\WOW6432Node\Microsoft\Microsoft SDKs\MPI",
            "HKLM:\SOFTWARE\Microsoft\MPI\SDK",
            "HKLM:\SOFTWARE\WOW6432Node\Microsoft\MPI\SDK"
        ) `
        -Names @("InstallationFolder", "InstallRoot", "Path", "InstallDir")

    if ($sdkRoot) {
        $header = Find-File @(
            (Join-Path -Path $sdkRoot -ChildPath "Include\mpi.h"),
            (Join-Path -Path $sdkRoot -ChildPath "mpi.h")
        )
        if ($header) {
            return [pscustomobject]@{
                Installed = $true
                Header = $header
                Source = "registry"
            }
        }
    }

    return [pscustomobject]@{
        Installed = $false
        Header = $null
        Source = $null
    }
}

function Resolve-MsMpiDownloadLinks {
    Write-Host "Resolving MS-MPI $Version download links from Microsoft Download Center..."

    $pages = @($ConfirmationPage, $DownloadPage)
    $links = @{}

    foreach ($page in $pages) {
        try {
            $response = Invoke-WebRequest -Uri $page -UseBasicParsing
            $content = [System.Net.WebUtility]::HtmlDecode($response.Content)

            $dlcMatch = [regex]::Match(
                $content,
                "window\.__DLCDetails__=(?<json>\{.*?\});",
                [System.Text.RegularExpressions.RegexOptions]::Singleline
            )
            if ($dlcMatch.Success) {
                try {
                    $details = $dlcMatch.Groups["json"].Value | ConvertFrom-Json
                    foreach ($fileInfo in $details.dlcDetailsView.downloadFile) {
                        if ($fileInfo.name -in @($RuntimeFile, $SdkFile) -and $fileInfo.url) {
                            $links[$fileInfo.name] = [string]$fileInfo.url
                        }
                    }
                }
                catch {
                    Write-Verbose "Could not parse Microsoft Download Center metadata from $page : $($_.Exception.Message)"
                }
            }

            foreach ($file in @($RuntimeFile, $SdkFile)) {
                if ($links.ContainsKey($file)) {
                    continue
                }

                $pattern = "https://download\.microsoft\.com/[^`"'\s<>]+/$file"
                $match = [regex]::Match($content, $pattern, [System.Text.RegularExpressions.RegexOptions]::IgnoreCase)
                if ($match.Success) {
                    $links[$file] = $match.Value
                }
            }
        }
        catch {
            Write-Verbose "Could not inspect $page : $($_.Exception.Message)"
        }
    }

    if (-not $links.ContainsKey($RuntimeFile) -or -not $links.ContainsKey($SdkFile)) {
        if ($KnownLinks.ContainsKey($Version)) {
            foreach ($file in @($RuntimeFile, $SdkFile)) {
                if (-not $links.ContainsKey($file) -and $KnownLinks[$Version].ContainsKey($file)) {
                    $links[$file] = $KnownLinks[$Version][$file]
                }
            }
        }
    }

    if (-not $links.ContainsKey($RuntimeFile) -or -not $links.ContainsKey($SdkFile)) {
        throw @"
Could not resolve direct MS-MPI installer links automatically.
Open $DownloadPage manually, download $RuntimeFile and $SdkFile into:
  $CacheDir
Then rerun this script.
"@
    }

    return $links
}

function Save-Installer {
    param(
        [string]$FileName,
        [string]$Uri
    )

    New-Item -ItemType Directory -Force -Path $CacheDir | Out-Null
    $target = Join-Path -Path $CacheDir -ChildPath $FileName

    if ((Test-Path -LiteralPath $target) -and -not $ForceDownload) {
        Write-Host "Using cached $FileName"
        return $target
    }

    Write-Host "Downloading $FileName..."
    Invoke-WebRequest -Uri $Uri -OutFile $target -UseBasicParsing
    return $target
}

function Get-InstallerPath {
    param(
        [string]$FileName,
        [hashtable]$Links
    )

    $cached = Join-Path -Path $CacheDir -ChildPath $FileName
    if ((Test-Path -LiteralPath $cached) -and -not $ForceDownload) {
        return (Resolve-Path -LiteralPath $cached).Path
    }

    return Save-Installer -FileName $FileName -Uri $Links[$FileName]
}

function Invoke-Install {
    param(
        [string]$FilePath,
        [string]$Kind
    )

    if ($NoInstall) {
        Write-Host "NoInstall set; skipping $Kind install: $FilePath"
        return
    }

    if (-not (Test-IsAdmin)) {
        throw "Installing $Kind requires an elevated PowerShell session."
    }

    if ($Kind -eq "runtime") {
        Write-Host "Installing MS-MPI runtime..."
        $process = Start-Process -FilePath $FilePath -ArgumentList "/unattend" -Wait -PassThru -WindowStyle Hidden
    }
    elseif ($Kind -eq "sdk") {
        Write-Host "Installing MS-MPI SDK..."
        $arguments = @("/i", "`"$FilePath`"", "/quiet", "/norestart")
        $process = Start-Process -FilePath "msiexec.exe" -ArgumentList $arguments -Wait -PassThru -WindowStyle Hidden
    }
    else {
        throw "Unknown installer kind: $Kind"
    }

    if ($process.ExitCode -ne 0) {
        throw "$Kind installer failed with exit code $($process.ExitCode)."
    }
}

Write-Host "Checking Microsoft MPI $Version..."

$runtime = Get-MsMpiRuntime
$sdk = Get-MsMpiSdk

Write-Host ("Runtime: {0}" -f ($(if ($runtime.Installed) { "found ($($runtime.Path))" } else { "missing" })))
Write-Host ("SDK:     {0}" -f ($(if ($sdk.Installed) { "found ($($sdk.Header))" } else { "missing" })))

if ($runtime.Installed -and $sdk.Installed) {
    Write-Host "MS-MPI runtime and SDK are already installed."
    exit 0
}

$links = $null
if (-not $runtime.Installed -or -not $sdk.Installed) {
    $links = Resolve-MsMpiDownloadLinks
}

if (-not $runtime.Installed) {
    $runtimeInstaller = Get-InstallerPath -FileName $RuntimeFile -Links $links
    Invoke-Install -FilePath $runtimeInstaller -Kind "runtime"
}

if (-not $sdk.Installed) {
    $sdkInstaller = Get-InstallerPath -FileName $SdkFile -Links $links
    Invoke-Install -FilePath $sdkInstaller -Kind "sdk"
}

if ($NoInstall) {
    Write-Host ""
    Write-Host "NoInstall set; installers are available in:"
    Write-Host "  $CacheDir"
    exit 0
}

$runtime = Get-MsMpiRuntime
$sdk = Get-MsMpiSdk

Write-Host ""
Write-Host ("Runtime: {0}" -f ($(if ($runtime.Installed) { "found ($($runtime.Path))" } else { "missing" })))
Write-Host ("SDK:     {0}" -f ($(if ($sdk.Installed) { "found ($($sdk.Header))" } else { "missing" })))

if (-not $runtime.Installed -or -not $sdk.Installed) {
    Write-Host "MS-MPI install may have succeeded, but this shell does not see the updated environment yet."
    Write-Host "Open a new PowerShell session and rerun this script to verify."
    exit 2
}

Write-Host "MS-MPI runtime and SDK are ready."
exit 0
