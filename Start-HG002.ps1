<#
Avvia la pipeline Bash in Ubuntu WSL2 come attivita pianificata.
La finestra PowerShell puo poi essere chiusa; non fare logout da Windows.

Esempi:
  .\Start-HG002.ps1 -PreflightOnly
  .\Start-HG002.ps1 -SmokeTest
  .\Start-HG002.ps1
#>
[CmdletBinding()]
param(
    [switch] $SmokeTest,
    [switch] $PreflightOnly,
    [string] $Root = 'C:\HG002 Parabricks Experiment',
    [string] $Distro = 'Ubuntu-24.04',
    [string] $TaskName = 'HG002-Parabricks'
)

$ErrorActionPreference = 'Stop'

$task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
if ($task -and $task.State -eq 'Running') {
    throw "La pipeline e' gia in esecuzione. Usa .\Stop-HG002.ps1 solo se vuoi interromperla."
}

$fullRoot = [IO.Path]::GetFullPath($Root)
if ($fullRoot -notmatch '^([A-Za-z]):\\(.*)$') {
    throw "Root deve essere un percorso Windows locale, per esempio C:\HG002 Parabricks Experiment."
}
$drive = $Matches[1].ToLowerInvariant()
$tail = $Matches[2] -replace '\\','/'
$wslRoot = "/mnt/$drive/$tail"

& wsl.exe -d $Distro -- test -f "$wslRoot/run_parabricks_hg002.sh"
if ($LASTEXITCODE -ne 0) {
    throw "Ubuntu WSL non vede lo script in $wslRoot."
}

$pipelineArgs = @()
if ($SmokeTest) { $pipelineArgs += '--smoke' }
if ($PreflightOnly) { $pipelineArgs += '--preflight' }
$extra = $pipelineArgs -join ' '

$arguments = "-d $Distro --cd `"$wslRoot`" -- /usr/bin/bash ./run_parabricks_hg002.sh $extra"
$action = New-ScheduledTaskAction `
    -Execute (Join-Path $env:WINDIR 'System32\wsl.exe') `
    -Argument $arguments `
    -WorkingDirectory $Root
$principal = New-ScheduledTaskPrincipal -UserId $env:USERNAME -LogonType Interactive
$settings = New-ScheduledTaskSettingsSet `
    -AllowStartIfOnBatteries `
    -DontStopIfGoingOnBatteries `
    -ExecutionTimeLimit ([TimeSpan]::Zero) `
    -MultipleInstances IgnoreNew

Register-ScheduledTask -TaskName $TaskName -Action $action `
    -Principal $principal -Settings $settings -Force | Out-Null
Start-ScheduledTask -TaskName $TaskName
Start-Sleep -Seconds 2

$state = (Get-ScheduledTask -TaskName $TaskName).State
Write-Host "Pipeline Bash avviata: $state" -ForegroundColor Green
Write-Host "Modalita: $(if ($SmokeTest) {'smoke test'} else {'40x completa'})"
Write-Host "Monitor:  .\Get-PipelineStatus.ps1 -Watch"
Write-Host "Stop:     .\Stop-HG002.ps1"
