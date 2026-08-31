<#
Interrompe la pipeline HG002. I checkpoint gia validati restano nel volume
Docker e saranno riutilizzati al prossimo avvio.
#>
[CmdletBinding()]
param(
    [string] $TaskName = 'HG002-Parabricks'
)

$ErrorActionPreference = 'Continue'

$task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
if ($task -and $task.State -eq 'Running') {
    Stop-ScheduledTask -TaskName $TaskName
    Write-Host "Attivita pianificata fermata." -ForegroundColor Yellow
}

$containers = docker ps -q --filter 'name=hg002_' 2>$null
if ($containers) {
    docker rm -f $containers 2>$null | Out-Null
    Write-Host "Container della pipeline fermati." -ForegroundColor Yellow
}
else {
    Write-Host "Nessun container HG002 attivo." -ForegroundColor DarkGray
}

Write-Host "I read, i log e i checkpoint validi non sono stati cancellati." -ForegroundColor Green

