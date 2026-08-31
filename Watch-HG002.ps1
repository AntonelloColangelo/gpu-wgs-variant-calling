<#
.SYNOPSIS
    Segue in diretta il log della pipeline HG002 in corso.

.DESCRIPTION
    E' un'operazione di sola lettura. La pipeline gira come attivita pianificata
    di Windows e scrive su file, quindi si possono attaccare quanti lettori si
    vuole, in parallelo, senza disturbarla. Ctrl+C stacca il lettore e non ferma
    l'analisi: per fermarla davvero serve .\Stop-HG002.ps1

    Senza parametri segue il log complessivo della run registrata come corrente
    in logs\current_run.env, e se non ce n'e' una prende la piu recente presente
    su disco.

.PARAMETER Root
    Cartella del progetto. Serve solo se lo script viene eseguito da altrove.

.PARAMETER RunId
    Identificatore della run, per esempio 20260729_232130_full. Se omesso viene
    dedotto da logs\current_run.env, altrimenti si usa la run piu recente.

.PARAMETER Step
    Segue il log di un singolo step invece di quello complessivo. Basta una
    parte del nome: 'fq2bam', 'bqsr', 'haplotype'.

.PARAMETER CurrentStep
    Segue il log dello step in esecuzione in questo momento, letto da
    logs\current_step.env. Utile per vedere l'avanzamento fine di Parabricks.

.PARAMETER List
    Non segue niente: elenca gli step della run con dimensione del log e ora
    dell'ultimo aggiornamento, poi esce.

.PARAMETER Tail
    Quante righe finali mostrare prima di iniziare a seguire. Default 40.

.EXAMPLE
    .\Watch-HG002.ps1
    Il log complessivo della pipeline, in diretta.

.EXAMPLE
    .\Watch-HG002.ps1 -CurrentStep
    Il log dettagliato dello step in esecuzione adesso.

.EXAMPLE
    .\Watch-HG002.ps1 -Step fq2bam -Tail 100
    Le ultime 100 righe di fq2bam, poi prosegue in diretta.

.EXAMPLE
    .\Watch-HG002.ps1 -List
    Quali step esistono nella run corrente.

.LINK
    Stop-HG002.ps1
.LINK
    Get-PipelineStatus.ps1
#>
[CmdletBinding()]
param(
    [string] $Root = 'C:\HG002 Parabricks Experiment',
    [string] $RunId,
    [string] $Step,
    [switch] $CurrentStep,
    [switch] $List,
    [int]    $Tail = 40
)

$ErrorActionPreference = 'Stop'

function Read-EnvFile([string] $Path) {
    $values = @{}
    if (Test-Path -LiteralPath $Path) {
        Get-Content -LiteralPath $Path | ForEach-Object {
            if ($_ -match '^([^=]+)=(.*)$') { $values[$Matches[1]] = $Matches[2] }
        }
    }
    return $values
}

# La run da seguire: quella indicata, altrimenti quella registrata come
# corrente, altrimenti la piu recente presente su disco.
$runsDir = Join-Path $Root 'logs\runs'
if (-not (Test-Path -LiteralPath $runsDir)) {
    throw "Non trovo $runsDir. La pipeline non e' mai stata avviata da questa cartella."
}

$run = Read-EnvFile (Join-Path $Root 'logs\current_run.env')
if (-not $RunId) { $RunId = $run['RUN_ID'] }
if (-not $RunId -or -not (Test-Path -LiteralPath (Join-Path $runsDir $RunId))) {
    $ultima = Get-ChildItem -LiteralPath $runsDir -Directory |
              Sort-Object LastWriteTime -Descending | Select-Object -First 1
    if (-not $ultima) { throw "Nessuna run in $runsDir." }
    $RunId = $ultima.Name
}
$runDir = Join-Path $runsDir $RunId

if ($List) {
    Write-Host "Step disponibili in $RunId :" -ForegroundColor Cyan
    Get-ChildItem -LiteralPath $runDir -Filter '*.log' |
        Where-Object { $_.Name -ne 'pipeline.log' } |
        Sort-Object Name |
        Select-Object @{n='Step';e={$_.BaseName}},
                      @{n='KB';e={[math]::Round($_.Length/1KB,1)}},
                      LastWriteTime |
        Format-Table -AutoSize
    return
}

# Quale file seguire.
if ($CurrentStep) {
    $stato = Read-EnvFile (Join-Path $Root 'logs\current_step.env')
    $chiave = $stato['STEP_KEY']
    if (-not $chiave) { throw "Nessuno step corrente registrato in logs\current_step.env" }
    $Step = $chiave
}

if ($Step) {
    $candidati = Get-ChildItem -LiteralPath $runDir -Filter '*.log' |
                 Where-Object { $_.BaseName -like "*$Step*" -and $_.Name -ne 'pipeline.log' } |
                 Sort-Object LastWriteTime -Descending
    if (-not $candidati) {
        Write-Host "Nessuno step corrisponde a '$Step'. Disponibili:" -ForegroundColor Yellow
        Get-ChildItem -LiteralPath $runDir -Filter '*.log' |
            Where-Object { $_.Name -ne 'pipeline.log' } |
            Sort-Object Name | ForEach-Object { "  $($_.BaseName)" }
        return
    }
    $file = $candidati[0].FullName
} else {
    $file = Join-Path $runDir 'pipeline.log'
}

# Lo step in corso puo non avere ancora scritto il proprio log.
if (-not (Test-Path -LiteralPath $file)) {
    Write-Host "Attendo che $([IO.Path]::GetFileName($file)) venga creato..." -ForegroundColor DarkGray
    while (-not (Test-Path -LiteralPath $file)) { Start-Sleep -Milliseconds 500 }
}

$manifest = Join-Path $runDir 'run_manifest.json'
$riferimento = '(non dichiarato)'
if (Test-Path -LiteralPath $manifest) {
    $riferimento = (Get-Content -LiteralPath $manifest -Raw |
                    ConvertFrom-Json).ingredients.reference
}

$avviata = $run['STARTED']
$trascorso = ''
if ($avviata) {
    try {
        $d = [datetime]::Parse($avviata)
        $t = (Get-Date) - $d
        $trascorso = "  |  in corso da {0:hh\:mm\:ss}" -f $t
    } catch { }
}

Write-Host ''
Write-Host '============================================================' -ForegroundColor Cyan
Write-Host " Run:         $RunId$trascorso"
Write-Host " Prefisso:    $($run['PREFIX'])"
Write-Host " Chiave:      $($run['RUN_KEY'])   Volume: $($run['WORK_VOLUME'])"
Write-Host " Riferimento: $riferimento"
Write-Host " Seguo:       $([IO.Path]::GetFileName($file))"
Write-Host '============================================================' -ForegroundColor Cyan
Write-Host ' Ctrl+C stacca solo questa finestra. Per fermare: .\Stop-HG002.ps1' -ForegroundColor DarkGray
Write-Host ''

Get-Content -LiteralPath $file -Wait -Tail $Tail
