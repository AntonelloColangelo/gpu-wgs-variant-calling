<#
Mostra lo stato della nuova pipeline Bash.

Esempi:
  .\Get-PipelineStatus.ps1
  .\Get-PipelineStatus.ps1 -Watch
#>
[CmdletBinding()]
param(
    [switch] $Watch,
    [int] $IntervalSec = 30,
    [int] $Tail = 12,
    [string] $Root = 'C:\HG002 Parabricks Experiment',
    [string] $TaskName = 'HG002-Parabricks'
)

$ErrorActionPreference = 'Continue'

$Steps = @(
    @{ Key='1_fastq_integrity'; Label='integrita FASTQ e pairing' },
    @{ Key='2_fq2bam';          Label='fq2bam (allineamento, sort, duplicati)' },
    @{ Key='2b_bam_check';      Label='integrita BAM e read group' },
    @{ Key='3_bqsr';            Label='BQSR separato' },
    @{ Key='3b_bqsr_check';     Label='integrita report BQSR' },
    @{ Key='4_haplotypecaller'; Label='HaplotypeCaller (gVCF)' },
    @{ Key='4b_gvcf_check';     Label='integrita gVCF' },
    @{ Key='5_genotypegvcf';    Label='GenotypeGVCF' },
    @{ Key='5b_vcf_check';      Label='integrita VCF' },
    @{ Key='6_hardfilter';      Label='hard filtering' },
    @{ Key='7_annotation';      Label='annotazione snpEff' },
    @{ Key='8_qc';              Label='QC completo' },
    @{ Key='8b_export';         Label='esportazione su Windows' },
    @{ Key='9_report';          Label='report HTML' }
)

function Read-EnvFile([string] $Path) {
    $values = @{}
    if (Test-Path -LiteralPath $Path) {
        Get-Content -LiteralPath $Path | ForEach-Object {
            if ($_ -match '^([^=]+)=(.*)$') { $values[$Matches[1]] = $Matches[2] }
        }
    }
    return $values
}

function Show-Status {
    if ($Watch) { Clear-Host }
    Write-Host "=== Pipeline HG002 - $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') ===" -ForegroundColor Cyan

    $run = Read-EnvFile (Join-Path $Root 'logs\current_run.env')
    $state = Read-EnvFile (Join-Path $Root 'logs\current_step.env')
    $prefix = $run['PREFIX']
    $mode = $run['MODE']
    $volume = $run['WORK_VOLUME']
    $runId = $run['RUN_ID']
    $runDir = if ($runId) { Join-Path $Root "logs\runs\$runId" } else { $null }

    Write-Host "`n--- Run selezionata ---" -ForegroundColor Yellow
    if ($prefix) {
        "  Modalita: $mode"
        "  Prefisso: $prefix"
        "  Run ID:   $runId"
        "  Volume:   $volume"
    } else {
        Write-Host "  Nessuna nuova run registrata." -ForegroundColor DarkGray
    }

    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    if ($task) {
        $info = Get-ScheduledTaskInfo -TaskName $TaskName -ErrorAction SilentlyContinue
        "  Task:     $($task.State) (ultimo codice: $($info.LastTaskResult))"
    }

    $containerRows = docker ps --filter 'name=hg002_' --format '{{.Names}}|{{.Status}}' 2>$null
    $containers = @()
    Write-Host "`n--- Container attivi ---" -ForegroundColor Yellow
    if ($containerRows) {
        $containerRows | ForEach-Object {
            $parts = $_ -split '\|', 2
            $containers += $parts[0]
            "  {0,-31} {1}" -f $parts[0], $parts[1]
        }
    } else {
        Write-Host "  nessuno" -ForegroundColor DarkGray
    }

    $timingFile = if ($runDir -and $prefix) {
        Join-Path $runDir "$prefix.step_times.tsv"
    } else { $null }
    $timings = if ($timingFile -and (Test-Path -LiteralPath $timingFile)) {
        @(Import-Csv -LiteralPath $timingFile -Delimiter "`t")
    } else { @() }

    Write-Host "`n--- Step ---" -ForegroundColor Yellow
    foreach ($step in $Steps) {
        $row = $timings | Where-Object { $_.Step -eq $step.Key } | Select-Object -Last 1
        if ($row -and $row.Status -in @('PASS','SKIPPED_VALID')) {
            Write-Host ("  [x] {0}" -f $step.Label) -ForegroundColor DarkGreen
        } elseif ($row -and $row.Status -eq 'FAILED') {
            Write-Host ("  [!] {0}  FALLITO" -f $step.Label) -ForegroundColor Red
        } elseif ($state['STEP_KEY'] -eq $step.Key -and $state['STEP_STATUS'] -eq 'RUNNING') {
            Write-Host ("  [>] {0}  in corso" -f $step.Label) -ForegroundColor Green
        } else {
            "  [ ] {0}" -f $step.Label
        }
    }
    if ($state['STEP_KEY'] -eq 'complete' -and $state['STEP_STATUS'] -eq 'PASS') {
        Write-Host "  PIPELINE COMPLETATA E VALIDATA" -ForegroundColor Green
    } elseif ($state['STEP_STATUS'] -eq 'FAILED') {
        Write-Host "  Ultimo stato: FALLITO - consultare il log indicato sotto." -ForegroundColor Red
    }

    # Progresso numerico di fq2bam.
    $fqLog = if ($runDir) { Join-Path $runDir '2_fq2bam.log' } else { $null }
    if ($fqLog -and (Test-Path -LiteralPath $fqLog)) {
        $fqRecord = $timings | Where-Object { $_.Step -eq '2_fq2bam' } | Select-Object -Last 1
        $fqCompleted = $fqRecord -and $fqRecord.Status -in @('PASS','SKIPPED_VALID')
        $samples = @(
            Get-Content -LiteralPath $fqLog -Tail 600 | ForEach-Object {
                if ($_ -match 'pool:\s+\d+\s+(\d+)\s+bases/GPU/minute:\s+([0-9.]+)') {
                    [pscustomobject]@{ Bases=[double]$Matches[1]; Speed=[double]$Matches[2] }
                }
            }
        )
        if ($samples.Count -gt 0) {
            $last = $samples[-1]
            $recent = @($samples | Select-Object -Last 6 | ForEach-Object Speed)
            $speed = ($recent | Measure-Object -Average).Average
            $pairs = if ($mode -eq 'smoke') { 1000000 } else { 474384500 }
            $totalBases = [double]$pairs * 2 * 151
            $pct = [math]::Min(100, 100 * $last.Bases / $totalBases)
            Write-Host "`n--- Progresso fq2bam ---" -ForegroundColor Yellow
            if ($fqCompleted) {
                "  Stato: allineamento, sort e marcatura duplicati completati"
                if ($fqRecord.Seconds) { "  Durata fq2bam: $($fqRecord.Seconds) secondi" }
            } elseif ($speed -gt 0) {
                "  Allineamento: {0:N1} %" -f $pct
                "  Basi:         {0:N2} / {1:N2} miliardi" -f ($last.Bases/1e9), ($totalBases/1e9)
                "  Velocita':     {0:N2} miliardi di basi/min" -f ($speed/1e9)
                $minutes = ($totalBases - $last.Bases) / $speed
                "  Fine stimata allineamento: {0:HH:mm:ss} (circa {1:N0} min)" -f (Get-Date).AddMinutes($minutes), $minutes
            } else {
                "  Stato: inizializzazione delle code di allineamento"
            }
        }
    }

    # Progresso HaplotypeCaller dal ProgressMeter.
    $hcLog = if ($runDir) { Join-Path $runDir '4_haplotypecaller.log' } else { $null }
    if ($hcLog -and (Test-Path -LiteralPath $hcLog)) {
        $hc = @(
            Get-Content -LiteralPath $hcLog -Tail 1200 | ForEach-Object {
                if ($_ -match '\]\s+(.+):(\d+)\s+([0-9.]+)\s+(\d+)\s+([0-9.]+)\s*$') {
                    [pscustomobject]@{
                        Contig=$Matches[1]; Pos=[int64]$Matches[2]
                        Minutes=[double]$Matches[3]; Regions=[int64]$Matches[4]
                        Rate=[double]$Matches[5]
                    }
                }
            }
        )
        if ($hc.Count -gt 0) {
            $last = $hc[-1]
            Write-Host "`n--- Progresso HaplotypeCaller ---" -ForegroundColor Yellow
            "  Posizione:    $($last.Contig):$($last.Pos)"
            "  Trascorso:    {0:N1} min" -f $last.Minutes
            "  Regioni:      {0:N0}" -f $last.Regions
            "  Velocita':    {0:N0} regioni/min" -f $last.Rate

            $progressPct = $null
            if ($mode -eq 'smoke' -and $last.Contig -eq 'chr20') {
                $progressPct = 100 * [math]::Min($last.Pos, 1000000) / 1000000
            } else {
                $fai = Join-Path $Root 'ref\Homo_sapiens_assembly38.fasta.fai'
                if (Test-Path -LiteralPath $fai) {
                    [double]$before = 0
                    [double]$total = 0
                    [double]$currentLength = 0
                    $found = $false
                    Get-Content -LiteralPath $fai | ForEach-Object {
                        $f = $_ -split "`t"
                        $length = [double]$f[1]
                        if (-not $found -and $f[0] -eq $last.Contig) {
                            $before = $total
                            $currentLength = $length
                            $found = $true
                        }
                        $total += $length
                    }
                    if ($found -and $total -gt 0) {
                        $progressPct = 100 * ($before + [math]::Min($last.Pos, $currentLength)) / $total
                    }
                }
            }
            if ($null -ne $progressPct) {
                "  Genoma/intervallo percorso: {0:N1} %" -f $progressPct
                if ($progressPct -ge 1 -and $progressPct -lt 99.9) {
                    $left = $last.Minutes * (100-$progressPct) / $progressPct
                    "  Fine stimata: {0:HH:mm:ss} (circa {1:N0} min)" -f (Get-Date).AddMinutes($left), $left
                }
            }
        }
    }

    # Dimensioni dei checkpoint nel volume Linux.
    if ($volume -and $prefix) {
        $volumeInfo = docker run --rm --mount "type=volume,source=$volume,target=/work,readonly" `
            alpine:3.20 sh -c "for f in /work/$prefix.bam /work/$prefix.partial.bam /work/$prefix.g.vcf.gz /work/$prefix.vcf.gz; do if [ -e `"`$f`" ]; then stat -c '%n|%s' `"`$f`"; fi; done" 2>$null
        if ($volumeInfo) {
            Write-Host "`n--- Checkpoint nel volume Linux ---" -ForegroundColor Yellow
            $volumeInfo | ForEach-Object {
                $p = $_ -split '\|'
                "  {0,-45} {1,9:N2} GiB" -f ([IO.Path]::GetFileName($p[0])), ([double]$p[1]/1GB)
            }
        }
    }

    Write-Host "`n--- Risorse ---" -ForegroundColor Yellow
    $gpu = nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw --format=csv,noheader,nounits 2>$null
    if ($gpu) {
        $g = $gpu -split ',' | ForEach-Object { $_.Trim() }
        "  GPU:  $($g[0]) % | VRAM $($g[1])/$($g[2]) MiB | $($g[3]) C | $($g[4]) W"
    }
    if ($containers.Count -gt 0) {
        $stats = docker stats --no-stream --format '{{.Name}}|{{.CPUPerc}}|{{.MemUsage}}' $containers[0] 2>$null
        if ($stats) {
            $s = $stats -split '\|'
            "  Container: CPU $($s[1]) | RAM $($s[2])"
        }
    }
    $free = [math]::Round((Get-PSDrive -Name $Root[0]).Free/1GB, 1)
    "  Disco $($Root[0]): $free GB liberi"

    $currentLog = if ($runDir -and $state['STEP_KEY']) {
        Join-Path $runDir "$($state['STEP_KEY']).log"
    } else { $null }
    if (-not $currentLog -or -not (Test-Path -LiteralPath $currentLog)) {
        $currentLog = if ($runDir) {
            Get-ChildItem -LiteralPath $runDir -Filter '*.log' -ErrorAction SilentlyContinue |
                Sort-Object LastWriteTime | Select-Object -Last 1 -ExpandProperty FullName
        } else { $null }
    }
    if ($currentLog -and (Test-Path -LiteralPath $currentLog)) {
        Write-Host "`n--- $(Split-Path $currentLog -Leaf) (ultime $Tail righe) ---" -ForegroundColor Yellow
        Get-Content -LiteralPath $currentLog -Tail $Tail | ForEach-Object { "  $_" }
    }
}

if ($Watch) {
    while ($true) {
        Show-Status
        Write-Host "`n(aggiornamento ogni ${IntervalSec}s - Ctrl+C per uscire)" -ForegroundColor DarkGray
        Start-Sleep -Seconds $IntervalSec
    }
} else {
    Show-Status
}
