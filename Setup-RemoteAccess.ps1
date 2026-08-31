<#
.SYNOPSIS
    Configura l'accesso remoto al PC di analisi: SSH + Tailscale.

.DESCRIPTION
    DA ESEGUIRE COME AMMINISTRATORE.
    Tasto destro su PowerShell -> "Esegui come amministratore", poi:

        cd "C:\HG002 Parabricks Experiment"
        .\Setup-RemoteAccess.ps1

    Cosa fa:
      1. avvia OpenSSH Server (gia' installato su questa macchina) e lo mette
         in avvio automatico
      2. apre la porta 22 sul firewall
      3. imposta PowerShell come shell predefinita di SSH al posto di cmd
      4. prepara authorized_keys con i permessi corretti
      5. installa Tailscale

    Tailscale serve perche' il portatile sara' su reti diverse: crea una rete
    privata fra le tue macchine senza aprire porte sul router e senza esporre
    SSH su internet.

.PARAMETER PublicKey
    Chiave pubblica SSH del portatile Linux (contenuto di ~/.ssh/id_ed25519.pub).
    Se omessa, resta attiva l'autenticazione con password.

.PARAMETER SkipTailscale
    Salta l'installazione di Tailscale (se usi gia' un'altra VPN).
#>
[CmdletBinding()]
param(
    [string] $PublicKey,
    [switch] $SkipTailscale
)

$ErrorActionPreference = 'Stop'

# --- 0) controllo privilegi ---------------------------------------------------
$id = [Security.Principal.WindowsIdentity]::GetCurrent()
if (-not (New-Object Security.Principal.WindowsPrincipal($id)).IsInRole(
          [Security.Principal.WindowsBuiltInRole]::Administrator)) {
    throw "Questo script richiede privilegi di amministratore. Riapri PowerShell con 'Esegui come amministratore'."
}
Write-Host "Eseguito come amministratore ($($id.Name))." -ForegroundColor Green

# --- 1) OpenSSH Server --------------------------------------------------------
Write-Host "`n[1/5] OpenSSH Server" -ForegroundColor Cyan
$cap = Get-WindowsCapability -Online -Name 'OpenSSH.Server*'
if ($cap.State -ne 'Installed') {
    Write-Host "  installo la funzionalita' ..."
    Add-WindowsCapability -Online -Name $cap.Name | Out-Null
} else {
    Write-Host "  gia' installata"
}
Set-Service -Name sshd -StartupType Automatic
Start-Service -Name sshd
Write-Host "  servizio sshd: $((Get-Service sshd).Status), avvio automatico" -ForegroundColor Green

# --- 2) firewall --------------------------------------------------------------
Write-Host "`n[2/5] Firewall" -ForegroundColor Cyan
if (-not (Get-NetFirewallRule -Name 'sshd' -ErrorAction SilentlyContinue)) {
    New-NetFirewallRule -Name 'sshd' -DisplayName 'OpenSSH Server (sshd)' `
        -Enabled True -Direction Inbound -Protocol TCP -Action Allow -LocalPort 22 | Out-Null
    Write-Host "  regola creata per la porta 22" -ForegroundColor Green
} else {
    Write-Host "  regola gia' presente"
}

# --- 3) shell predefinita -----------------------------------------------------
# Senza questo SSH apre cmd.exe, dove gli script .ps1 non partono direttamente.
Write-Host "`n[3/5] Shell predefinita di SSH" -ForegroundColor Cyan
New-ItemProperty -Path 'HKLM:\SOFTWARE\OpenSSH' -Name DefaultShell `
    -Value 'C:\Windows\System32\WindowsPowerShell\v1.0\powershell.exe' `
    -PropertyType String -Force | Out-Null
Write-Host "  impostata su PowerShell" -ForegroundColor Green

# --- 4) chiave pubblica -------------------------------------------------------
Write-Host "`n[4/5] Autenticazione a chiave" -ForegroundColor Cyan
if ($PublicKey) {
    # Windows OpenSSH tratta diversamente gli utenti amministratori: per loro
    # legge una chiave globale, non quella nella home.
    $isAdminUser = ([Security.Principal.WindowsPrincipal]$id).IsInRole('Administrators')
    if ($isAdminUser) {
        $keyFile = 'C:\ProgramData\ssh\administrators_authorized_keys'
    } else {
        $sshDir  = Join-Path $env:USERPROFILE '.ssh'
        New-Item -ItemType Directory -Force -Path $sshDir | Out-Null
        $keyFile = Join-Path $sshDir 'authorized_keys'
    }

    Add-Content -Path $keyFile -Value $PublicKey.Trim() -Encoding ASCII
    # ACL restrittive: sshd rifiuta il file se e' scrivibile da altri
    icacls $keyFile /inheritance:r /grant 'SYSTEM:F' /grant 'Administrators:F' | Out-Null
    Write-Host "  chiave aggiunta in $keyFile" -ForegroundColor Green
} else {
    Write-Host "  nessuna chiave fornita: resta attiva la password" -ForegroundColor Yellow
    Write-Host "  per aggiungerla poi: .\Setup-RemoteAccess.ps1 -PublicKey 'ssh-ed25519 AAAA...'"
}

# --- 5) Tailscale -------------------------------------------------------------
Write-Host "`n[5/5] Tailscale" -ForegroundColor Cyan
if ($SkipTailscale) {
    Write-Host "  saltato su richiesta"
} elseif (Get-Command tailscale -ErrorAction SilentlyContinue) {
    Write-Host "  gia' installato: $(tailscale version | Select-Object -First 1)"
} else {
    $msi = Join-Path $env:TEMP 'tailscale-setup.msi'
    Write-Host "  scarico l'installer ..."
    curl.exe -L --fail --no-progress-meter -o $msi 'https://pkgs.tailscale.com/stable/tailscale-setup-latest-amd64.msi'
    if ($LASTEXITCODE -ne 0) { throw "Download di Tailscale fallito (exit $LASTEXITCODE)" }
    Write-Host "  installo ..."
    Start-Process msiexec.exe -ArgumentList '/i', "`"$msi`"", '/quiet', '/norestart' -Wait
    Write-Host "  installato" -ForegroundColor Green
}

# --- riepilogo ----------------------------------------------------------------
Write-Host "`n=== Fatto. Passi rimanenti, da fare a mano ===" -ForegroundColor Cyan
Write-Host @"

  1. Autentica Tailscale su QUESTO pc (apre il browser):

         tailscale up

     poi annota il nome/IP che ti assegna:

         tailscale ip -4

  2. Sul portatile Linux Mint:

         curl -fsSL https://tailscale.com/install.sh | sh
         sudo tailscale up

  3. Collegati dal portatile:

         ssh $env:USERNAME@<ip-o-nome-tailscale>

  4. Una volta dentro:

         cd 'C:\HG002 Parabricks Experiment'
         .\Get-PipelineStatus.ps1

"@
