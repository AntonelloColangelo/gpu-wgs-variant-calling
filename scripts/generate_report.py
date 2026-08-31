#!/usr/bin/env python3
# Report finale della pipeline HG002.
"""Validazione FASTQ e report automatico per la pipeline HG002 Parabricks.

Usa esclusivamente la libreria standard Python per poter girare dentro
bioinfo-codeserver senza installare dipendenze.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import os
import re
from collections import Counter
from datetime import datetime
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", errors="replace"
    )


def normalized_read_id(header: str) -> str:
    token = header.strip().split()[0]
    if token.startswith("@"):
        token = token[1:]
    return re.sub(r"/[12]$", "", token)


def validate_fastq(args: argparse.Namespace) -> None:
    paths = [Path(args.r1), Path(args.r2)]
    result = {
        "r1": str(paths[0]),
        "r2": str(paths[1]),
        "r1_bytes": paths[0].stat().st_size,
        "r2_bytes": paths[1].stat().st_size,
        "checked_pairs": 0,
        "full_gzip_integrity_checked": bool(args.full_gzip_checked),
        "total_pairs": args.total_pairs,
        "paired_ids_match": True,
        "fastq_structure_valid": True,
        "read_length_min": None,
        "read_length_max": None,
        "status": "PASS",
    }
    lengths: list[int] = []
    with open_text(paths[0]) as r1, open_text(paths[1]) as r2:
        for pair_index in range(args.pairs):
            records = []
            for stream in (r1, r2):
                record = [stream.readline() for _ in range(4)]
                if not record[0]:
                    records = []
                    break
                if any(line == "" for line in record) or not record[0].startswith("@") or not record[2].startswith("+"):
                    result["fastq_structure_valid"] = False
                    result["status"] = "FAIL"
                    raise SystemExit(f"FASTQ non valido alla coppia {pair_index + 1}")
                if len(record[1].rstrip()) != len(record[3].rstrip()):
                    result["fastq_structure_valid"] = False
                    result["status"] = "FAIL"
                    raise SystemExit(f"Sequenza e qualità hanno lunghezze diverse alla coppia {pair_index + 1}")
                records.append(record)
            if not records:
                break
            if normalized_read_id(records[0][0]) != normalized_read_id(records[1][0]):
                result["paired_ids_match"] = False
                result["status"] = "FAIL"
                raise SystemExit(f"Identificativi R1/R2 non corrispondenti alla coppia {pair_index + 1}")
            lengths.extend((len(records[0][1].rstrip()), len(records[1][1].rstrip())))
            result["checked_pairs"] += 1
    if not result["checked_pairs"]:
        raise SystemExit("FASTQ vuoti")
    result["read_length_min"] = min(lengths)
    result["read_length_max"] = max(lengths)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result))


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_flagstat(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    if not path.exists():
        return values
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "total": r"^(\d+) \+ \d+ in total",
        "primary": r"^(\d+) \+ \d+ primary$",
        "mapped": r"^(\d+) \+ \d+ mapped \(([0-9.]+)%",
        "properly_paired": r"^(\d+) \+ \d+ properly paired \(([0-9.]+)%",
        "singletons": r"^(\d+) \+ \d+ singletons \(([0-9.]+)%",
        "duplicates": r"^(\d+) \+ \d+ duplicates$",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            values[key] = int(match.group(1))
            if len(match.groups()) > 1:
                values[f"{key}_pct"] = float(match.group(2))
    if values.get("primary") and "duplicates" in values:
        values["duplicates_pct"] = 100.0 * values["duplicates"] / values["primary"]
    return values


def parse_duplicate_metrics(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    lines = [line.rstrip("\n") for line in path.open(encoding="utf-8", errors="replace")]
    for index, line in enumerate(lines):
        if line.startswith("LIBRARY") and index + 1 < len(lines):
            keys = line.split("\t")
            vals = lines[index + 1].split("\t")
            return dict(zip(keys, vals))
    return {}


def variant_stats(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "records": 0,
        "snp_alleles": 0,
        "indel_complex_alleles": 0,
        "multiallelic": 0,
        "transitions": 0,
        "transversions": 0,
        "filters": Counter(),
        "called_genotypes": 0,
        "het": 0,
        "hom_alt": 0,
        "dp_sum": 0,
        "dp_count": 0,
        "sample": "",
    }
    if not path.exists():
        return result
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                result["sample"] = fields[9] if len(fields) > 9 else ""
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 8:
                continue
            result["records"] += 1
            ref, alts, filt = fields[3], fields[4].split(","), fields[6]
            result["filters"][filt] += 1
            if len(alts) > 1:
                result["multiallelic"] += 1
            for alt in alts:
                if len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                    result["snp_alleles"] += 1
                    if (ref, alt) in transitions:
                        result["transitions"] += 1
                    else:
                        result["transversions"] += 1
                else:
                    result["indel_complex_alleles"] += 1
            if len(fields) > 9:
                fmt = fields[8].split(":")
                sample = fields[9].split(":")
                data = dict(zip(fmt, sample))
                gt = data.get("GT", "./.").replace("|", "/")
                if gt not in (".", "./."):
                    result["called_genotypes"] += 1
                    alleles = gt.split("/")
                    if len(set(alleles)) > 1:
                        result["het"] += 1
                    elif alleles[0] != "0":
                        result["hom_alt"] += 1
                dp = data.get("DP", ".")
                if dp.isdigit():
                    result["dp_sum"] += int(dp)
                    result["dp_count"] += 1
    tv = int(result["transversions"])
    result["titv"] = round(int(result["transitions"]) / tv, 3) if tv else None
    result["mean_dp"] = round(int(result["dp_sum"]) / int(result["dp_count"]), 2) if result["dp_count"] else None
    result["filters"] = dict(result["filters"])
    return result


def parse_timing(path: Path) -> list[dict[str, str]]:
    return read_tsv(path)


def pct(value: object, digits: int = 2) -> str:
    if value is None:
        return "n/d"
    return f"{float(value):.{digits}f}%".replace(".", ",")


def number(value: object, digits: int = 0) -> str:
    if value is None:
        return "n/d"
    if digits:
        return f"{float(value):,.{digits}f}".replace(",", "X").replace(".", ",").replace("X", ".")
    return f"{int(value):,}".replace(",", ".")


def reference_label() -> str:
    """Nome del riferimento usato dalla run, letto dall'ambiente.

    Lo passa run_tools in pipeline_functions.sh. Scriverlo a mano nel report
    significherebbe dichiarare un riferimento che potrebbe non essere quello
    con cui le varianti sono state chiamate.
    """
    reference = os.environ.get("REFERENCE", "")
    return Path(reference).name if reference else "non dichiarato"


def bar_chart(items: list[tuple[str, float]], maximum: float | None = None, suffix: str = "") -> str:
    if not items:
        return "<p>Dati non disponibili.</p>"
    max_value = maximum or max(v for _, v in items) or 1
    rows = []
    for label, value in items:
        width = max(0.5, min(100.0, 100 * value / max_value))
        rows.append(
            f'<div class="bar-row"><span>{html.escape(label)}</span><div class="track">'
            f'<i style="width:{width:.2f}%"></i></div><b>{value:,.2f}{suffix}</b></div>'
        )
    return "".join(rows)


def table(headers: list[str], rows: list[list[object]]) -> str:
    head = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row) + "</tr>" for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def generate_report(args: argparse.Namespace) -> None:
    root = Path(args.root)
    out = root / "output"
    reports = root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix
    mode_label = "smoke test (chr20:1–1.000.000)" if args.mode == "smoke" else "WGS NovaSeq PCR-free 40×"

    required = [
        out / f"{prefix}.bam",
        out / f"{prefix}.bam.bai",
        out / f"{prefix}.duplicate_metrics.txt",
        out / f"{prefix}.recal.txt",
        out / f"{prefix}.g.vcf.gz",
        out / f"{prefix}.g.vcf.gz.tbi",
        out / f"{prefix}.vcf.gz",
        out / f"{prefix}.vcf.gz.tbi",
        out / f"{prefix}.hardfiltered.vcf.gz",
        out / f"{prefix}.hardfiltered.vcf.gz.tbi",
        out / f"{prefix}_annotated.vcf.gz",
        out / f"{prefix}_annotated.vcf.gz.tbi",
        out / f"{prefix}_high_impact_variants.tsv",
        out / f"{prefix}.annotation_qc.tsv",
        out / f"{prefix}.flagstat.txt",
        out / f"{prefix}.coverage.tsv",
        out / f"{prefix}.wgs_metrics.txt",
        out / f"{prefix}.depth_thresholds.tsv",
        out / f"{prefix}.step_times.tsv",
        out / f"{prefix}.fastq_validation.json",
    ]
    missing = [str(path) for path in required if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise SystemExit("Report non generato: output obbligatori mancanti o vuoti:\n- " + "\n- ".join(missing))

    flagstat = parse_flagstat(out / f"{prefix}.flagstat.txt")
    duplicate = parse_duplicate_metrics(out / f"{prefix}.duplicate_metrics.txt")
    coverage = read_tsv(out / f"{prefix}.coverage.tsv")
    depth_rows = read_tsv(out / f"{prefix}.depth_thresholds.tsv")
    timing = parse_timing(out / f"{prefix}.step_times.tsv")
    raw = variant_stats(out / f"{prefix}.vcf.gz")
    filtered = variant_stats(out / f"{prefix}.hardfiltered.vcf.gz")
    annotated = variant_stats(out / f"{prefix}_annotated.vcf.gz")
    fastq_validation_path = out / f"{prefix}.fastq_validation.json"
    fastq_validation = json.loads(fastq_validation_path.read_text(encoding="utf-8"))
    annotation_qc_rows = read_tsv(out / f"{prefix}.annotation_qc.tsv")
    annotation_qc = {row.get("metric", ""): row.get("value", "") for row in annotation_qc_rows}

    primary = [row for row in coverage if re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)", row.get("#rname", ""))]
    autosomes = [row for row in primary if re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])", row.get("#rname", ""))]
    mean_depth = (
        sum(float(row["meandepth"]) * int(row["endpos"]) for row in primary)
        / sum(int(row["endpos"]) for row in primary)
        if primary
        else None
    )
    mean_breadth = (
        sum(float(row["coverage"]) * int(row["endpos"]) for row in primary)
        / sum(int(row["endpos"]) for row in primary)
        if primary
        else None
    )
    rates = []
    for row in primary:
        length_mb = int(row["endpos"]) / 1_000_000
        rates.append((row["#rname"], int(row["numreads"]) / length_mb if length_mb else 0))
    auto_rates = [value for chrom, value in rates if chrom not in ("chrX", "chrY")]
    auto_mean = sum(auto_rates) / len(auto_rates) if auto_rates else None
    auto_cv = (
        (sum((value - auto_mean) ** 2 for value in auto_rates) / len(auto_rates)) ** 0.5 / auto_mean * 100
        if auto_rates and auto_mean
        else None
    )

    depth = depth_rows[0] if depth_rows else {}
    timing_rows = []
    for row in timing:
        seconds = float(row.get("Seconds", 0) or 0)
        timing_rows.append([row.get("Step", ""), f"{seconds / 60:.1f}", row.get("Status", "")])

    filters = raw.get("filters", {})
    raw_pass = sum(v for k, v in filters.items() if k in ("PASS", "."))
    hard_filters = filtered.get("filters", {})
    filtered_pass = sum(v for k, v in hard_filters.items() if k in ("PASS", "."))
    total_filtered = int(filtered.get("records", 0))
    pass_pct = 100 * filtered_pass / total_filtered if total_filtered else None

    ann_errors = int(annotation_qc.get("chromosome_not_found", 0) or 0)
    ann_log = out / f"{prefix}.snpeff.log"
    if not ann_errors and ann_log.exists():
        ann_text = ann_log.read_text(encoding="utf-8", errors="replace")
        match = re.search(r"ERROR_CHROMOSOME_NOT_FOUND\s+(\d+)", ann_text)
        ann_errors = int(match.group(1)) if match else 0

    technical_errors = []
    if fastq_validation.get("status") != "PASS":
        technical_errors.append("controllo FASTQ non superato")
    if not flagstat.get("total"):
        technical_errors.append("flagstat vuoto")
    if not depth:
        technical_errors.append("metriche di copertura mancanti")
    if int(raw.get("records", 0)) != int(filtered.get("records", -1)):
        technical_errors.append("record persi durante l'hard filtering")
    if int(filtered.get("records", 0)) != int(annotated.get("records", -1)):
        technical_errors.append("record persi durante l'annotazione")
    if any(row.get("Status") == "FAILED" for row in timing):
        technical_errors.append("uno step risulta fallito")
    if technical_errors:
        raise SystemExit("Report non generato: " + "; ".join(technical_errors))

    status = "COMPLETATO CON LIMITI DI ANNOTAZIONE" if ann_errors else "COMPLETATO"
    chrom_chart = bar_chart(rates, suffix="")
    depth_chart = bar_chart(
        [
            ("≥1×", float(depth.get("pct_ge1", 0) or 0)),
            ("≥10×", float(depth.get("pct_ge10", 0) or 0)),
            ("≥20×", float(depth.get("pct_ge20", 0) or 0)),
            ("≥30×", float(depth.get("pct_ge30", 0) or 0)),
        ],
        maximum=100,
        suffix="%",
    )
    runtime_chart = bar_chart(
        [(row.get("Step", ""), float(row.get("Seconds", 0) or 0) / 60) for row in timing],
        suffix=" min",
    )

    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sample": args.sample,
        "prefix": prefix,
        "mode": args.mode,
        "status": status,
        "flagstat": flagstat,
        "duplicate_metrics": duplicate,
        "mean_depth_primary": mean_depth,
        "breadth_primary": mean_breadth,
        "autosomal_read_per_mb_cv_pct": auto_cv,
        "depth_thresholds": depth,
        "raw_vcf": raw,
        "hardfiltered_vcf": filtered,
        "annotated_vcf": annotated,
        "pass_high_records": int(annotation_qc.get("pass_high_records", 0) or 0),
        "annotation_chromosome_not_found": ann_errors,
        "giab_benchmark_performed": False,
    }
    (reports / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    css = """
    :root{--blue:#22577a;--dark:#173b53;--gold:#b78318;--pale:#f3f7f9;--line:#d8e1e7;--ink:#1e2933}
    *{box-sizing:border-box}body{font:15px/1.5 Arial,sans-serif;color:var(--ink);margin:0;background:#eef3f6}
    main{max-width:1120px;margin:32px auto;background:white;padding:48px 56px;box-shadow:0 5px 25px #2335}
    h1{font-size:34px;color:var(--dark);margin:0 0 8px}h2{color:var(--blue);border-bottom:2px solid var(--line);padding-bottom:8px;margin-top:38px}
    h3{color:var(--dark)}.subtitle{color:#66737d;font-size:18px}.verdict{padding:16px 20px;background:#dfedf5;border-left:7px solid var(--blue);font-weight:bold;margin:24px 0}
    .warning{background:#f8edcf;border-left-color:var(--gold)}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}
    .card{background:var(--pale);padding:16px;border-radius:6px}.card b{display:block;font-size:24px;color:var(--dark)}
    table{border-collapse:collapse;width:100%;margin:12px 0 22px}th{background:var(--dark);color:white;text-align:left}th,td{padding:8px 10px;border:1px solid var(--line)}
    .bar-row{display:grid;grid-template-columns:150px 1fr 110px;gap:10px;align-items:center;margin:7px 0}.track{height:18px;background:#e8eef2;border-radius:4px;overflow:hidden}
    .track i{display:block;height:100%;background:var(--blue)}.bar-row:nth-child(2n) i{background:#39799d}.sources{font-family:Consolas,monospace;font-size:12px}
    @media(max-width:800px){main{margin:0;padding:24px}.cards{grid-template-columns:1fr 1fr}.bar-row{grid-template-columns:90px 1fr 85px}}
    @media print{body{background:white}main{box-shadow:none;margin:0;max-width:none}h2{break-after:avoid}table,.bar-row{break-inside:avoid}}
    """
    key_cards = [
        ("Read mappate", pct(flagstat.get("mapped_pct"))),
        ("Profondità media", f"{mean_depth:.2f}×" if mean_depth is not None else "n/d"),
        ("Breadth ≥1×", pct(float(depth.get("pct_ge1", 0))) if depth else pct(mean_breadth)),
        ("Varianti PASS", pct(pass_pct)),
    ]
    card_html = "".join(f'<div class="card"><span>{html.escape(k)}</span><b>{html.escape(v)}</b></div>' for k, v in key_cards)
    cov_rows = [
        [
            row["#rname"],
            number(row["numreads"]),
            pct(row["coverage"], 2),
            number(row["meandepth"], 2),
            number(dict(rates).get(row["#rname"], 0), 1),
        ]
        for row in primary
    ]
    variant_rows = [
        ["Record", number(raw["records"]), number(filtered["records"])],
        ["Alleli SNP", number(raw["snp_alleles"]), number(filtered["snp_alleles"])],
        ["Alleli indel/complessi", number(raw["indel_complex_alleles"]), number(filtered["indel_complex_alleles"])],
        ["Siti multiallelici", number(raw["multiallelic"]), number(filtered["multiallelic"])],
        ["Ti/Tv", number(raw["titv"], 3), number(filtered["titv"], 3)],
        ["DP medio dei genotipi", number(raw["mean_dp"], 2), number(filtered["mean_dp"], 2)],
        ["Record PASS", number(raw_pass), number(filtered_pass)],
    ]
    limitations = [
        "Il benchmark GIAB non è stato eseguito: precision, recall e F1 non sono misurati.",
        "Gli hard filter sono controlli tecnici conservativi e non sostituiscono VQSR o una validazione clinica.",
        "BQSR usa il dbSNP disponibile nel progetto; non è un bundle Broad completo di known-sites.",
        f"snpEff ha segnalato {number(ann_errors)} record su contig non riconosciuti." if ann_errors else "snpEff non ha segnalato contig non riconosciuti.",
        "Le varianti annotate non costituiscono interpretazione clinica.",
    ]
    report = f"""<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
    <title>Report {html.escape(prefix)}</title><style>{css}</style></head><body><main>
    <p class="subtitle">HG002 · {html.escape(mode_label)} · NVIDIA Parabricks</p>
    <h1>Report automatico della pipeline germinale</h1>
    <p>Generato il {datetime.now().strftime("%d/%m/%Y alle %H:%M")}. Riferimento: {html.escape(reference_label())}.</p>
    <div class="verdict {'warning' if ann_errors else ''}">ESITO TECNICO: {status}. Benchmark GIAB escluso per scelta progettuale.</div>
    <div class="cards">{card_html}</div>

    <h2>1. Input e preflight</h2>
    {table(["Controllo", "Risultato"], [
        ["Coppie FASTQ totali", number(fastq_validation.get("total_pairs"))],
        ["Coppie FASTQ campionate per ID", number(fastq_validation.get("checked_pairs"))],
        ["Integrità gzip completa", "verificata" if fastq_validation.get("full_gzip_integrity_checked") else "saltata"],
        ["Identificativi R1/R2", "coerenti" if fastq_validation.get("paired_ids_match") else "non verificati"],
        ["Struttura FASTQ", "valida" if fastq_validation.get("fastq_structure_valid") else "non verificata"],
        ["Lunghezza read osservata", f'{fastq_validation.get("read_length_min", "n/d")}–{fastq_validation.get("read_length_max", "n/d")} bp'],
    ])}

    <h2>2. Allineamento, pairing e duplicati</h2>
    {table(["Metrica", "Valore"], [
        ["Read totali", number(flagstat.get("total"))],
        ["Read mappate", f'{number(flagstat.get("mapped"))} ({pct(flagstat.get("mapped_pct"))})'],
        ["Properly paired", f'{number(flagstat.get("properly_paired"))} ({pct(flagstat.get("properly_paired_pct"))})'],
        ["Singleton", f'{number(flagstat.get("singletons"))} ({pct(flagstat.get("singletons_pct"))})'],
        ["Duplicati", f'{number(flagstat.get("duplicates"))} ({pct(flagstat.get("duplicates_pct"))})'],
        ["Percent duplication Picard", duplicate.get("PERCENT_DUPLICATION", "n/d")],
        ["Optical duplicate pairs", duplicate.get("READ_PAIR_OPTICAL_DUPLICATES", "n/d")],
    ])}

    <h2>3. Copertura e uniformità</h2>
    <p>Profondità media pesata sui cromosomi principali: <b>{number(mean_depth, 2)}×</b>.
    Breadth: <b>{pct(mean_breadth)}</b>. CV autosomico delle read/Mb: <b>{pct(auto_cv)}</b>.</p>
    <h3>Soglie di copertura</h3>{depth_chart}
    <h3>Read mappate per megabase</h3>{chrom_chart}
    {table(["Cromosoma", "Read", "Breadth", "Depth", "Read/Mb"], cov_rows)}

    <h2>4. Variant calling e hard filtering</h2>
    {table(["Metrica", "VCF grezzo", "VCF hard-filtered"], variant_rows)}
    <p>Il VCF grezzo è conservato integralmente. Nel VCF hard-filtered i record non conformi sono marcati nel campo FILTER e non eliminati.</p>

    <h2>5. Annotazione snpEff</h2>
    <p>Prodotti VCF annotato, indice, statistiche snpEff e tabella delle varianti HIGH.
    Record con contig non riconosciuto: <b>{number(ann_errors)}</b>.
    Varianti PASS con impatto HIGH: <b>{number(annotation_qc.get("pass_high_records", 0))}</b>.</p>

    <h2>Tempi di esecuzione</h2>{runtime_chart}
    {table(["Step", "Minuti", "Stato"], timing_rows)}

    <h2>Limiti e interpretazione</h2>
    <ul>{''.join(f'<li>{html.escape(item)}</li>' for item in limitations)}</ul>

    <h2>Output e tracciabilità</h2>
    <div class="sources"><p>output/{prefix}.bam (+ .bai)</p>
    <p>output/{prefix}.g.vcf.gz (+ .tbi)</p>
    <p>output/{prefix}.vcf.gz (+ .tbi)</p>
    <p>output/{prefix}.hardfiltered.vcf.gz (+ .tbi)</p>
    <p>output/{prefix}_annotated.vcf.gz (+ .tbi)</p>
    <p>output/{prefix}_high_impact_variants.tsv</p>
    <p>reports/{prefix}_summary.json</p></div>
    </main></body></html>"""
    report_path = reports / f"{prefix}_Report.html"
    report_path.write_text(report, encoding="utf-8")
    print(report_path)


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    validate = sub.add_parser("validate-fastq")
    validate.add_argument("--r1", required=True)
    validate.add_argument("--r2", required=True)
    validate.add_argument("--pairs", type=int, default=100_000)
    validate.add_argument("--full-gzip-checked", action="store_true")
    validate.add_argument("--total-pairs", type=int)
    validate.add_argument("--output", required=True)
    validate.set_defaults(func=validate_fastq)
    report = sub.add_parser("report")
    report.add_argument("--root", required=True)
    report.add_argument("--prefix", required=True)
    report.add_argument("--sample", default="HG002")
    report.add_argument("--mode", choices=("full", "smoke"), default="full")
    report.set_defaults(func=generate_report)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
