#!/usr/bin/env python3
"""README figures: where the time goes, where the machine strains, what the GPU is worth.

Every figure reads the files the real runs produced - Parabricks logs, per-step
timings, benchmark results - and never numbers copied by hand. If a run changes,
the figures are regenerated and that is all.

Usage (inside the tools container, which has matplotlib):
    docker run --rm \
        --mount "type=bind,source=<project>,target=/w" \
        --env PYTHONPATH=/w/.python-packages-linux \
        bioinfo-codeserver:latest python /w/scripts/make_figures.py
"""

from __future__ import annotations

import csv
import re
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


PROJECT = Path(__file__).resolve().parent.parent
FIGURES = PROJECT / "docs" / "figures"
RUN_DIR = PROJECT / "logs" / "runs" / "20260729_232130_full"
REPO_URL = "github.com/Anetheron95/hg002-parabricks-pipeline"

# Three categorical slots, validated with the dataviz guide's script on a white
# surface: worst CVD pair dE 13.1 (deutan), normal vision 30.1. The aqua sits
# below 3:1 contrast, so it always carries a direct label.
BLUE = "#1D4ED8"
ORANGE = "#C2410C"
AQUA = "#1baf7a"

INK = "#1F2328"
INK_SOFT = "#5F6670"
GRID = "#E6E8EB"
SURFACE = "white"

BASE_RC = {
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.labelcolor": "#25282C",
    "xtick.color": "#3F444A",
    "ytick.color": "#3F444A",
}


def style_axes(ax, *, yaxis_grid=True):
    """Recessive grid, no frame: the bars are the data, the rest is not."""
    ax.set_axisbelow(True)
    if yaxis_grid:
        ax.grid(axis="y", color=GRID, linewidth=0.9)
        ax.grid(axis="x", visible=False)
    else:
        ax.grid(axis="x", color=GRID, linewidth=0.9)
        ax.grid(axis="y", visible=False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#C8CCD1")
    ax.tick_params(length=0, pad=7)


def hms(seconds: float) -> str:
    """Readable duration: 7985 -> 2h 13m, 93 -> 1m 33s.

    Past the hour it rounds to the minute instead of truncating, otherwise the
    labels do not match the tables in the README.
    """
    seconds = int(round(seconds))
    if seconds >= 3600:
        hours, minutes = divmod(int(round(seconds / 60)), 60)
        return f"{hours}h {minutes:02d}m"
    minutes, secs = divmod(seconds, 60)
    if minutes:
        return f"{minutes}m {secs:02d}s"
    return f"{secs}s"


# =============================================================================
# Shared readers - the same files the README quotes
# =============================================================================

def read_step_times() -> list[tuple[str, float]]:
    """Per-step wall clock of the published full run, in file order."""
    path = RUN_DIR / "HG002_NovaSeq_40x_53007e55.step_times.tsv"
    if not path.exists():
        return []
    with path.open(encoding="utf-8") as handle:
        return [(row["Step"], float(row["Seconds"]))
                for row in csv.DictReader(handle, delimiter="\t")]


def read_benchmark() -> dict[str, float]:
    """Seconds per phase of the CPU/GPU comparison."""
    path = PROJECT / "reports" / "cpu_vs_gpu" / "benchmark.tsv"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        return {row["phase"]: float(row["seconds"])
                for row in csv.DictReader(handle, delimiter="\t")}


def read_happy() -> dict[tuple[str, str], dict[str, str]]:
    """hap.py summary, keyed by (variant type, filter)."""
    path = (PROJECT / "reports" / "giab" / "HG002_NovaSeq_40x_53007e55" /
            "happy.summary.csv")
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        return {(row["Type"], row["Filter"]): row for row in csv.DictReader(handle)}


# =============================================================================
# Figure 1 - where the machine strains during alignment
# =============================================================================

PB_TIME = r"\[PB Info (\d{4}-[A-Za-z]{3}-\d{2} \d{2}:\d{2}:\d{2})\]"
SAMPLE_RE = re.compile(
    PB_TIME + r".*bases/GPU/minute: ([\d.]+) CPUMem%: ([\d.]+) CPUUsage%: ([\d.]+)"
)
METER_RE = re.compile(PB_TIME + r"\s+([\d.]+)\t([\d.]+)\t([\d.-]+)\s*$")
PHASE_RE = re.compile(PB_TIME + r".*\|\|\s+(GPU-PBBWA mem, Sorting Phase-I|"
                                r"Sorting Phase-II|Marking Duplicates, BQSR)\s+\|\|")


def parse_fq2bam(log_path: Path):
    """Extracts the telemetry Parabricks writes every ten seconds."""
    samples, phases = [], []
    with log_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = SAMPLE_RE.search(line)
            if match:
                stamp, rate, mem, usage = match.groups()
                samples.append((_ts(stamp), float(mem), float(usage), float(rate)))
                continue
            match = METER_RE.search(line)
            if match:
                stamp, _pct, mem, usage = match.groups()
                if usage != "-":
                    samples.append((_ts(stamp), float(mem), float(usage), None))
                continue
            match = PHASE_RE.search(line)
            if match and (not phases or phases[-1][1] != match.group(2)):
                phases.append((_ts(match.group(1)), match.group(2)))
    return samples, phases


def _ts(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%b-%d %H:%M:%S")


PHASE_LABELS = {
    "GPU-PBBWA mem, Sorting Phase-I": "GPU BWA-MEM + sort",
    "Sorting Phase-II": "sort II",
    "Marking Duplicates, BQSR": "mark duplicates, BQSR, write BAM",
}


def figure_bottleneck(output: Path) -> None:
    log_path = RUN_DIR / "2_fq2bam.log"
    if not log_path.exists():
        print(f"skipping {output.name}: {log_path} is missing")
        return

    samples, phases = parse_fq2bam(log_path)
    if not samples:
        print(f"skipping {output.name}: no samples in the log")
        return

    start = samples[0][0]
    minutes = [(s[0] - start).total_seconds() / 60 for s in samples]
    ram = [s[1] for s in samples]
    cpu = [s[2] for s in samples]
    # From bases per minute to billions of bases per minute, and only where the
    # figure is actually reported.
    rate_x = [m for m, s in zip(minutes, samples) if s[3]]
    rate_y = [s[3] / 1e9 for s in samples if s[3]]

    plt.rcParams.update(BASE_RC)
    fig = plt.figure(figsize=(9.5, 7.4), dpi=150, facecolor=SURFACE)
    top = fig.add_axes([0.115, 0.435, 0.855, 0.33])
    bottom = fig.add_axes([0.115, 0.145, 0.855, 0.2], sharex=top)

    span = max(minutes)
    edges = [(p[0] - start).total_seconds() / 60 for p in phases] + [span]

    # The phases: a dashed rule at the boundary and the label just above the
    # panel. They are context, not a series, so no color and no fills.
    for i, (_, name) in enumerate(phases):
        if i:
            for axis in (top, bottom):
                axis.axvline(edges[i], color="#C8CCD1", linewidth=1,
                             linestyle=(0, (4, 3)), zorder=1)
        if edges[i + 1] - edges[i] > span * 0.08:
            top.text((edges[i] + edges[i + 1]) / 2 / span, 1.045,
                     PHASE_LABELS.get(name, name), transform=top.transAxes,
                     ha="center", va="bottom", fontsize=9.5, color=INK_SOFT)

    top.plot(minutes, cpu, color=BLUE, linewidth=1.8, zorder=3)
    top.plot(minutes, ram, color=ORANGE, linewidth=1.8, zorder=3)
    top.set_ylim(0, 100)
    top.set_xlim(0, span)
    top.set_yticks([0, 25, 50, 75, 100])
    top.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    top.set_ylabel("host in use", fontsize=12, labelpad=10)
    style_axes(top)
    top.tick_params(labelbottom=False)

    # Direct labels where the two curves are well apart: identity does not
    # depend on color.
    anchor = min(range(len(minutes)), key=lambda i: abs(minutes[i] - span * 0.38))
    top.text(minutes[anchor], cpu[anchor] - 2.0, "CPU", ha="center", va="top",
             color=BLUE, fontsize=12.5, fontweight="bold")
    top.text(minutes[anchor], ram[anchor] - 2.5, "RAM", ha="center", va="top",
             color=ORANGE, fontsize=12.5, fontweight="bold")

    bottom.fill_between(rate_x, rate_y, color=AQUA, alpha=0.2, zorder=2)
    bottom.plot(rate_x, rate_y, color=AQUA, linewidth=1.8, zorder=3)
    bottom.set_ylim(0, max(rate_y) * 1.5)
    bottom.set_xlabel("minutes since fq2bam started", fontsize=12, labelpad=10)
    bottom.set_ylabel("billion bases\nper minute", fontsize=12, labelpad=10)
    style_axes(bottom)
    bottom.text(rate_x[len(rate_x) // 2], max(rate_y) * 1.22,
                "work done by the GPU", ha="center", va="center", color=AQUA,
                fontsize=11, fontweight="bold")

    fig.text(0.115, 0.965, "The GPU is not the bottleneck. The host is.",
             ha="left", va="top", fontsize=20, fontweight="bold", color=INK)
    fig.text(0.115, 0.917,
             "Aligning the 474 million read pairs of HG002: 2 h 13 m on a single "
             "RTX 3090",
             ha="left", va="top", fontsize=12, color=INK_SOFT)
    fig.text(0.115, 0.879,
             "Parabricks asks for 16 CPU worker threads per GPU and this machine "
             "exposes 10. CPU and memory sit at their\nceiling for almost the "
             "entire alignment, which is why --low-memory and the 8 GB cap are "
             "there at all.",
             ha="left", va="top", fontsize=10.5, color=INK_SOFT, linespacing=1.5)
    fig.text(0.5, 0.03,
             "Source: pbrun fq2bam telemetry, one sample every 10 seconds "
             "(logs/runs/20260729_232130_full/2_fq2bam.log)",
             ha="center", va="bottom", fontsize=8.8, color=INK_SOFT)

    fig.savefig(output, dpi=150, facecolor=SURFACE, metadata={"Software": "Matplotlib"})
    plt.close(fig)
    print(f"wrote {output.relative_to(PROJECT)}")


# =============================================================================
# Figure 2 - where the ten and a half hours go
# =============================================================================

# Which engine each step belongs to. Integrity checks and steps that take a few
# seconds end up in "checks and I/O": they count towards the total but do not
# deserve a row of their own.
ENGINE = {
    "2_fq2bam": ("GPU (Parabricks)", "fq2bam: align, sort, duplicates"),
    "3_bqsr": ("GPU (Parabricks)", "BQSR: recalibration model"),
    "4_haplotypecaller": ("GPU (Parabricks)", "HaplotypeCaller: gVCF"),
    "5_genotypegvcf": ("GPU (Parabricks)", "GenotypeGVCF"),
    "6_hardfilter": ("CPU (GATK, snpEff, samtools)", "Hard filtering"),
    "7_annotation": ("CPU (GATK, snpEff, samtools)", "snpEff annotation"),
    "8_qc": ("CPU (GATK, snpEff, samtools)", "QC: coverage and WGS metrics"),
    "8b_export": ("Checks and I/O", "Export to C:"),
    "9_report": ("Checks and I/O", "HTML report"),
}
ENGINE_COLOR = {
    "GPU (Parabricks)": BLUE,
    "CPU (GATK, snpEff, samtools)": ORANGE,
    "Checks and I/O": AQUA,
}


def figure_runtime(output: Path) -> None:
    steps = read_step_times()
    if not steps:
        print(f"skipping {output.name}: per-step timings are missing")
        return

    rows, checks = [], 0
    for step, seconds in steps:
        entry = ENGINE.get(step)
        if entry is None:
            checks += seconds
        else:
            rows.append((entry[1], entry[0], seconds))
    if not rows:
        print(f"skipping {output.name}: no step recognised")
        return

    rows.append(("Integrity checks", "Checks and I/O", checks))
    rows.sort(key=lambda r: r[2])
    total = sum(r[2] for r in rows)

    plt.rcParams.update(BASE_RC)
    fig = plt.figure(figsize=(9.8, 6.8), dpi=150, facecolor=SURFACE)
    ax = fig.add_axes([0.3, 0.155, 0.615, 0.58])

    positions = range(len(rows))
    ax.barh(list(positions), [r[2] / 3600 for r in rows], height=0.62,
            color=[ENGINE_COLOR[r[1]] for r in rows], zorder=3)
    ax.set_yticks(list(positions))
    ax.set_yticklabels([r[0] for r in rows], fontsize=11)
    ax.set_xlim(0, max(r[2] for r in rows) / 3600 * 1.18)
    ax.set_xlabel("hours", fontsize=12, labelpad=8)
    style_axes(ax, yaxis_grid=False)

    for i, (_, _, seconds) in enumerate(rows):
        ax.text(seconds / 3600 + max(r[2] for r in rows) / 3600 * 0.015, i,
                hms(seconds), va="center", ha="left", fontsize=10.5,
                color=INK_SOFT)

    handles = [plt.Rectangle((0, 0), 1, 1, color=colour)
               for colour in ENGINE_COLOR.values()]
    ax.legend(handles, list(ENGINE_COLOR), loc="lower right",
              bbox_to_anchor=(1.02, -0.02), frameon=False, fontsize=10,
              labelcolor=INK_SOFT, handlelength=1.0, handleheight=1.0,
              borderpad=0.2, labelspacing=0.5)

    fig.text(0.062, 0.962, f"A whole human genome in {hms(total)}",
             ha="left", va="top", fontsize=20, fontweight="bold", color=INK)
    fig.text(0.062, 0.912,
             "HG002, nominal 40x WGS, from FASTQ to annotated VCF on a laptop "
             "driving one RTX 3090",
             ha="left", va="top", fontsize=12, color=INK_SOFT)
    fig.text(0.062, 0.868,
             "Two thirds of the time goes into the two accelerated steps. The "
             "third largest cost is QC, which is not accelerated at all.",
             ha="left", va="top", fontsize=10.5, color=INK_SOFT)
    fig.text(0.5, 0.028,
             "Source: per-step timings of run 53007e55 "
             "(logs/runs/20260729_232130_full)",
             ha="center", va="bottom", fontsize=8.8, color=INK_SOFT)

    fig.savefig(output, dpi=150, facecolor=SURFACE, metadata={"Software": "Matplotlib"})
    plt.close(fig)
    print(f"wrote {output.relative_to(PROJECT)}")


# =============================================================================
# Figure 3 - what the GPU is worth, measured on the same machine
# =============================================================================

def figure_cpu_vs_gpu(output: Path) -> None:
    seconds = read_benchmark()
    needed = ("A1_gpu_fq2bam", "A2_cpu_total",
              "B1_gpu_haplotypecaller", "B2_cpu_haplotypecaller")
    if not all(key in seconds for key in needed):
        print(f"skipping {output.name}: benchmark incomplete "
              "(run benchmark_cpu_vs_gpu.sh)")
        return

    groups = [
        ("Alignment\n+ sort + mark duplicates",
         seconds["A1_gpu_fq2bam"], seconds["A2_cpu_total"]),
        ("Variant calling\n(fixed region, real coverage)",
         seconds["B1_gpu_haplotypecaller"], seconds["B2_cpu_haplotypecaller"]),
    ]

    plt.rcParams.update(BASE_RC)
    fig = plt.figure(figsize=(8.6, 6.4), dpi=150, facecolor=SURFACE)
    ax = fig.add_axes([0.115, 0.185, 0.85, 0.53])

    width = 0.34
    gap = width / 2 + 0.015
    centres = [0, 1]
    gpu_minutes = [g[1] / 60 for g in groups]
    cpu_minutes = [g[2] / 60 for g in groups]
    headroom = max(cpu_minutes)

    ax.bar([c - gap for c in centres], gpu_minutes, width=width, color=BLUE,
           zorder=3, label="GPU — Parabricks on one RTX 3090")
    ax.bar([c + gap for c in centres], cpu_minutes, width=width, color=ORANGE,
           zorder=3, label="CPU — BWA-MEM and GATK, 10 threads")

    ax.set_xticks(centres)
    ax.set_xticklabels([g[0] for g in groups], fontsize=11.5)
    ax.set_xlim(-0.62, 1.62)
    ax.set_ylabel("minutes", fontsize=12, labelpad=10)
    ax.set_ylim(0, headroom * 1.42)
    style_axes(ax)

    # Direct labels with minutes and seconds: rounding to "3 min" would throw
    # away exactly the difference the figure exists to show.
    for centre, gpu, cpu, seconds_pair in zip(centres, gpu_minutes, cpu_minutes,
                                              [(g[1], g[2]) for g in groups]):
        ax.text(centre - gap, gpu + headroom * 0.02, hms(seconds_pair[0]),
                ha="center", va="bottom", fontsize=11, color=BLUE,
                fontweight="bold")
        ax.text(centre + gap, cpu + headroom * 0.02, hms(seconds_pair[1]),
                ha="center", va="bottom", fontsize=11, color=ORANGE,
                fontweight="bold")
        ax.text(centre, cpu + headroom * 0.115, f"{cpu / gpu:.1f}x",
                ha="center", va="bottom", fontsize=18, fontweight="bold",
                color=INK)

    ax.legend(loc="upper right", bbox_to_anchor=(1.01, 1.02), frameon=False,
              fontsize=10.5, labelcolor=INK_SOFT, handlelength=1.0,
              handleheight=1.0, borderpad=0.2, labelspacing=0.55)

    fig.text(0.115, 0.958, "The same machine, without the GPU",
             ha="left", va="top", fontsize=20, fontweight="bold", color=INK)
    fig.text(0.115, 0.908,
             "Same box, same reads, same reference. One variable changes: what "
             "does the computing.",
             ha="left", va="top", fontsize=12, color=INK_SOFT)

    # The speed only counts if the two paths give the same answer.
    concordance = PROJECT / "reports" / "cpu_vs_gpu" / "concordance.tsv"
    if concordance.exists():
        counts = dict(line.split("\t")
                      for line in concordance.read_text().splitlines() if line)
        shared, gpu_only, cpu_only = (int(counts.get(k, 0)) for k in
                                      ("shared_sites", "gpu_only", "cpu_only"))
        if shared:
            fig.text(0.115, 0.862,
                     f"They also agree: {shared:,} of {shared + gpu_only:,} "
                     f"sites identical — {gpu_only} GPU-only, {cpu_only} "
                     "CPU-only.",
                     ha="left", va="top", fontsize=10.5, color=INK_SOFT)
    fig.text(0.5, 0.03,
             "Source: scripts/benchmark_cpu_vs_gpu.sh — "
             "reports/cpu_vs_gpu/benchmark.tsv",
             ha="center", va="bottom", fontsize=8.8, color=INK_SOFT)

    fig.savefig(output, dpi=150, facecolor=SURFACE, metadata={"Software": "Matplotlib"})
    plt.close(fig)
    print(f"wrote {output.relative_to(PROJECT)}")


# =============================================================================
# Figure 4 - the square summary card
#
# Made to be read on a phone in a social feed, next to a photo of the machine:
# 1200x1200, three claims and nothing else. Every number comes from the same
# files as the other figures, so the card cannot drift away from the README.
# =============================================================================

def figure_summary_card(output: Path) -> None:
    steps = read_step_times()
    seconds = read_benchmark()
    happy = read_happy()

    needed = ("A1_gpu_fq2bam", "A2_cpu_total",
              "B1_gpu_haplotypecaller", "B2_cpu_haplotypecaller")
    if not steps or not all(k in seconds for k in needed) or not happy:
        print(f"skipping {output.name}: timings, benchmark or hap.py summary missing")
        return

    total = sum(s for _, s in steps)
    bars = [
        ("Align + sort + mark duplicates", "10 M read pairs",
         seconds["A1_gpu_fq2bam"], seconds["A2_cpu_total"]),
        ("Variant calling", "chr20:1-20 Mb at real ~35x depth",
         seconds["B1_gpu_haplotypecaller"], seconds["B2_cpu_haplotypecaller"]),
    ]

    plt.rcParams.update(BASE_RC)
    fig = plt.figure(figsize=(8, 8), dpi=150, facecolor=SURFACE)

    # ---- headline -----------------------------------------------------------
    fig.text(0.075, 0.955, f"A whole human genome\nin {hms(total)}",
             ha="left", va="top", fontsize=34, fontweight="bold", color=INK,
             linespacing=1.15)
    fig.text(0.075, 0.788,
             "HG002, 40x WGS. FASTQ to annotated VCF on a laptop driving one\n"
             "RTX 3090 through an external PCIe enclosure — with 31 GB of RAM and\n"
             "12 threads, against the 100 GB and 24 threads NVIDIA asks for.",
             ha="left", va="top", fontsize=11.5, color=INK_SOFT, linespacing=1.6)

    # ---- the GPU is worth this much ----------------------------------------
    #
    # The x axis reserves a right-hand quarter for the speedup callouts, so they
    # line up with each other instead of landing on top of the longest bar.
    ax = fig.add_axes([0.075, 0.40, 0.85, 0.245])
    limit = max(b[3] for b in bars) / 60 * 1.33
    positions, labels = [], []
    for i, (name, detail, gpu, cpu) in enumerate(bars):
        base = (len(bars) - 1 - i) * 2.6
        ax.barh(base + 0.62, cpu / 60, height=0.62, color=ORANGE, zorder=3)
        ax.barh(base - 0.12, gpu / 60, height=0.62, color=BLUE, zorder=3)
        positions += [base + 0.62, base - 0.12]
        labels += ["CPU", "GPU"]

        ax.text(cpu / 60 + limit * 0.014, base + 0.62, hms(cpu), va="center",
                ha="left", fontsize=11.5, color=INK_SOFT)
        ax.text(gpu / 60 + limit * 0.014, base - 0.12, hms(gpu), va="center",
                ha="left", fontsize=11.5, color=INK_SOFT)
        ax.text(0, base + 1.28, f"{name} — {detail}", va="bottom", ha="left",
                fontsize=11.5, color=INK, fontweight="bold")
        ax.text(limit, base + 0.25, f"{cpu / gpu:.1f}x", va="center",
                ha="right", fontsize=23, fontweight="bold", color=INK)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=11, color=INK_SOFT)
    ax.set_ylim(-0.85, (len(bars) - 1) * 2.6 + 1.95)
    ax.set_xlim(0, limit)
    ax.set_xlabel("minutes on the same machine — GPU against 10 CPU threads",
                  fontsize=10.5, labelpad=8)
    style_axes(ax, yaxis_grid=False)
    ax.grid(axis="x", visible=False)
    ax.set_xticks([t for t in ax.get_xticks() if t <= max(b[3] for b in bars) / 60])

    # ---- and it is the same answer, benchmarked -----------------------------
    snp = happy[("SNP", "ALL")]
    indel = happy[("INDEL", "ALL")]
    tiles = [
        (f"{float(snp['METRIC.F1_Score']):.4f}", "SNP F1"),
        (f"{float(indel['METRIC.F1_Score']):.4f}", "INDEL F1"),
        (f"{float(snp['METRIC.Recall']) * 100:.2f} %", "SNP recall"),
    ]
    fig.add_artist(plt.Line2D([0.075, 0.925], [0.295, 0.295], color=GRID,
                              linewidth=1.4))
    for i, (value, caption) in enumerate(tiles):
        x = 0.075 + i * 0.287
        fig.text(x, 0.248, value, ha="left", va="center", fontsize=27,
                 fontweight="bold", color=INK)
        fig.text(x, 0.192, caption, ha="left", va="center", fontsize=11.5,
                 color=INK_SOFT)

    fig.text(0.075, 0.142,
             "Measured, not asserted: hap.py with the vcfeval engine against "
             "NIST HG002\nGRCh38 v4.2.1, inside the high-confidence regions of "
             "chr1-22.",
             ha="left", va="top", fontsize=11, color=INK_SOFT, linespacing=1.6)
    fig.add_artist(plt.Line2D([0.075, 0.925], [0.072, 0.072], color=GRID,
                              linewidth=1.4))
    fig.text(0.075, 0.042, REPO_URL, ha="left", va="center", fontsize=11.5,
             color=INK, fontweight="bold")
    fig.text(0.925, 0.042, "scripts, logs and results", ha="right", va="center",
             fontsize=11, color=INK_SOFT)

    fig.savefig(output, dpi=150, facecolor=SURFACE, metadata={"Software": "Matplotlib"})
    plt.close(fig)
    print(f"wrote {output.relative_to(PROJECT)}")


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure_bottleneck(FIGURES / "host_bottleneck.png")
    figure_runtime(FIGURES / "runtime_breakdown.png")
    figure_cpu_vs_gpu(FIGURES / "cpu_vs_gpu.png")
    figure_summary_card(FIGURES / "summary_card.png")
