"""Create the LinkedIn chart for HG002 SNP false negatives by chromosome.

Two runs, one variable changed between them: the reference. Iteration 1 used
GRCh38 with ALT contigs and no companion .alt file; iteration 2 used the no-ALT
analysis set. Numbers come from hap.py, counting BD=FN against BD=TP on the
truth sample of reports/giab/<prefix>/happy.vcf.gz.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


CHROMOSOMES = list(range(1, 23))

# % of true SNPs missed, chr1-22, inside the GIAB high-confidence BED.
BEFORE = [
    1.64, 1.34, 1.39, 2.12, 2.07, 10.39, 2.09, 3.50, 2.48, 1.45, 2.28,
    2.21, 1.81, 2.95, 6.47, 3.29, 8.80, 2.49, 5.26, 1.48, 2.91, 4.01,
]
AFTER = [
    0.91, 0.70, 0.50, 0.64, 0.72, 0.64, 0.70, 0.65, 1.54, 0.78, 0.63,
    0.52, 0.58, 0.42, 1.50, 1.05, 0.73, 0.51, 0.56, 0.56, 0.92, 0.76,
]
GENOME_BEFORE = 3.03
GENOME_AFTER = 0.74

# Two categorical slots, validated with the dataviz palette validator:
# worst adjacent pair dE 29.9 (protan), 36.6 (normal vision), all checks pass.
ORANGE = "#C2410C"   # iteration 1, the defect
BLUE = "#1D4ED8"     # iteration 2, after the fix

INK = "#1F2328"
INK_SOFT = "#5F6670"
GRID = "#E6E8EB"
SURFACE = "white"


def build_chart(output_path: Path) -> None:
    """Render the chart to a 1200 x 1200 pixel PNG."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 13,
            "axes.labelcolor": "#25282C",
            "xtick.color": "#3F444A",
            "ytick.color": "#3F444A",
        }
    )

    fig = plt.figure(figsize=(8, 8), dpi=150, facecolor=SURFACE)
    ax = fig.add_axes([0.105, 0.175, 0.855, 0.585])

    # Thin marks with a visible gap between the two bars of a pair: identity
    # never rests on colour alone.
    width = 0.38
    offset = width / 2 + 0.015
    left = [c - offset for c in CHROMOSOMES]
    right = [c + offset for c in CHROMOSOMES]

    ax.bar(left, BEFORE, width=width, color=ORANGE, zorder=3,
           label="Iteration 1 — GRCh38 with ALT contigs, no .alt file")
    ax.bar(right, AFTER, width=width, color=BLUE, zorder=3,
           label="Iteration 2 — no-ALT analysis set")

    # Il margine destro extra ospita le etichette delle medie: cosi' non hanno
    # nessuna barra dietro e non serve mascherarle.
    ax.set_xlim(0.35, 25.9)
    ax.set_ylim(0, 12.6)
    ax.set_xticks(CHROMOSOMES)
    ax.set_xticklabels([str(c) for c in CHROMOSOMES], fontsize=11)
    ax.set_ylabel("% of true SNPs missed", fontsize=14, labelpad=12)
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.tick_params(axis="x", length=0, pad=8)
    ax.tick_params(axis="y", length=0, labelsize=11, pad=7)

    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, linewidth=0.9)
    ax.grid(axis="x", visible=False)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#C8CCD1")

    # Genome averages: contesto recessivo, non una serie. Le etichette stanno
    # nel margine destro, oltre l'ultima barra.
    for value, colour, label in (
        (GENOME_BEFORE, ORANGE, f"genome avg\n{GENOME_BEFORE}%"),
        (GENOME_AFTER, BLUE, f"genome avg\n{GENOME_AFTER}%"),
    ):
        ax.axhline(value, color=colour, linewidth=1.1, linestyle=(0, (5, 4)),
                   alpha=0.5, zorder=2, xmax=0.885)
        ax.text(23.0, value, label, ha="left", va="center",
                color=colour, fontsize=9.5, linespacing=1.3, zorder=5)

    # Il titolare della storia, etichettato direttamente su entrambe le barre.
    ax.text(6 - offset, BEFORE[5] + 0.3, "10.4%", ha="center", va="bottom",
            color=ORANGE, fontsize=19, fontweight="bold", zorder=6)
    ax.text(6 + offset, 1.45, "0.6%", ha="center", va="bottom",
            color=BLUE, fontsize=11.5, fontweight="bold", zorder=6)

    # Blocco esplicativo sopra chr8-chr14, l'unica banda larga davvero libera:
    # li' nessuna barra supera il 3,5 %.
    ax.text(8.0, 10.1,
            "chr6 carries the MHC region,\n"
            "with seven alternate\n"
            "haplotypes. Without the .alt\n"
            "file, those reads get MAPQ 0\n"
            "and are never called.",
            ha="left", va="top", color=INK_SOFT, fontsize=10,
            linespacing=1.55, zorder=6)

    ax.legend(loc="upper right", bbox_to_anchor=(1.005, 1.045), frameon=False,
              fontsize=10, labelcolor=INK_SOFT, handlelength=1.0,
              handleheight=1.0, borderpad=0.2, labelspacing=0.55)

    fig.text(0.105, 0.925, "Where my pipeline was losing variants",
             ha="left", va="top", fontsize=21, fontweight="bold", color=INK)
    fig.text(0.105, 0.877,
             "True SNPs missed against the GIAB v4.2.1 truth set, by chromosome",
             ha="left", va="top", fontsize=12.5, color=INK_SOFT)
    fig.text(0.105, 0.845,
             "One variable changed between the two runs: the reference genome.",
             ha="left", va="top", fontsize=11, color=INK_SOFT)

    fig.text(0.5, 0.045,
             "HG002 · 40x WGS · NVIDIA Parabricks + GATK on a single RTX 3090 · "
             "hap.py/vcfeval, high-confidence regions",
             ha="center", va="bottom", fontsize=8.8, color=INK_SOFT)

    fig.savefig(output_path, dpi=150, facecolor=SURFACE,
                metadata={"Software": "Matplotlib"})
    plt.close(fig)


if __name__ == "__main__":
    build_chart(Path(__file__).with_name("hg002_snp_false_negatives.png"))
