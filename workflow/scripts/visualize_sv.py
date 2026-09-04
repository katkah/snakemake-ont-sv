#!/usr/bin/env python3
"""
Circular structural variant visualisation using pyCirclize.

Reads a VCF file (typically unique-to-mutant variants from bcftools isec)
and produces a circular genome plot showing:
  - BND (breakends/translocations): inter-chromosomal links in red
  - INV (inversions >= min_size): intra-chromosomal arcs in blue
  - DUP (duplications >= min_size): intra-chromosomal arcs in green

Run by Snakemake through the script: directive — see visualize_sv and
visualize_sv_truvari (visualize.smk) and visualize_sv_joint (sv_sniffles2.smk).
"""

import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from pycirclize import Circos


def choose_tick_interval(max_chrom_size, target_ticks=6):
    """
    Pick a round tick spacing so the largest chromosome gets ~target_ticks marks.

    Keeps the axis readable regardless of organism: a hardcoded interval that
    suits C. elegans (5 Mb) would give 50 overlapping ticks on human chr1 and a
    single tick on a yeast chromosome.
    """
    raw = max(max_chrom_size / target_ticks, 1)
    magnitude = 10 ** math.floor(math.log10(raw))
    for step in (1, 2, 5, 10):
        if step * magnitude >= raw:
            return int(step * magnitude)
    return int(10 * magnitude)


def format_tick(pos, interval):
    """Label a tick in Mb or kb, whichever suits the spacing."""
    if interval >= 1_000_000:
        return f"{pos / 1_000_000:g}Mb"
    if interval >= 1_000:
        return f"{pos / 1_000:g}kb"
    return f"{pos:g}bp"


def parse_vcf(vcf_path, chromosomes, min_inv_size, min_dup_size):
    """
    Parse VCF and extract SVs for plotting.

    Returns:
        bnd_links: list of ((chr1, pos1), (chr2, pos2)) for BND events
        inv_arcs:  list of (chr, start, end) for large inversions
        dup_arcs:  list of (chr, start, end) for large duplications
    """
    valid_chroms = set(chromosomes.keys())
    bnd_links, inv_arcs, dup_arcs = [], [], []

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos   = int(fields[1])
            alt   = fields[4]
            info  = fields[7]

            if chrom not in valid_chroms:
                continue

            # Parse INFO field into dict
            info_dict = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            svtype = info_dict.get("SVTYPE", "")
            try:
                end = int(info_dict.get("END", pos))
                svlen = abs(int(info_dict.get("SVLEN", 0)))
            except ValueError:
                end, svlen = pos, 0

            # BND — inter-chromosomal breakends
            if svtype == "BND" or "<BND>" in alt:
                match = re.search(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]", alt)
                if match:
                    chr2 = match.group(1)
                    pos2 = int(match.group(2))
                    if chr2 in valid_chroms and chr2 != chrom:
                        bnd_links.append(((chrom, pos), (chr2, pos2)))

            # INV — inversions above size threshold
            elif svtype == "INV" and svlen >= min_inv_size:
                inv_arcs.append((chrom, pos, end))

            # DUP — duplications above size threshold
            elif svtype == "DUP" and svlen >= min_dup_size:
                dup_arcs.append((chrom, pos, end))

    return bnd_links, inv_arcs, dup_arcs


def clamp_pos(chrom, pos, chromosomes):
    """Clamp a coordinate into its contig, warning if it started outside."""
    size = chromosomes[chrom]
    if pos > size:
        print(f"WARNING: {chrom}:{pos} lies {pos - size} bp past the end of "
              f"{chrom} ({size} bp); clamping for the plot", file=sys.stderr)
    return min(pos, size)


def make_plot(chromosomes, bnd_links, inv_arcs, dup_arcs, sample, output,
              min_inv_size, min_dup_size, caller="", method=""):
    """Build the circular SV plot with pyCirclize."""

    circos = Circos(chromosomes, space=4)
    # Plots from different callers and different comparison methods look alike,
    # so name both in the figure itself — a PNG often travels without its path.
    title = sample
    subtitle = " · ".join(part for part in (caller, method) if part)
    if subtitle:
        title += f"\n{subtitle}"
    circos.text(title, size=13, r=20)

    # Tick spacing is derived from the largest chromosome so the scale stays
    # readable for any organism, not just C. elegans.
    tick_interval = choose_tick_interval(max(chromosomes.values()))

    # Draw chromosome sectors with labels, plus coordinate ticks on the outer edge.
    # The ring sits at (86, 93) rather than (93, 100) so the outward-facing ticks
    # and their labels fit inside the canvas (pyCirclize max radius = 100).
    for sector in circos.sectors:
        track = sector.add_track((86, 93))
        track.axis(fc="#4C72B0", lw=0)
        track.text(sector.name, color="white", size=10, adjust_rotation=True)

        major_ticks = list(range(0, int(sector.size), tick_interval))
        track.xticks(
            major_ticks,
            labels=[format_tick(t, tick_interval) for t in major_ticks],
            outer=True,
            label_size=7,
            tick_length=2,
            label_orientation="vertical",
        )

    # --- Links ---
    # BND: red inter-chromosomal links
    for (chr1, pos1), (chr2, pos2) in bnd_links:
        # cuteSV can emit a POS one base past the contig end (seen on III at
        # 13,783,802), which pyCirclize rejects. Clamp, then anchor the 50 kb
        # window to the contig end rather than to the clamped start, which
        # would collapse the link to zero width and hide the variant.
        pos1 = clamp_pos(chr1, pos1, chromosomes)
        pos2 = clamp_pos(chr2, pos2, chromosomes)
        end1 = min(pos1 + 50_000, chromosomes[chr1])
        end2 = min(pos2 + 50_000, chromosomes[chr2])
        pos1 = max(0, end1 - 50_000)
        pos2 = max(0, end2 - 50_000)
        circos.link((chr1, pos1, end1),
                    (chr2, pos2, end2),
                    color="red", alpha=0.6, lw=0.5)

    # INV: blue intra-chromosomal arcs
    for chrom, start, end in inv_arcs:
        end = min(end, chromosomes[chrom])
        start = min(clamp_pos(chrom, start, chromosomes), end)
        circos.link((chrom, start, end),
                    (chrom, start, end),
                    color="steelblue", alpha=0.5, lw=0.5)

    # DUP: green intra-chromosomal arcs
    for chrom, start, end in dup_arcs:
        end = min(end, chromosomes[chrom])
        start = min(clamp_pos(chrom, start, chromosomes), end)
        circos.link((chrom, start, end),
                    (chrom, start, end),
                    color="seagreen", alpha=0.5, lw=0.5)

    # --- Render (must happen before accessing circos.ax) ---
    fig = circos.plotfig()

    # --- Legend ---
    legend_elements = []
    if bnd_links:
        legend_elements.append(plt.Line2D([0], [0], color="red",    lw=2, label=f"BND ({len(bnd_links)})"))
    if inv_arcs:
        legend_elements.append(plt.Line2D([0], [0], color="steelblue", lw=2, label=f"INV ≥{min_inv_size//1000}kb ({len(inv_arcs)})"))
    if dup_arcs:
        legend_elements.append(plt.Line2D([0], [0], color="seagreen",  lw=2, label=f"DUP ≥{min_dup_size//1000}kb ({len(dup_arcs)})"))
    if legend_elements:
        circos.ax.legend(handles=legend_elements, loc="lower right",
                         fontsize=9, framealpha=0.8)

    # --- Save ---
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_path))

    # Also save PNG
    png_path = output_path.with_suffix(".png")
    fig.savefig(str(png_path), dpi=150)

    print(f"Saved: {output_path}")
    print(f"Saved: {png_path}")
    print(f"  BND links:  {len(bnd_links)}")
    print(f"  INV arcs:   {len(inv_arcs)}")
    print(f"  DUP arcs:   {len(dup_arcs)}")


# logging
sys.stderr = open(snakemake.log[0], "w")

p = snakemake.params
chromosomes = p.chromosomes

bnd_links, inv_arcs, dup_arcs = parse_vcf(
    snakemake.input.vcf, chromosomes, p.min_inv_size, p.min_dup_size
)

if not any([bnd_links, inv_arcs, dup_arcs]):
    print(f"Warning: no plottable SVs found in {snakemake.input.vcf} "
          f"(after size thresholds INV>={p.min_inv_size}bp, DUP>={p.min_dup_size}bp)")
    # Still produce an empty plot so Snakemake output is satisfied

make_plot(chromosomes, bnd_links, inv_arcs, dup_arcs,
          snakemake.wildcards.sample, snakemake.output.svg,
          p.min_inv_size, p.min_dup_size, p.caller, p.method)
