#!/usr/bin/env python3
"""
Circular structural variant visualisation using pyCirclize.

Reads a VCF file (typically unique-to-mutant variants from bcftools isec)
and produces a circular genome plot showing:
  - BND (breakends/translocations): inter-chromosomal links in red
  - INV (inversions >= min_size): intra-chromosomal arcs in blue
  - DUP (duplications >= min_size): intra-chromosomal arcs in green

Usage:
    python visualize_sv.py \\
        --vcf results/compare/mutant_A_vs_wt/unique_to_mutant_A.vcf \\
        --output results/visualize/mutant_A_sv.svg \\
        --chromosomes I:15072434 II:15279421 III:13783801 IV:17493829 V:20924180 X:17718942 \\
        --sample mutant_A \\
        --min-inv-size 1000 \\
        --min-dup-size 1000
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from pycirclize import Circos


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf",      required=True,  help="Input VCF file")
    parser.add_argument("--output",   required=True,  help="Output SVG path")
    parser.add_argument("--chromosomes", required=True, nargs="+",
                        metavar="CHR:SIZE",
                        help="Chromosome names and sizes, e.g. I:15072434 II:15279421")
    parser.add_argument("--sample",   default="sample", help="Sample name for plot title")
    parser.add_argument("--min-inv-size", type=int, default=1000,
                        help="Minimum inversion size to plot (bp) [default: 1000]")
    parser.add_argument("--min-dup-size", type=int, default=1000,
                        help="Minimum duplication size to plot (bp) [default: 1000]")
    return parser.parse_args()


def parse_chromosomes(chrom_args):
    """Parse 'CHR:SIZE' strings into an ordered dict."""
    chromosomes = {}
    for item in chrom_args:
        try:
            name, size = item.split(":")
            chromosomes[name] = int(size)
        except ValueError:
            sys.exit(f"Error: chromosome argument must be CHR:SIZE, got: {item}")
    return chromosomes


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


def make_plot(chromosomes, bnd_links, inv_arcs, dup_arcs, sample, output):
    """Build the circular SV plot with pyCirclize."""

    circos = Circos(chromosomes, space=4)
    circos.text(f"{sample}\nStructural Variants", size=13, r=20)

    # Draw chromosome sectors with labels
    for sector in circos.sectors:
        track = sector.add_track((93, 100))
        track.axis(fc="#4C72B0", lw=0)
        track.text(sector.name, color="white", size=10, adjust_rotation=True)

    # Draw Mb tick marks
    for sector in circos.sectors:
        tick_track = sector.add_track((90, 93))
        tick_track.axis(fc="none", lw=0)
        major_ticks = range(0, sector.size, 5_000_000)
        for t in major_ticks:
            tick_track.xticks([t], labels=[f"{t // 1_000_000}Mb"], label_size=7,
                               tick_length=1.5, label_orientation="vertical")

    # --- Links ---
    # BND: red inter-chromosomal links
    for (chr1, pos1), (chr2, pos2) in bnd_links:
        end1 = min(pos1 + 50_000, chromosomes[chr1])
        end2 = min(pos2 + 50_000, chromosomes[chr2])
        circos.link((chr1, pos1, end1),
                    (chr2, pos2, end2),
                    color="red", alpha=0.6, lw=0.5)

    # INV: blue intra-chromosomal arcs
    for chrom, start, end in inv_arcs:
        circos.link((chrom, start, end),
                    (chrom, start, end),
                    color="steelblue", alpha=0.5, lw=0.5)

    # DUP: green intra-chromosomal arcs
    for chrom, start, end in dup_arcs:
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
        legend_elements.append(plt.Line2D([0], [0], color="steelblue", lw=2, label=f"INV ≥{args.min_inv_size//1000}kb ({len(inv_arcs)})"))
    if dup_arcs:
        legend_elements.append(plt.Line2D([0], [0], color="seagreen",  lw=2, label=f"DUP ≥{args.min_dup_size//1000}kb ({len(dup_arcs)})"))
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


def main():
    global args
    args = parse_args()

    chromosomes = parse_chromosomes(args.chromosomes)
    bnd_links, inv_arcs, dup_arcs = parse_vcf(
        args.vcf, chromosomes, args.min_inv_size, args.min_dup_size
    )

    if not any([bnd_links, inv_arcs, dup_arcs]):
        print(f"Warning: no plottable SVs found in {args.vcf} "
              f"(after size thresholds INV>={args.min_inv_size}bp, DUP>={args.min_dup_size}bp)")
        # Still produce an empty plot so Snakemake output is satisfied

    make_plot(chromosomes, bnd_links, inv_arcs, dup_arcs, args.sample, args.output)


if __name__ == "__main__":
    main()
