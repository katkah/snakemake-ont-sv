#!/usr/bin/env python3
"""
Extract inter-chromosomal split-read link coordinates from a SAM stream.

Reads SAM format from stdin (supplementary alignments only), parses the SA
tag to find the partner chromosome and position, and writes a 6-column BED-like
file suitable for Circos-style visualisation.

Chromosomes listed in --exclude are filtered out from both ends of each link
(typically used to remove mitochondrial DNA, which produces many split reads
that are not informative for nuclear SV detection).

Usage:
    samtools view split.bam | python extract_split_links.py \\
        --exclude MtDNA \\
        --output interchromosomal_links.txt
"""

import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exclude", nargs="*", default=[],
                        metavar="CHROM",
                        help="Chromosomes to exclude from both ends of each link")
    parser.add_argument("--output", required=True,
                        help="Output links file (tab-separated, 6 columns)")
    return parser.parse_args()


def main():
    args = parse_args()
    exclude = set(args.exclude)

    with open(args.output, "w") as out:
        for line in sys.stdin:
            if line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) < 12:
                continue

            chr1   = fields[2]
            start1 = int(fields[3])

            if chr1 in exclude:
                continue

            # Parse SA tag (supplementary alignments)
            for field in fields[11:]:
                if not field.startswith("SA:Z:"):
                    continue
                for sa in field[5:].rstrip(";").split(";"):
                    if not sa:
                        continue
                    parts = sa.split(",")
                    if len(parts) < 2:
                        continue
                    chr2  = parts[0]
                    try:
                        start2 = int(parts[1])
                    except ValueError:
                        continue

                    if chr1 != chr2 and chr2 not in exclude:
                        out.write(
                            f"{chr1}\t{start1}\t{start1 + 100}\t"
                            f"{chr2}\t{start2}\t{start2 + 100}\n"
                        )


if __name__ == "__main__":
    main()
