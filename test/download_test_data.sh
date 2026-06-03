#!/usr/bin/env bash
# =============================================================================
# Download and prepare test data for snakemake-ont-sv pipeline
#
# Data source:
#   Adams et al. 2024, PeerJ
#   "Identifying transgene insertions in C. elegans genomes with
#    Oxford Nanopore sequencing"
#   https://peerj.com/articles/18100/
#
# Samples:
#   WT  — BY250 (pPdat-1::GFP), SRA: SRR11790534, BioProject: PRJNA627737
#   MUT — UA44  (alpha-synuclein + GFP), SRA: SRR11808611, BioProject: PRJNA627736
#
# Reference genome:
#   C. elegans WBcel235 (ce11), Ensembl release 112
#   Chromosome names: I, II, III, IV, V, X, MtDNA (no chr prefix)
#
# Usage:
#   bash test/download_test_data.sh
#
# Output directory: test/test_data/
# Approximate disk usage: ~3 GB
# =============================================================================

set -euo pipefail   # exit on error, unset variable, or pipe failure

OUTDIR="test/test_data"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== Downloading test data for snakemake-ont-sv ==="
echo "Output directory: $(pwd)"
echo ""

# -----------------------------------------------------------------------------
# 1. Mutant sample — UA44 (alpha-synuclein transgenic)
#    Full run SRR11808611, ~1.6 GB, ~16x genome coverage
# -----------------------------------------------------------------------------
echo "[1/3] Downloading UA44 mutant (SRR11808611) ..."
wget -q --show-progress \
     ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR118/011/SRR11808611/SRR11808611_1.fastq.gz \
     -O mutant_UA44.fastq.gz
echo "      Done: $(du -sh mutant_UA44.fastq.gz | cut -f1)"

# -----------------------------------------------------------------------------
# 2. WT sample — BY250 (GFP dopaminergic marker)
#    Subsample first 400,000 reads (~12x genome coverage) from SRR11790534
#    Full file is 12 GB — streaming avoids downloading it entirely
# -----------------------------------------------------------------------------
echo ""
echo "[2/3] Streaming BY250 WT (SRR11790534), taking first 400,000 reads ..."
echo "      (wget will report a broken pipe at the end — this is expected)"
wget -qO- \
     ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/034/SRR11790534/SRR11790534_1.fastq.gz \
  | zcat | head -n 1600000 | gzip > wt_BY250.fastq.gz || true
echo "      Done: $(du -sh wt_BY250.fastq.gz | cut -f1)"

# -----------------------------------------------------------------------------
# 3. Reference genome — C. elegans WBcel235, Ensembl release 112
#    ~100 MB compressed, ~300 MB uncompressed
# -----------------------------------------------------------------------------
echo ""
echo "[3/3] Downloading C. elegans WBcel235 reference genome ..."
wget -q --show-progress \
     https://ftp.ensembl.org/pub/release-112/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz \
     -O ce11.fa.gz
gunzip ce11.fa.gz
echo "      Done: $(du -sh ce11.fa | cut -f1)"

# -----------------------------------------------------------------------------
# 4. Gap file — empty placeholder (only used for circular visualisation)
#    Replace with a real BED file of assembly gaps if needed
# -----------------------------------------------------------------------------
touch gaps.bed

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=== Download complete ==="
echo ""
echo "Files:"
ls -lh
echo ""
echo "Next steps:"
echo "  1. Copy the config template:"
echo "       cp config/config.yaml.template config/config.yaml"
echo "  2. Edit config/config.yaml — set:"
echo "       genome:   test/test_data/ce11.fa"
echo "       gap_file: test/test_data/gaps.bed"
echo "       samples:  config/samples.tsv"
echo "  3. Create config/samples.tsv — see config/samples.tsv.example"
echo "  4. Run the pipeline:"
echo "       snakemake --snakefile workflow/Snakefile --cores 8 --use-conda"
