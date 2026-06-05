# snakemake-ont-sv

[![Pipeline dry-run CI](https://github.com/katkah/snakemake-ont-sv/actions/workflows/test.yml/badge.svg)](https://github.com/katkah/snakemake-ont-sv/actions/workflows/test.yml)

A Snakemake pipeline for structural variant (SV) detection from Oxford Nanopore long-read sequencing data. Supports three SV callers and produces per-sample statistics, cross-sample comparisons, and circular genome visualisations.

## Features

- **Three SV callers** — choose one per run: [Sniffles2](https://github.com/fritzsedlazeck/Sniffles), [NanoVar](https://github.com/cytham/nanovar), or [cuteSV](https://github.com/tjiangHIT/cuteSV)
- **Alignment** — minimap2 with ONT preset (`map-ont`)
- **QC** — NanoPlot read-quality reports, tinycov coverage plots, MultiQC summary
- **Cross-sample comparison** — mutant vs. wild-type using bcftools isec
- **SV statistics** — per-sample TSV summaries (type counts, size distributions)
- **Split-read detection** — supplementary alignment extraction, inter-chromosomal link coordinates
- **Circular visualisation** — SVG/PNG Circos-style plots via [pyCirclize](https://github.com/moshi4/pyCirclize) (BND, INV, DUP)
- **Sniffles2 joint calling** — optional population-level joint VCF from per-sample `.snf` files
- **CI dry-run** — GitHub Actions tests the full DAG for all three callers on every push

## Requirements

### Conda (one-time setup)

**[Miniforge](https://github.com/conda-forge/miniforge) is strongly recommended** over Miniconda or Anaconda.
It ships with the fast `libmamba` solver by default, uses `conda-forge` instead of `defaults`, and avoids licence issues.

```bash
# Install Miniforge (Linux x86-64)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

If you already have Miniconda/Anaconda, install the libmamba solver and set flexible channel priority:

```bash
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
conda config --set channel_priority flexible
```

> **Why `channel_priority: flexible`?**
> Bioinformatics packages on bioconda (htslib, samtools, bcftools, sniffles) depend on older
> builds of shared libraries (e.g. libdeflate 1.25). When conda-forge ships a newer version
> (1.26+), `strict` priority excludes the bioconda builds entirely. `flexible` allows bioconda's
> builds to be used without conflicting with conda-forge packages.

### Snakemake

```bash
conda create -n snakemake -c bioconda -c conda-forge snakemake>=9
conda activate snakemake
```

All tool dependencies (minimap2, samtools, Sniffles2, NanoPlot, bcftools, etc.) are installed
**automatically** by Snakemake into isolated conda environments on first run — you do not install
them manually.

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/katkah/snakemake-ont-sv.git
cd snakemake-ont-sv

# 2. Create your config files from the templates
cp config/config.yaml.template config/config.yaml
cp config/samples.tsv.example  config/samples.tsv

# 3. Edit config/config.yaml — set your genome path and sample sheet
# Edit config/samples.tsv    — add your sample names, conditions, and FASTQ paths

# 4. Create the logs directory (gitignored, must exist before running)
mkdir -p logs

# 5. Run (from the repository root)
snakemake --snakefile workflow/Snakefile \
          --cores 8 \
          --software-deployment-method conda
```

## Test data

A download script is provided to fetch public ONT data from SRA (C. elegans, ~3 GB):

```bash
bash test/download_test_data.sh
```

This downloads:
- **Mutant**: SRR11808611 — UA44 alpha-synuclein transgenic strain (~1.6 GB, ~16× coverage)
- **WT**: SRR11790534 — BY250 strain, first 400k reads subsampled (~300 MB, ~12× coverage)
- **Reference**: *C. elegans* WBcel235 genome from Ensembl release 112

After downloading, create `config/samples.tsv`:

```
sample_name	condition	fastq
wt_BY250	wt	test/test_data/wt_BY250.fastq.gz
mutant_UA44	mutant	test/test_data/mutant_UA44.fastq.gz
```

And set in `config/config.yaml`:
```yaml
genome:   test/test_data/ce11.fa
gap_file: test/test_data/gaps.bed
```

## Configuration

### `config/config.yaml`

| Key | Description | Example |
|---|---|---|
| `samples` | Path to sample sheet | `config/samples.tsv` |
| `genome` | Path to reference FASTA | `/data/genome/ce11.fa` |
| `gap_file` | BED file of assembly gaps (for NanoVar) | `/data/genome/gaps.bed` |
| `chromosomes` | Dict of chromosome names and sizes | see template |
| `sv_caller` | Which caller to use: `sniffles2`, `nanovar`, or `cutesv` | `sniffles2` |
| `min_sv_size` | Minimum SV size to report (bp) | `50` |
| `min_support_reads` | Minimum supporting reads | `3` |
| `min_inv_size` | Minimum inversion size to visualise (bp) | `1000` |
| `min_dup_size` | Minimum duplication size to visualise (bp) | `1000` |
| `sniffles2_joint_call` | Enable multi-sample joint calling (Sniffles2 only) | `false` |

See `config/config.yaml.template` for all options including per-tool thread counts.

### `config/samples.tsv`

Tab-separated, three columns:

```
sample_name	condition	fastq
WT_01	wt	/data/fastq/wt_01.fastq.gz
mutant_A	mutant	/data/fastq/mutant_a.fastq.gz
```

- `condition` must be `wt` or `mutant`
- At least one `wt` sample is required (used as reference in comparisons)
- Multiple mutant samples are supported

## Outputs

```
results/
├── qc/
│   ├── {sample}/NanoPlot-report.html     # per-sample read QC
│   └── {sample}/{sample}_coverage.png    # coverage plot
├── align/
│   └── {sample}/{sample}.bam             # sorted, indexed alignments
├── sv_calls/
│   └── {sample}/{sample}.vcf.gz          # per-sample SV calls
├── sv_joint/
│   └── joint.vcf.gz                      # joint VCF (Sniffles2 only)
├── compare/
│   └── {sample}_vs_wt/
│       ├── unique_to_{sample}.vcf        # SVs private to mutant
│       ├── unique_to_wt.vcf              # SVs private to WT
│       └── shared.vcf                    # shared SVs
├── sv_stats/
│   └── {sample}_sv_summary.tsv           # SV type counts and sizes
├── split_reads/
│   └── {sample}/{sample}_split.bam       # supplementary alignments
├── visualize/
│   └── {sample}_sv.svg                   # circular SV visualisation
└── multiqc/
    └── multiqc_report.html               # aggregated QC report
```

## SV caller comparison

| Caller | Strengths | Best for |
|---|---|---|
| **Sniffles2** | Highest precision (94%), joint calling, fast | Default choice; population studies |
| **NanoVar** | Works at low coverage (≥4×), neural network scoring | Low-coverage samples |
| **cuteSV** | Highest F1 score (82.5%), best sensitivity | Discovery; maximise recall |

## Known limitations

**Single WT sample.** The comparison step (`bcftools isec`) always uses the first WT sample
in `samples.tsv` as the reference. If you list multiple WT samples, all of them are aligned
and SV-called, but only the first is used in the mutant-vs-WT comparison. The others are
included in QC and MultiQC but silently excluded from the comparison.

If you have multiple WT samples the recommended workaround for now is to merge their VCFs
before running, or simply designate one representative WT sample and list it first.

## Running on HPC (PBS/SLURM)

```bash
snakemake --snakefile workflow/Snakefile \
          --profile profiles/pbs \
          --software-deployment-method conda
```

A PBS profile template will be added in a future release.

## Organism-specific notes

The pipeline is organism-agnostic. The `config.yaml.template` ships with chromosome sizes for
*C. elegans* (WBcel235/ce11) as an example. Replace the `chromosomes:` block with your
organism's assembly.

