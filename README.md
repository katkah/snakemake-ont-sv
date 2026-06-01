# snakemake-ont-sv

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

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Miniforge)
- Snakemake ≥ 7.32 (`conda install -c bioconda snakemake`)
- All tool dependencies are installed automatically via `--use-conda`

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/katkah/snakemake-ont-sv.git
cd snakemake-ont-sv

# 2. Create your config files from the templates
cp config/config.yaml.template config/config.yaml
cp config/samples.tsv.example  config/samples.tsv

# 3. Edit config/config.yaml — set your genome path, gap file, and sv_caller
# Edit config/samples.tsv    — add your sample names, conditions, and FASTQ paths

# 4. Run (from the repository root)
snakemake --snakefile workflow/Snakefile --cores 8 --use-conda
```

## Configuration

### `config/config.yaml`

| Key | Description | Example |
|---|---|---|
| `samples` | Path to sample sheet | `config/samples.tsv` |
| `genome` | Path to reference FASTA | `/data/genome/ce11.fa` |
| `gap_file` | BED file of assembly gaps (for visualisation) | `/data/genome/gaps.bed` |
| `chromosomes` | Dict of chromosome names and sizes | see template |
| `sv_caller` | Which caller to use: `sniffles2`, `nanovar`, or `cutesv` | `sniffles2` |
| `min_sv_size` | Minimum SV size to report (bp) | `50` |
| `min_support_reads` | Minimum supporting reads | `3` |
| `min_inv_size` | Minimum inversion size to visualise (bp) | `1000` |
| `min_dup_size` | Minimum duplication size to visualise (bp) | `1000` |
| `sniffles2_joint_call` | Enable multi-sample joint calling (Sniffles2 only) | `false` |

See `config/config.yaml.template` for all available options including thread counts and memory limits.

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
│   └── {sample}_sv_summary.tsv          # SV type counts and sizes
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

## Running on HPC (PBS/SLURM)

```bash
snakemake --snakefile workflow/Snakefile \
          --profile profiles/pbs \
          --use-conda
```

A PBS profile template will be added in a future release.

## Organism-specific notes

The pipeline is organism-agnostic. The `config.yaml.template` ships with chromosome sizes for *C. elegans* (WBcel235/ce11) as an example. Replace the `chromosomes:` block with your organism's assembly.

## CI status

[![Pipeline dry-run CI](https://github.com/katkah/snakemake-ont-sv/actions/workflows/test.yml/badge.svg)](https://github.com/katkah/snakemake-ont-sv/actions/workflows/test.yml)
