"""
Quality control: NanoPlot (raw reads + alignment) + tinycov genome coverage.
"""

rule nanoplot:
    input:
        bam = "results/align/{sample}/{sample}_sorted.bam",
        bai = "results/align/{sample}/{sample}_sorted.bam.bai"
    output:
        report = "results/qc/{sample}/NanoPlot-report.html"
    params:
        outdir = "results/qc/{sample}"
    log:
        "logs/qc/{sample}_nanoplot.log"
    threads: 8
    resources:
        mem_mb = config["mem_mb"]["nanoplot"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoPlot --bam {input.bam} \
                 --outdir {params.outdir} \
                 --threads {threads} \
                 --no-N50 2> {log}
        """


rule nanoplot_raw:
    # QC of the raw reads *before* alignment (includes reads that fail to map,
    # which the --bam report above excludes). Written to a raw/ subdir so it
    # does not collide with the alignment-based NanoPlot report.
    input:
        fastq = get_unit_fastq
    output:
        report = "results/qc/{sample}/raw/{unit}/NanoPlot-report.html"
    params:
        outdir = "results/qc/{sample}/raw/{unit}"
    log:
        "logs/qc/{sample}-{unit}_nanoplot_raw.log"
    threads: 8
    resources:
        mem_mb = config["mem_mb"]["nanoplot"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoPlot --fastq {input.fastq} \
                 --outdir {params.outdir} \
                 --threads {threads} 2> {log}
        """


rule coverage:
    input:
        bam = "results/align/{sample}/{sample}_sorted.bam",
        bai = "results/align/{sample}/{sample}_sorted.bam.bai"
    output:
        plot = "results/qc/{sample}/{sample}_coverage.png",
        hist = "results/qc/{sample}/{sample}_coverage_hist.png",
        depth = "results/qc/{sample}/{sample}_depth.txt"
    log:
        "logs/qc/{sample}_coverage.log"
    resources:
        mem_mb = config["mem_mb"]["coverage"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        tinycov covplot {input.bam} \
                -o results/qc/{wildcards.sample}/{wildcards.sample}_coverage \
                2>> {log}
        tinycov covhist {input.bam} \
                -o results/qc/{wildcards.sample}/{wildcards.sample}_coverage_hist \
                2>> {log}
        samtools depth -a {input.bam} > {output.depth} 2>> {log}
        """
