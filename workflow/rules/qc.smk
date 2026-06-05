"""
Quality control: NanoPlot read quality + tinycov genome coverage.
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
    conda:
        "../envs/qc.yaml"
    shell:
        """
        NanoPlot --bam {input.bam} \
                 --outdir {params.outdir} \
                 --threads {threads} \
                 --no-N50 2> {log}
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
