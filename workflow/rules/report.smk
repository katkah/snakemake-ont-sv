"""
Aggregate QC report with MultiQC.
"""

rule multiqc:
    input:
        expand("results/qc/{sample}/NanoPlot-report.html",    sample=ALL_SAMPLES),
        expand("results/align/{sample}/{sample}_flagstat.txt", sample=ALL_SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    params:
        outdir   = "results/multiqc",
        dirs     = "results/qc results/align"
    log:
        "logs/multiqc.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.dirs} \
                --outdir {params.outdir} \
                --force 2> {log}
        """
