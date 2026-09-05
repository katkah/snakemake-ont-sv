"""
Aggregate QC report with MultiQC.
"""

rule multiqc:
    input:
        expand("results/qc/{sample}/NanoPlot-report.html",     sample=ALL_SAMPLES),
        [f"results/qc/{u.sample_name}/raw/{u.unit_name}/NanoPlot-report.html"
         for u in units_df.itertuples()],
        expand("results/align/{sample}/{sample}_flagstat.txt", sample=ALL_SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    params:
        # Derived from the output rather than hardcoded: a literal prefix
        # breaks on systems without a shared filesystem.
        outdir   = lambda w, output: os.path.dirname(output[0]),
        dirs     = "results/qc results/align"
    log:
        "logs/multiqc.log"
    resources:
        mem_mb = config["mem_mb"]["multiqc"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.dirs} \
                --outdir {params.outdir} \
                --force 2> {log}
        """
