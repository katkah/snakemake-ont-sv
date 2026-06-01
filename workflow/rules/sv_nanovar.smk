"""
SV calling: NanoVar
Designed for low-coverage ONT data (4x homozygous, 8x heterozygous).
Uses a neural network classifier for SV confidence scoring.
"""

rule sv_call_nanovar:
    input:
        bam     = "results/align/{sample}/{sample}_sorted.bam",
        bai     = "results/align/{sample}/{sample}_sorted.bam.bai",
        genome  = config["genome"],
        gap_file = config["gap_file"]
    output:
        vcf    = "results/sv_calls/{sample}/{sample}.vcf",
        outdir = directory("results/sv_calls/{sample}/nanovar_out")
    log:
        "logs/sv_calls/{sample}_nanovar.log"
    threads:
        config["nanovar_threads"]
    conda:
        "../envs/nanovar.yaml"
    shell:
        """
        nanovar -t {threads} \
                -f {input.gap_file} \
                {input.bam} \
                {input.genome} \
                {output.outdir} 2> {log}

        # Rename NanoVar output VCF to expected path
        cp {output.outdir}/*.nanovar.pass.vcf {output.vcf}
        """
