"""
SV calling: cuteSV
Highest overall F1 score and recall among benchmarked ONT SV callers.
Best choice when sensitivity is the priority.
"""

rule sv_call_cutesv:
    input:
        bam    = "results/align/{sample}/{sample}_sorted.bam",
        bai    = "results/align/{sample}/{sample}_sorted.bam.bai",
        genome = config["genome"]
    output:
        vcf    = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf",
        workdir = directory(f"results/sv_calls/{SV_CALLER}/{{sample}}/cutesv_work")
    log:
        "logs/sv_calls/{sample}_cutesv.log"
    threads:
        config["cutesv_threads"]
    conda:
        "../envs/cutesv.yaml"
    shell:
        """
        mkdir -p {output.workdir}
        cuteSV {input.bam} \
               {input.genome} \
               {output.vcf} \
               {output.workdir} \
               --threads {threads} \
               --min_size {config[min_sv_size]} \
               --min_support {config[min_support_reads]} \
               --genotype 2> {log}
        """
