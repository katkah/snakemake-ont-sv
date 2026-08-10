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
        f"logs/sv_calls/{SV_CALLER}/{{sample}}_cutesv.log"
    threads:
        config["cutesv_threads"]
    resources:
        mem_mb = config["mem_mb"]["cutesv"]
    shadow: "minimal"
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
               --max_cluster_bias_INS {config[cutesv_max_cluster_bias_ins]} \
               --diff_ratio_merging_INS {config[cutesv_diff_ratio_merging_ins]} \
               --max_cluster_bias_DEL {config[cutesv_max_cluster_bias_del]} \
               --diff_ratio_merging_DEL {config[cutesv_diff_ratio_merging_del]} \
               --genotype 2> {log}
        """
