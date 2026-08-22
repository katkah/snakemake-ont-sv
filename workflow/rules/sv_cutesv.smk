"""
SV calling: cuteSV
Highest overall F1 score and recall among benchmarked ONT SV callers.
Best choice when sensitivity is the priority.
"""

rule sv_call_cutesv:
    input:
        bam    = "results/align/{sample}/{sample}_sorted.bam",
        bai    = "results/align/{sample}/{sample}_sorted.bam.bai",
        genome = config["genome"],
        # cuteSV opens the reference from several worker processes at once when
        # writing output. With no index present they race to build one, and a
        # worker reading it half-written fails with "No corresponding contig in
        # reference with X". Declared, not generated: the pipeline writes only
        # to results/ and logs/, never next to the genome.
        fai    = config["genome"] + ".fai"
    output:
        vcf    = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf",
        workdir = directory(f"results/sv_calls/{SV_CALLER}/{{sample}}/cutesv_work")
    params:
        min_size       = config["min_sv_size"],
        min_support    = config["min_support_reads"],
        bias_ins       = config["cutesv_max_cluster_bias_ins"],
        diff_ratio_ins = config["cutesv_diff_ratio_merging_ins"],
        bias_del       = config["cutesv_max_cluster_bias_del"],
        diff_ratio_del = config["cutesv_diff_ratio_merging_del"]
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
               --min_size {params.min_size} \
               --min_support {params.min_support} \
               --max_cluster_bias_INS {params.bias_ins} \
               --diff_ratio_merging_INS {params.diff_ratio_ins} \
               --max_cluster_bias_DEL {params.bias_del} \
               --diff_ratio_merging_DEL {params.diff_ratio_del} \
               --genotype 2> {log}
        """
