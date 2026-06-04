"""
SV calling: Sniffles2
Two-step workflow:
  1. Per-sample: BAM → VCF + SNF (Sniffles2 native format for joint calling)
  2. Joint call: all SNFs → joint VCF (optional, controlled by config)
"""

rule sv_call_sniffles2:
    input:
        bam = "results/align/{sample}/{sample}_sorted.bam",
        bai = "results/align/{sample}/{sample}_sorted.bam.bai",
        genome = config["genome"]
    output:
        vcf = "results/sv_calls/{sample}/{sample}.vcf",
        snf = "results/sv_calls/{sample}/{sample}.snf"
    log:
        "logs/sv_calls/{sample}_sniffles2.log"
    threads:
        config["sniffles2_threads"]
    conda:
        "../envs/sniffles2.yaml"
    shell:
        """
        sniffles --input {input.bam} \
                 --vcf {output.vcf} \
                 --snf {output.snf} \
                 --reference {input.genome} \
                 --threads {threads} \
                 --minsvlen {config[min_sv_size]} \
                 --minsupport {config[min_support_reads]} 2> {log}
        """


rule sv_joint_sniffles2:
    input:
        snfs = expand("results/sv_calls/{sample}/{sample}.snf", sample=ALL_SAMPLES)
    output:
        vcf = "results/sv_joint/joint.vcf"
    log:
        "logs/sv_joint/sniffles2_joint.log"
    threads:
        config["sniffles2_threads"]
    conda:
        "../envs/sniffles2.yaml"
    shell:
        """
        sniffles --input {input.snfs} \
                 --vcf {output.vcf} \
                 --threads {threads} \
                 --minsvlen {config[min_sv_size]} \
                 --minsupport {config[min_support_reads]} 2> {log}
        """
