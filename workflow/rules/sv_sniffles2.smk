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
        vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf",
        snf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.snf"
    log:
        f"logs/sv_calls/{SV_CALLER}/{{sample}}_sniffles2.log"
    threads:
        config["sniffles2_threads"]
    resources:
        mem_mb = 16384
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
        snfs = expand(f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.snf", sample=ALL_SAMPLES)
    output:
        vcf = f"results/sv_joint/{SV_CALLER}/joint.vcf"
    log:
        f"logs/sv_joint/{SV_CALLER}_joint.log"
    threads:
        config["sniffles2_threads"]
    resources:
        mem_mb = 32768
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
