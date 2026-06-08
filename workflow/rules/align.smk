"""
Alignment: minimap2 map-ont → sorted, indexed BAM.
"""

rule align:
    input:
        fastq = get_fastq,
        genome = config["genome"]
    output:
        bam   = "results/align/{sample}/{sample}_sorted.bam",
        stats = "results/align/{sample}/{sample}_flagstat.txt"
    log:
        "logs/align/{sample}.log"
    threads:
        config["minimap2_threads"]
    resources:
        mem_mb = config["mem_mb"]["align"]
    conda:
        "../envs/align.yaml"
    shell:
        """
        minimap2 -ax {config[minimap2_preset]} \
                 -t {threads} \
                 {input.genome} {input.fastq} 2>> {log} \
          | samtools sort -@ {threads} -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        """
