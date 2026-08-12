"""
Alignment: minimap2 map-ont → sorted, indexed BAM.

Units (sequencing runs) are aligned separately so per-run QC stays visible,
then merged into one BAM per sample. Everything downstream reads only the
merged BAM, so a sample with one run and a sample with three look identical
from rule sv_call onwards.
"""

rule align_unit:
    input:
        fastq  = get_unit_fastq,
        genome = config["genome"]
    output:
        bam   = "results/align/{sample}/units/{sample}-{unit}.bam",
        stats = "results/align/{sample}/units/{sample}-{unit}_flagstat.txt"
    log:
        "logs/align/{sample}-{unit}.log"
    threads:
        config["minimap2_threads"]
    resources:
        mem_mb = config["mem_mb"]["align"]
    shadow: "minimal"
    conda:
        "../envs/align.yaml"
    shell:
        # The read group records which run each read came from, so the merged
        # BAM stays traceable back to a unit.
        """
        minimap2 -ax {config[minimap2_preset]} \
                 -t {threads} \
                 -R '@RG\\tID:{wildcards.unit}\\tSM:{wildcards.sample}' \
                 {input.genome} {input.fastq} 2>> {log} \
          | samtools sort -@ {threads} -o {output.bam} 2>> {log}
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        """


rule merge_bams:
    # A single unit is merged too, so the downstream path is the same either way.
    input:
        get_unit_bams
    output:
        bam   = "results/align/{sample}/{sample}_sorted.bam",
        bai   = "results/align/{sample}/{sample}_sorted.bam.bai",
        stats = "results/align/{sample}/{sample}_flagstat.txt"
    log:
        "logs/align/{sample}_merge.log"
    threads: 8
    resources:
        mem_mb = config["mem_mb"]["align"]
    conda:
        "../envs/align.yaml"
    shell:
        """
        samtools merge -@ {threads} -f -o {output.bam} {input} 2> {log}
        samtools index {output.bam} 2>> {log}
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        """
