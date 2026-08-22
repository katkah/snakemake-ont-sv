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
        # Deleted once merge_bams has consumed it; the flagstat is kept, so
        # per-run mapping stats survive without storing the reads twice.
        bam   = temp("results/align/{sample}/units/{sample}-{unit}.bam"),
        stats = "results/align/{sample}/units/{sample}-{unit}_flagstat.txt"
    params:
        preset = config["minimap2_preset"]
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
        minimap2 -ax {params.preset} \
                 -t {threads} \
                 -R '@RG\\tID:{wildcards.unit}\\tSM:{wildcards.sample}' \
                 {input.genome} {input.fastq} 2>> {log} \
          | samtools sort -@ {threads} -o {output.bam} 2>> {log}
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        """


rule merge_bams:
    # A single unit is merged too, so the downstream path is the same either way.
    # Restricted to real samples: the pooled controls share this output pattern
    # but are built by pool_controls instead.
    wildcard_constraints:
        sample="|".join(ALL_SAMPLES)
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


rule pool_controls:
    """
    Merge a group's reference samples into one pooled control BAM.

    Written to the same path shape as a real sample, so the existing SV-calling
    rules produce results/sv_calls/{caller}/{group}_controls/... with no special
    casing. Kept (not temp) so a locus can be inspected in IGV across all
    controls at once — the per-unit @RG tags survive the merge, so reads stay
    traceable to the animal they came from.

    Only isec and Truvari use it; those methods read presence/absence, never
    genotypes, which is what makes pooling biological replicates valid here.
    Joint calling deliberately does not use it, so no read is counted twice.
    """
    input:
        group_control_bams
    output:
        bam   = "results/align/{group}_controls/{group}_controls_sorted.bam",
        bai   = "results/align/{group}_controls/{group}_controls_sorted.bam.bai",
        stats = "results/align/{group}_controls/{group}_controls_flagstat.txt"
    log:
        "logs/align/{group}_controls_pool.log"
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
