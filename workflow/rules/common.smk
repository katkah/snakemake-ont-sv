"""
Common helper rules: index BAM and VCF files.
"""

rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "{prefix}.bam.bai.log"
    conda:
        "../envs/align.yaml"
    shell:
        "samtools index {input} 2> {log}"


rule bgzip_vcf:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    log:
        "{prefix}.vcf.gz.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output} 2> {log}"


rule index_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    log:
        "{prefix}.vcf.gz.tbi.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index --tbi {input} 2> {log}"
