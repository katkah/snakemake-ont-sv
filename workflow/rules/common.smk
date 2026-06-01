"""
Common helper rules: index BAM and VCF files.
"""

rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    conda:
        "../envs/align.yaml"
    shell:
        "samtools index {input}"


rule bgzip_vcf:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"


rule index_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index --tbi {input}"
