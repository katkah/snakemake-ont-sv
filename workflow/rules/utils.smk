"""
Generic utility rules: index BAM files, compress and index VCF files.

"""

# merge_bams already produces .bai as a declared output (samtools index runs inline).
# index_bam is kept for other BAM files (e.g. split reads).
# ruleorder tells Snakemake to prefer merge_bams when both could produce the same .bai.
ruleorder: merge_bams > index_bam
ruleorder: pool_controls > index_bam


rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "{prefix}.bam.bai.log"
    wrapper:
        "v9.16.0/bio/samtools/index"


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
