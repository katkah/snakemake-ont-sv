"""
Cross-sample comparison and split-read detection.

compare_to_wt: bcftools isec to find SVs unique to each mutant vs WT.
split_reads:   Extract supplementary alignments and inter-chromosomal
               split-read coordinates for visualisation.
sv_stats:      Per-sample SV statistics using analyze_vcf_variants.py.
"""

rule compare_to_wt:
    input:
        mutant_vcf = "results/sv_calls/{sample}/{sample}.vcf.gz",
        mutant_tbi = "results/sv_calls/{sample}/{sample}.vcf.gz.tbi",
        wt_vcf     = get_wt_vcf,
        wt_tbi     = get_wt_vcf_tbi
    output:
        unique_mutant = "results/compare/{sample}_vs_wt/unique_to_{sample}.vcf",
        unique_wt     = "results/compare/{sample}_vs_wt/unique_to_wt.vcf",
        shared        = "results/compare/{sample}_vs_wt/shared.vcf"
    params:
        outdir = "results/compare/{sample}_vs_wt"
    log:
        "logs/compare/{sample}_vs_wt.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools isec {input.mutant_vcf} {input.wt_vcf} \
                      -p {params.outdir} 2> {log}

        # bcftools isec output:
        #   0000.vcf = unique to mutant (first input)
        #   0001.vcf = unique to WT (second input)
        #   0002.vcf = shared (in mutant)
        #   0003.vcf = shared (in WT)
        cp {params.outdir}/0000.vcf {output.unique_mutant}
        cp {params.outdir}/0001.vcf {output.unique_wt}
        cp {params.outdir}/0002.vcf {output.shared}
        """


rule sv_stats:
    input:
        vcf = "results/sv_calls/{sample}/{sample}.vcf.gz"
    output:
        tsv    = "results/sv_stats/{sample}_sv_summary.tsv",
        report = "results/sv_stats/{sample}_sv_report.txt"
    params:
        vcf_dir = "results/sv_calls/{sample}"
    log:
        "logs/sv_stats/{sample}.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/analyze_vcf_variants.py \
               {params.vcf_dir} \
               -p "*.vcf.gz" \
               -r {output.report} \
               -c {output.tsv} 2> {log}
        """


rule split_reads:
    input:
        bam = "results/align/{sample}/{sample}_sorted.bam",
        bai = "results/align/{sample}/{sample}_sorted.bam.bai"
    output:
        split_bam     = "results/split_reads/{sample}/{sample}_split.bam",
        split_bai     = "results/split_reads/{sample}/{sample}_split.bam.bai",
        circos_links  = "results/split_reads/{sample}/{sample}_interchromosomal_links.txt"
    log:
        "logs/split_reads/{sample}.log"
    threads: 8
    conda:
        "../envs/align.yaml"
    shell:
        """
        # Extract supplementary (split) alignments
        samtools view -@ {threads} -f 2048 -b {input.bam} \
          > {output.split_bam} 2>> {log}
        samtools index {output.split_bam} 2>> {log}

        # Generate inter-chromosomal link coordinates for visualisation
        samtools view {output.split_bam} | awk '
        {{
            for (j=12; j<=NF; j++) {{
                if ($j ~ /^SA:Z:/) {{
                    split($j, sa, ":");
                    split(sa[3], coords, ",");
                    chr1   = $3;
                    start1 = $4;
                    chr2   = coords[1];
                    start2 = coords[2];
                    if (chr1 != chr2 && chr1 != "MtDNA" && chr2 != "MtDNA") {{
                        print chr1"\t"start1"\t"(start1+100)"\t"chr2"\t"start2"\t"(start2+100);
                    }}
                }}
            }}
        }}' > {output.circos_links} 2>> {log}
        """
