"""
Cross-sample comparison and split-read detection.

compare_to_wt: bcftools isec to find SVs unique to each mutant vs WT.
split_reads:   Extract supplementary alignments and inter-chromosomal
               split-read coordinates for visualisation.
sv_stats:      Per-sample SV statistics using analyze_vcf_variants.py.
"""

rule compare_to_wt:
    input:
        mutant_vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz",
        mutant_tbi = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz.tbi",
        wt_vcf     = get_wt_vcf,
        wt_tbi     = get_wt_vcf_tbi
    output:
        unique_mutant = f"results/compare/{SV_CALLER}/{{sample}}_vs_wt/unique_to_{{sample}}.vcf",
        unique_wt     = f"results/compare/{SV_CALLER}/{{sample}}_vs_wt/unique_to_wt.vcf",
        shared        = f"results/compare/{SV_CALLER}/{{sample}}_vs_wt/shared.vcf"
    params:
        outdir = f"results/compare/{SV_CALLER}/{{sample}}_vs_wt"
    log:
        f"logs/compare/{SV_CALLER}/{{sample}}_vs_wt.log"
    resources:
        mem_mb = config["mem_mb"]["compare"]
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
        vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz"
    output:
        tsv    = f"results/sv_stats/{SV_CALLER}/{{sample}}_sv_summary.tsv",
        report = f"results/sv_stats/{SV_CALLER}/{{sample}}_sv_report.txt"
    params:
        vcf_dir = f"results/sv_calls/{SV_CALLER}/{{sample}}"
    log:
        f"logs/sv_stats/{SV_CALLER}/{{sample}}.log"
    resources:
        mem_mb = config["mem_mb"]["sv_stats"]
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/analyze_vcf_variants.py \
               {params.vcf_dir} \
               -p "*.vcf.gz" \
               -r {output.report} \
               -c {output.tsv} > /dev/null 2> {log}
        """


rule split_reads:
    input:
        bam = "results/align/{sample}/{sample}_sorted.bam",
        bai = "results/align/{sample}/{sample}_sorted.bam.bai"
    output:
        split_bam    = "results/split_reads/{sample}/{sample}_split.bam",
        split_bai    = "results/split_reads/{sample}/{sample}_split.bam.bai",
        circos_links = "results/split_reads/{sample}/{sample}_interchromosomal_links.txt"
    params:
        exclude = config.get("exclude_chroms", [])
    log:
        "logs/split_reads/{sample}.log"
    threads: 8
    resources:
        mem_mb = config["mem_mb"]["split_reads"]
    shadow: "minimal"
    conda:
        "../envs/align.yaml"
    shell:
        """
        # Extract supplementary (split) alignments
        samtools view -@ {threads} -f 2048 -b {input.bam} \
          > {output.split_bam} 2>> {log}
        samtools index {output.split_bam} 2>> {log}

        # Generate inter-chromosomal link coordinates for visualisation
        samtools view {output.split_bam} | \
        python workflow/scripts/extract_split_links.py \
               --exclude {params.exclude} \
               --output {output.circos_links} 2>> {log}
        """
