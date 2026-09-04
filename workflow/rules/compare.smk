"""
Cross-sample comparison and split-read detection.

compare_to_ctrl: bcftools isec to find SVs unique to each case sample.
split_reads:     Extract supplementary alignments and inter-chromosomal
                 split-read coordinates for visualisation.
sv_stats:        Per-sample SV statistics using analyze_vcf_variants.py.
"""

rule compare_to_ctrl:
    input:
        mutant_vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz",
        mutant_tbi = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz.tbi",
        ctrl_vcf   = get_ctrl_vcf,
        ctrl_tbi   = get_ctrl_vcf_tbi
    output:
        unique_mutant = f"results/compare/{SV_CALLER}/{{sample}}_vs_ctrl/unique_to_{{sample}}.vcf",
        unique_ctrl   = f"results/compare/{SV_CALLER}/{{sample}}_vs_ctrl/unique_to_ctrl.vcf",
        shared        = f"results/compare/{SV_CALLER}/{{sample}}_vs_ctrl/shared.vcf"
    params:
        outdir = f"results/compare/{SV_CALLER}/{{sample}}_vs_ctrl"
    log:
        f"logs/compare/{SV_CALLER}/{{sample}}_vs_ctrl.log"
    resources:
        mem_mb = config["mem_mb"]["compare"]
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools isec {input.mutant_vcf} {input.ctrl_vcf} \
                      -p {params.outdir} 2> {log}

        # bcftools isec output:
        #   0000.vcf = unique to mutant (first input)
        #   0001.vcf = unique to controls (second input)
        #   0002.vcf = shared (in mutant)
        #   0003.vcf = shared (in controls)
        cp {params.outdir}/0000.vcf {output.unique_mutant}
        cp {params.outdir}/0001.vcf {output.unique_ctrl}
        cp {params.outdir}/0002.vcf {output.shared}
        """


rule compare_to_ctrl_truvari:
    """
    Third case-vs-control comparison method.

    bcftools isec matches on exact POS/REF/ALT, so ONT breakpoint wobble makes
    shared variants look mutant-unique. Truvari matches within a breakpoint
    distance plus size/sequence similarity, and unlike the Sniffles joint path
    it works for every caller and also matches breakends (--bnddist).

    Truvari creates its own output directory and refuses to run if it already
    exists. Declaring it with directory() stops Snakemake pre-creating it —
    the same pattern as sv_call_nanovar.
    """
    input:
        mutant_vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz",
        mutant_tbi = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz.tbi",
        ctrl_vcf   = get_ctrl_vcf,
        ctrl_tbi   = get_ctrl_vcf_tbi,
        genome     = config["genome"]
    output:
        bench   = directory(f"results/compare_truvari/{SV_CALLER}/{{sample}}_vs_ctrl/truvari_out"),
        fp      = f"results/compare_truvari/{SV_CALLER}/{{sample}}_vs_ctrl/fp.vcf.gz",
        summary = f"results/compare_truvari/{SV_CALLER}/{{sample}}_vs_ctrl/summary.json",
        vcf     = f"results/compare_truvari/{SV_CALLER}/{{sample}}_vs_ctrl/unique_to_{{sample}}_truvari.vcf"
    params:
        refdist = config["truvari"]["refdist"],
        pctseq  = config["truvari"]["pctseq"],
        sizemax = config["truvari"]["sizemax"]
    log:
        f"logs/compare_truvari/{SV_CALLER}/{{sample}}_vs_ctrl.log"
    resources:
        mem_mb = config["mem_mb"]["compare"]
    shadow: "minimal"
    conda:
        "../envs/truvari.yaml"
    shell:
        # -b is the baseline (pooled controls), -c the comparison (mutant), so
        # Truvari's "FP" set — present in comp, absent from base — is the
        # mutant-unique callset. The precision/recall/F1 values in summary.json
        # are not meaningful here, because the controls are not a truth set;
        # only the partition is.
        """
        truvari bench -b {input.ctrl_vcf} \
                      -c {input.mutant_vcf} \
                      -o {output.bench} \
                      --reference {input.genome} \
                      --refdist {params.refdist} \
                      --pctseq {params.pctseq} \
                      --sizemax {params.sizemax} 2> {log}

        cp {output.bench}/fp.vcf.gz    {output.fp}
        cp {output.bench}/summary.json {output.summary}

        # visualize_sv.py reads uncompressed VCF
        zcat {output.fp} > {output.vcf} 2>> {log}
        """


rule sv_stats:
    input:
        vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf.gz"
    output:
        tsv    = f"results/sv_stats/{SV_CALLER}/{{sample}}_sv_summary.tsv",
        report = f"results/sv_stats/{SV_CALLER}/{{sample}}_sv_report.txt"
    log:
        f"logs/sv_stats/{SV_CALLER}/{{sample}}.log"
    resources:
        mem_mb = config["mem_mb"]["sv_stats"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_vcf_variants.py"


rule split_reads:
    input:
        bam = "results/align/{sample}/{sample}_sorted.bam",
        bai = "results/align/{sample}/{sample}_sorted.bam.bai"
    output:
        split_bam    = "results/split_reads/{sample}/{sample}_split.bam",
        split_bai    = "results/split_reads/{sample}/{sample}_split.bam.bai",
        circos_links = "results/split_reads/{sample}/{sample}_interchromosomal_links.txt"
    params:
        exclude = config.get("exclude_chroms", []),
        script  = workflow.basedir + "/scripts/extract_split_links.py"
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
          > {output.split_bam} 2> {log}
        samtools index {output.split_bam} 2>> {log}

        # Generate inter-chromosomal link coordinates for visualisation
        samtools view {output.split_bam} | \
        python {params.script} \
               --exclude {params.exclude} \
               --output {output.circos_links} 2>> {log}
        """
