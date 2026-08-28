"""
SV calling: Sniffles2
Two-step workflow:
  1. Per-sample: BAM → VCF + SNF (Sniffles2 native format for joint calling)
  2. Joint call: all SNFs → joint VCF (optional, controlled by config)

Joint-genotyping comparison (optional, Sniffles + joint calling only):
  3. joint VCF → per-mutant "unique vs WT" via genotypes → circular plot.
     Runs in parallel to the bcftools/isec comparison in compare.smk (it does
     not replace it); deriving mutant-unique SVs from the jointly genotyped
     calls reduces false positives from breakpoint wobble and missed WT calls.
"""

rule sv_call_sniffles2:
    input:
        # Input keys are the wrapper's API: it reads input.samples and input.ref.
        # bai is never read by the wrapper — it is declared so Snakemake builds
        # the index before Sniffles opens the BAM.
        samples = "results/align/{sample}/{sample}_sorted.bam",
        bai     = "results/align/{sample}/{sample}_sorted.bam.bai",
        ref     = config["genome"]
    output:
        vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf",
        snf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.snf"
    params:
        extra = (
            f"--sample-id {{sample}} "
            f"--minsvlen {config['min_sv_size']} "
            f"--minsupport {config['min_support_reads']}"
        )
    log:
        f"logs/sv_calls/{SV_CALLER}/{{sample}}_sniffles2.log"
    threads:
        config["sniffles2_threads"]
    resources:
        mem_mb = config["mem_mb"]["sniffles2"]
    shadow: "minimal"
    wrapper:
        "v9.16.0/bio/sniffles"


rule sv_joint_sniffles2:
    input:
        samples = expand(f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.snf", sample=ALL_SAMPLES)
    output:
        vcf = f"results/sv_joint/{SV_CALLER}/joint.vcf"
    params:
        # No wildcards in this rule, so a plain string is enough — the joint
        # VCF takes its sample IDs from the .snf files themselves.
        extra = (
            f"--minsvlen {config['min_sv_size']} "
            f"--minsupport {config['min_support_reads']}"
        )
    log:
        f"logs/sv_joint/{SV_CALLER}_joint.log"
    threads:
        config["sniffles2_threads"]
    resources:
        mem_mb = config["mem_mb"]["sniffles2_joint"]
    wrapper:
        "v9.16.0/bio/sniffles"


# ---------------------------------------------------------------------------
# Joint-genotyping-based comparison and visualisation (optional).
# Parallel to compare.smk's bcftools/isec path — does not replace it.
# ---------------------------------------------------------------------------
rule joint_unique_to_mutant:
    """Mutant-unique SVs from the joint VCF: mutant present AND WT absent."""
    input:
        vcf = f"results/sv_joint/{SV_CALLER}/joint.vcf.gz",
        tbi = f"results/sv_joint/{SV_CALLER}/joint.vcf.gz.tbi",
    output:
        vcf = f"results/compare_joint/{SV_CALLER}/{{sample}}_vs_wt/unique_to_{{sample}}_joint.vcf",
    params:
        samples   = joint_sample_list,
        gt_filter = joint_gt_filter,
    log:
        f"logs/compare_joint/{SV_CALLER}/{{sample}}_vs_wt.log"
    resources:
        mem_mb = config["mem_mb"]["compare"]
    conda:
        "../envs/bcftools.yaml"
    shell:
        # -s selects the case sample followed by every control in its group,
        # which fixes the GT[] indices: GT[0] is the case, GT[1..n] the
        # controls. Keep sites where the case is non-reference and *all*
        # controls are confidently 0/0.
        # The second view keeps only the case column (single-sample VCF,
        # matching the shape visualize_sv.py expects).
        r"""
        ( bcftools view -s {params.samples} \
                        -i '{params.gt_filter}' \
                        {input.vcf} \
          | bcftools view -s {wildcards.sample} - > {output.vcf} \
        ) 2> {log}
        """


rule visualize_sv_joint:
    """Circular plot of the joint-derived mutant-unique SVs."""
    input:
        vcf = f"results/compare_joint/{SV_CALLER}/{{sample}}_vs_wt/unique_to_{{sample}}_joint.vcf",
    output:
        svg = f"results/visualize_joint/{SV_CALLER}/{{sample}}_sv_joint.svg",
        png = f"results/visualize_joint/{SV_CALLER}/{{sample}}_sv_joint.png",
    params:
        chromosomes  = config["chromosomes"],
        caller       = SV_CALLER,
        method       = "joint genotyping",
        min_inv_size = config["min_inv_size"],
        min_dup_size = config["min_dup_size"],
    log:
        f"logs/visualize_joint/{SV_CALLER}/{{sample}}.log"
    resources:
        mem_mb = config["mem_mb"]["visualize"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/visualize_sv.py"
