"""
Visualisation: circular SV plots using pyCirclize.
Runs on case samples only (requires comparison to the pooled controls).
"""

rule visualize_sv:
    input:
        vcf = f"results/compare/{SV_CALLER}/{{sample}}_vs_ctrl/unique_to_{{sample}}.vcf"
    output:
        svg = f"results/visualize/{SV_CALLER}/{{sample}}_sv.svg",
        png = f"results/visualize/{SV_CALLER}/{{sample}}_sv.png"
    params:
        chromosomes  = config["chromosomes"],
        caller       = SV_CALLER,
        method       = "bcftools isec",
        min_inv_size = config["min_inv_size"],
        min_dup_size = config["min_dup_size"]
    log:
        f"logs/visualize/{SV_CALLER}/{{sample}}.log"
    resources:
        mem_mb = config["mem_mb"]["visualize"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/visualize_sv.py"


rule visualize_sv_truvari:
    """Circular plot of the Truvari-derived mutant-unique SVs."""
    input:
        vcf = f"results/compare_truvari/{SV_CALLER}/{{sample}}_vs_ctrl/unique_to_{{sample}}_truvari.vcf"
    output:
        svg = f"results/visualize_truvari/{SV_CALLER}/{{sample}}_sv_truvari.svg",
        png = f"results/visualize_truvari/{SV_CALLER}/{{sample}}_sv_truvari.png"
    params:
        chromosomes  = config["chromosomes"],
        caller       = SV_CALLER,
        method       = "truvari",
        min_inv_size = config["min_inv_size"],
        min_dup_size = config["min_dup_size"]
    log:
        f"logs/visualize_truvari/{SV_CALLER}/{{sample}}.log"
    resources:
        mem_mb = config["mem_mb"]["visualize"]
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/visualize_sv.py"
