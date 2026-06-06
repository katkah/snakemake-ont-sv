"""
Visualisation: circular SV plots using pyCirclize.
Runs on mutant samples only (requires comparison to WT).
"""

def chrom_args(wildcards):
    """Convert config chromosomes dict to CLI arguments for visualize_sv.py."""
    return " ".join(f"{k}:{v}" for k, v in config["chromosomes"].items())


rule visualize_sv:
    input:
        vcf = f"results/compare/{SV_CALLER}/{{sample}}_vs_wt/unique_to_{{sample}}.vcf"
    output:
        svg = f"results/visualize/{SV_CALLER}/{{sample}}_sv.svg",
        png = f"results/visualize/{SV_CALLER}/{{sample}}_sv.png"
    params:
        chromosomes  = chrom_args,
        min_inv_size = config["min_inv_size"],
        min_dup_size = config["min_dup_size"]
    log:
        f"logs/visualize/{SV_CALLER}/{{sample}}.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/visualize_sv.py \
               --vcf {input.vcf} \
               --output {output.svg} \
               --chromosomes {params.chromosomes} \
               --sample {wildcards.sample} \
               --min-inv-size {params.min_inv_size} \
               --min-dup-size {params.min_dup_size} 2> {log}
        """
