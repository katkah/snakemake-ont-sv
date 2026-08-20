"""
SV calling: NanoVar
Designed for low-coverage ONT data (4x homozygous, 8x heterozygous).
Uses a neural network classifier for SV confidence scoring.
"""

rule sv_call_nanovar:
    input:
        # Wrapper API: input.reads, input.ref, and optional input.bed (passed as -f).
        # bai is never read by the wrapper — declared so the index is built first.
        reads = "results/align/{sample}/{sample}_sorted.bam",
        bai   = "results/align/{sample}/{sample}_sorted.bam.bai",
        ref   = config["genome"],
        bed   = config["gap_file"]
    output:
        # No working directory to declare: the wrapper runs NanoVar in a
        # temporary directory and moves the VCF out, deriving the name from
        # NanoVar's own naming logic rather than globbing for it.
        vcf = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf"
    log:
        f"logs/sv_calls/{SV_CALLER}/{{sample}}_nanovar.log"
    threads:
        config["nanovar_threads"]
    resources:
        mem_mb = config["mem_mb"]["nanovar"]
    shadow: "minimal"
    wrapper:
        "v9.16.0/bio/nanovar"
