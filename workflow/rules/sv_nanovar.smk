"""
SV calling: NanoVar
Designed for low-coverage ONT data (4x homozygous, 8x heterozygous).
Uses a neural network classifier for SV confidence scoring.
"""

rule sv_call_nanovar:
    input:
        bam     = "results/align/{sample}/{sample}_sorted.bam",
        bai     = "results/align/{sample}/{sample}_sorted.bam.bai",
        genome  = config["genome"],
        gap_file = config["gap_file"]
    output:
        vcf    = f"results/sv_calls/{SV_CALLER}/{{sample}}/{{sample}}.vcf",
        outdir = directory(f"results/sv_calls/{SV_CALLER}/{{sample}}/nanovar_out")
    log:
        f"logs/sv_calls/{SV_CALLER}/{{sample}}_nanovar.log"
    threads:
        config["nanovar_threads"]
    resources:
        mem_mb = config["mem_mb"]["nanovar"]
    shadow: "minimal"
    conda:
        "../envs/nanovar.yaml"
    shell:
        """
        nanovar -t {threads} \
                -f {input.gap_file} \
                {input.bam} \
                {input.genome} \
                {output.outdir} 2> {log}

        # Move NanoVar output VCF to expected path
        # NanoVar names the file itself — find it and fail loudly if not exactly one
        vcf=$(ls {output.outdir}/*.nanovar.pass.vcf)
        n=$(echo "$vcf" | wc -l)
        if [ "$n" -ne 1 ]; then
            echo "Expected 1 NanoVar VCF in {output.outdir}, found $n" >&2
            exit 1
        fi
        cp "$vcf" {output.vcf}
        """
