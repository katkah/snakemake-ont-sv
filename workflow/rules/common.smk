"""
Input functions shared across rules.

"""

def get_unit_fastq(wildcards):
    """The fastq belonging to one sample+unit."""
    row = units_df[(units_df["sample_name"] == wildcards.sample) &
                   (units_df["unit_name"] == wildcards.unit)]
    if len(row) != 1:
        raise ValueError(
            f"units.tsv: expected exactly 1 row for "
            f"{wildcards.sample}/{wildcards.unit}, found {len(row)}"
        )
    return row["fastq"].tolist()

def get_unit_bams(wildcards):
    """Every per-unit BAM belonging to one sample."""
    units = units_df[units_df["sample_name"] == wildcards.sample]["unit_name"]
    return [f"results/align/{wildcards.sample}/units/{wildcards.sample}-{u}.bam"
            for u in units]

def get_controls(sample):
    """Reference samples belonging to the same group as `sample`."""
    grp = samples_df.loc[sample, "group"]
    controls = samples_df[(samples_df["group"] == grp) &
                          (samples_df[VARIABLE] == REFERENCE)]["sample_name"].tolist()
    if not controls:
        raise ValueError(
            f"samples.tsv: sample '{sample}' is in group '{grp}', which "
            f"contains no sample with {VARIABLE} == '{REFERENCE}'"
        )
    return controls

def group_control_bams(wildcards):
    """Merged BAMs of every reference sample in a group (input to pool_controls)."""
    controls = samples_df[(samples_df["group"] == wildcards.group) &
                          (samples_df[VARIABLE] == REFERENCE)]["sample_name"].tolist()
    return [f"results/align/{s}/{s}_sorted.bam" for s in controls]

def get_ctrl_vcf(wildcards):
    """Baseline VCF for bcftools isec and Truvari: the pooled controls of this
    sample's group.

    Pooling suits these two methods because they only ask whether a variant is
    *present* in the baseline — they never read its genotypes. It also
    multiplies control coverage, which is what makes the baseline callset
    comparable in depth to the case sample.

    The joint path does not use this: it needs per-sample genotype columns.
    """
    pooled = f"{samples_df.loc[wildcards.sample, 'group']}_controls"
    return f"results/sv_calls/{SV_CALLER}/{pooled}/{pooled}.vcf.gz"

def get_ctrl_vcf_tbi(wildcards):
    return get_ctrl_vcf(wildcards) + ".tbi"

def joint_sample_list(wildcards):
    """Case sample first, then its controls — this fixes the GT[] indices."""
    return ",".join([wildcards.sample] + get_controls(wildcards.sample))

def joint_gt_filter(wildcards):
    """Case is non-reference and every control is confidently 0/0.

    A control with a missing genotype (./.) is not "ref", so such sites are
    dropped rather than counted as case-specific.
    """
    n = len(get_controls(wildcards.sample))
    return " && ".join(['GT[0]="alt"'] +
                       [f'GT[{i}]="ref"' for i in range(1, n + 1)])
