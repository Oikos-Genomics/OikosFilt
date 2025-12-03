process FILT_MIN_DP {
    tag "$vcf"

    input:
    path vcf
    val min_dp

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_variants.count", emit: var_count

    script:
    filt_name = "min_dp_${min_dp}"
    """
    # Filter for minimum depth
    bcftools filter -i "INFO/DP>=${min_dp}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l >> ${params.prefix}_variants.count
    """
}