process FILT_MIN_DP {
    tag "$vcf"

    input:
    path vcf
    val min_dp

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    tuple val("${vcf.simpleName}"), path("${vcf.simpleName}_${filt_name}_variants.count"), emit: variant_counts

    script:
    filt_name = "min_dp_${min_dp}"
    """
    # Filter for minimum depth
    bcftools filter -i "INFO/DP>=${min_dp}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l > ${vcf.simpleName}_${filt_name}_variants.count

    # Clean up intermediate vcfs
    rm ${vcf}
    """
}