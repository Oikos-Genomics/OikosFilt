process FILT_MIN_GQ {
    tag "$vcf"

    input:
    path vcf
    val min_gq

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_variants.count", emit: var_count

    script:
    filt_name = "min_gq_${min_gq}"
    """
    # Filter for minimum depth
    bcftools filter -i "MIN(FORMAT/GQ)>${min_gq}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l >> ${params.prefix}_variants.count
    """
}