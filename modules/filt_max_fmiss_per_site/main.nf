process FILT_MAX_FMISS_PER_SITE{
    tag "$vcf"

    input:
    path vcf
    val max_fmiss_site

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_variants.count", emit: var_count

    script:
    filt_name = "max_fmiss_per_site_${max_fmiss_site.toString().replace('.','')}"
    """
    # Filter variants based on maximum amount of missing data
    bcftools view -i "F_MISSING < ${max_fmiss_site}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l >> ${params.prefix}_variants.count
    

    """
}