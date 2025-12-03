process GET_BI_SNPS {
    tag "$vcf"

    input:
    path vcf

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_variants.count", emit: var_count

    script:
    filt_name = "bi_snps"
    """
    # Filter for biallelic SNPs
    bcftools view -m2 -M2 -v snps ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l >> ${params.prefix}_variants.count

    """
}