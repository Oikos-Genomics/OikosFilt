process FILT_MIN_MAF{
    tag "$vcf"

    input:
    path vcf
    val min_maf

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    tuple val("${vcf.simpleName}"), path("${vcf.simpleName}_${filt_name}_variants.count"), emit: variant_counts

    script:
    filt_name = "min_maf_${min_maf.toString().replace('.','')}"
    """
    # Filter variants based on minor allele frequency
    bcftools view -i "MAF > ${min_maf}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l > ${vcf.simpleName}_${filt_name}_variants.count

    # Clean up intermediate vcfs
    rm ${vcf}
    """
}