process FILT_MIN_QUAL {
    tag "$vcf"

    input:
    path vcf
    path var_counts
    val min_qual

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt", emit: var_count

    script:
    filt_name = "min_qual_${min_qual}"
    """
    # Filter for minimum depth
    bcftools filter -i "MIN(QUAL)>${min_qual}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    
    # Count variants and save to file
    var_count=\$(bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l)
    echo -e "${vcf.simpleName}\t${filt_name}\t\${var_count}" >> ${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt
    """
}