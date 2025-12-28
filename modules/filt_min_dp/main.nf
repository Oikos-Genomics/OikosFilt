process FILT_MIN_DP {
    tag "$vcf"

    input:
    path vcf
    path var_counts
    val min_dp

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt", emit: var_count

    script:
    filt_name = "min_dp_${min_dp}"
    """
    # Filter for minimum depth
    bcftools filter -i "INFO/DP>=${min_dp}" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf
    # Count variants and save to file

    var_count=\$(bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l)
    echo -e "${vcf.simpleName}\t${filt_name}\t\${var_count}" >> ${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt
    """
}


