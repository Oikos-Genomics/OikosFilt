process GET_BI_SNPS {
    tag "$vcf"

    input:
    tuple path(vcf), path(var_count)

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt", emit: var_count

    script:
    filt_name = "bi_snps"
    """
    # Filter for biallelic SNPs
    bcftools view -m2 -M2 -v snps ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf

    # Count variants and save to file
    var_count=\$(bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l)
    echo -e "${vcf.simpleName}\t${filt_name}\t\${var_count}" >> ${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt
    """
}