process SPLIT_VCF {
    tag "$vcf"

    input:
    path vcf

    output:
    path "*.vcf", emit: individual_vcfs
    path "${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt", emit: var_count
    

    script:
    filt_name = "split_vcf"

    """
    # Get total number of variants (SNPS+indels)
    echo -e "vcf_name\tfilter_name\tvariant_count" >> ${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt
    var_count=\$(bcftools view -H ${vcf} | wc -l)
    echo -e "${vcf.simpleName}\t${filt_name}\t\${var_count}" >> ${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt

    # Get list of samples
    bcftools query -l ${vcf} > samples.txt

    # Split VCF by sample
    while read sample; do
        bcftools view -c1 -s \$sample -Ov ${vcf} > \${sample}.vcf
    done < samples.txt
    """
}