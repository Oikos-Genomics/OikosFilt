process SPLIT_VCF {
    tag "$vcf"

    input:
    path vcf

    output:
    path "*.vcf", emit: individual_vcfs
    path "${params.prefix}_variants.count", emit: var_rec

    script:
    filt_name = "split"
    """
    # Get total number of variants (SNPS+indels)
    var_count=\$(bcftools view -H ${vcf} | wc -l)
    echo -e "${vcf.simpleName}\t${filt_name}\t\$var_count" >> ${params.prefix}_variants.count

    # Get list of samples
    bcftools query -l ${vcf} > samples.txt

    # Split VCF by sample
    while read sample; do
        bcftools view -c1 -s \$sample -Ov ${vcf} > \${sample}.vcf
    done < samples.txt
    """
}