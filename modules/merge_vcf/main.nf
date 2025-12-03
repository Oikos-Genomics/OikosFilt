process MERGE_VCFS {
    tag "Merge individual VCFs"

    input:
    path vcfs
    val prefix

    output:
    path "${prefix}.vcf", emit: filt_vcf

    script:
    """
    # Compress all VCF files in parallel
    parallel bgzip ::: *.vcf
    
    # Index all compressed VCF files in parallel
    parallel bcftools index ::: *.vcf.gz
    
    # Create list of VCF files
    ls *.vcf.gz > vcf_list.txt

    # Merge VCF files
    bcftools merge -l vcf_list.txt -Ov -o ${prefix}.vcf

    """
}