process FILT_MAX_FMISS_PER_IND{
    tag "$vcf"

    input:
    path vcf
    val max_fmiss_ind

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_variants.count", emit: var_count

    script:
    filt_name = "max_fmiss_per_ind_${max_fmiss_ind.toString().replace('.','')}"
    """

    vcftools --vcf ${vcf} --missing-indv
    cat *imiss | awk '{print \$1" "\$5 }' | tail -n +2 > inds_and_depths.tsv

    # this constructs a string that bcftools uses to remove inds
    inds_to_remove='^'
    while read i; do
        ind=\$(echo \$i | awk '{print \$1}')
        f_miss=\$(echo \$i | awk '{print \$2}')
        if (( \$(echo "\${f_miss} > ${max_fmiss_ind}" | bc -l) ))
        then
            inds_to_remove=\${inds_to_remove}\${ind},
        fi
    done < inds_and_depths.tsv

    inds_to_remove=\${inds_to_remove%?}
    bcftools view ${vcf} -s \${inds_to_remove} -Ov -o ${vcf.simpleName}_${filt_name}.vcf

    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l >> ${params.prefix}_variants.count

    """
}