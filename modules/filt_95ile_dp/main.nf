process GET_95ILE_DP{
    tag "$vcf"

    input:
    path vcf

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_variants.count", emit: var_count 

    script:
    filt_name = "95ile_dp"
    """
    # Calculate 95% percentile for DP
    ind_cov=\$(grep -v '^#' ${vcf} | awk '{print \$10}' FS='\t' OFS='\n' | cut -d: -f3 | grep -oE '[0-9]+')
    lower_bound=\$(echo "\$ind_cov" | tr ' ' '\n' | sort -n | awk 'BEGIN{q=0.025} {a[NR]=\$1} END {print a[int(NR*q)];}')
    upper_bound=\$(echo "\$ind_cov" | tr ' ' '\n' | sort -n | awk 'BEGIN{q=0.975} {a[NR]=\$1} END {print a[int(NR*q)];}')

    # Filter variants based on DP 95% CI
    bcftools filter -e "FORMAT/DP<\$lower_bound || FORMAT/DP>\$upper_bound" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf

    # Count variants and save to file
    bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l >> ${params.prefix}_variants.count

    """
}