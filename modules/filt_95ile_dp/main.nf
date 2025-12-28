process FILT_95ILE_DP{
    tag "$vcf"

    input:
    path vcf
    path var_counts

    output:
    path "${vcf.simpleName}_${filt_name}.vcf", emit: filt_vcf
    path "${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt", emit: var_count

    script:
    filt_name = "95ile_dp"
    """
    # Calculate 95% percentile for DP
    ind_cov=\$(grep -v '^#' ${vcf} | awk '{print \$10}' FS='\t' OFS='\n' | cut -d: -f3 | grep -oE '[0-9]+')
    lower_bound=\$(echo "\$ind_cov" | tr ' ' '\n' | sort -n | awk 'BEGIN{q=0.025} {a[NR]=\$1} END {print a[int(NR*q)];}')
    upper_bound=\$(echo "\$ind_cov" | tr ' ' '\n' | sort -n | awk 'BEGIN{q=0.975} {a[NR]=\$1} END {print a[int(NR*q)];}')

    # Filter variants based on DP 95% CI
    bcftools filter -e "FORMAT/DP<\$lower_bound || FORMAT/DP>\$upper_bound" ${vcf} -Ov -o ${vcf.simpleName}_${filt_name}.vcf


    ## Find the var_count file with the largest file size (ie, the one with the most entries)
    #var_count_file=\$(du -a ./*_var_count.txt | sort -nr | head -n 1)

    # Count variants and save to file
    var_count=\$(bcftools view -H ${vcf.simpleName}_${filt_name}.vcf | wc -l)
    echo -e "${vcf.simpleName}\t${filt_name}\t\${var_count}" >> ${params.prefix}_${vcf.simpleName}_${filt_name}_var_count.txt
    """
}