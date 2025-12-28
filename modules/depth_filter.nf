include { FILT_95ILE_DP } from './filt_95ile_dp'
include { FILT_MIN_DP } from './filt_min_dp'

workflow DEPTH_FILTER {
    take:
    ch_vcf
    ch_var_count
    min_dp

    main:
    // If you do one of these processes, you need to do the other.
    FILT_95ILE_DP(ch_vcf, ch_var_count.collect())
    FILT_MIN_DP(FILT_95ILE_DP.out.filt_vcf, FILT_95ILE_DP.out.var_count.collect(), min_dp)

    emit:
    filt_vcf = FILT_MIN_DP.out.filt_vcf
    var_count = FILT_MIN_DP.out.var_count
}