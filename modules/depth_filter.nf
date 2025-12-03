include { GET_95ILE_DP } from './filt_95ile_dp'
include { FILT_MIN_DP } from './filt_min_dp'

workflow DEPTH_FILTER {
    take:
    ch_vcf
    min_dp

    main:
    GET_95ILE_DP(ch_vcf)
    FILT_MIN_DP(GET_95ILE_DP.out.filt_vcf, min_dp)

    emit:
    filt_vcf = FILT_MIN_DP.out.filt_vcf
    var_count_95ile = GET_95ILE_DP.out.var_count
    var_count_min_dp = FILT_MIN_DP.out.var_count
}