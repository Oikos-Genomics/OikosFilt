include { FILT_MIN_QUAL } from './filt_min_qual'
include { FILT_MIN_GQ } from './filt_min_gq'

workflow QUALITY_FILTER {
    take:
    ch_vcf
    min_qual
    min_gq

    main:
    FILT_MIN_QUAL(ch_vcf, min_qual)
    FILT_MIN_GQ(FILT_MIN_QUAL.out.filt_vcf, min_gq)

    emit:
    filt_vcf = FILT_MIN_QUAL.out.filt_vcf
    var_count_min_qual = FILT_MIN_QUAL.out.var_count
    var_count_min_gq = FILT_MIN_GQ.out.var_count
}
