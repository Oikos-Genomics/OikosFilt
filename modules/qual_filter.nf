include { FILT_MIN_QUAL } from './filt_min_qual'
include { FILT_MIN_GQ } from './filt_min_gq'

workflow QUAL_FILTER {
    take:
    ch_vcf
    ch_var_count
    min_qual
    min_gq

    main:
    FILT_MIN_QUAL(ch_vcf, ch_var_count.collect(), min_qual)
    FILT_MIN_GQ(FILT_MIN_QUAL.out.filt_vcf, FILT_MIN_QUAL.out.var_count.collect(), min_gq)

    emit:
    filt_vcf = FILT_MIN_QUAL.out.filt_vcf
    var_count = FILT_MIN_QUAL.out.var_count
}
