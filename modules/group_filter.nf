include { FILT_MIN_MAF } from './filt_min_maf'
include { FILT_MAX_FMISS_PER_IND } from './filt_max_fmiss_per_ind'
include { FILT_MAX_FMISS_PER_SITE } from './filt_max_fmiss_per_site'

workflow GROUP_FILTER {
    take:
    ch_vcf
    min_maf
    max_fmiss_ind
    max_fmiss_site

    main:
    def ch_maf
    if (min_maf != null) {
        ch_maf = FILT_MIN_MAF(ch_vcf, min_maf)
    } else {
        ch_maf = ch_vcf
    }

    def ch_fmiss_ind
    if (max_fmiss_ind != null) {
        ch_fmiss_ind = FILT_MAX_FMISS_PER_IND(ch_maf.filt_vcf, max_fmiss_ind)
    } else {
        ch_fmiss_ind = ch_maf
    }

    def ch_fmiss_site
    if (max_fmiss_site != null) {
        ch_fmiss_site = FILT_MAX_FMISS_PER_SITE(ch_fmiss_ind.filt_vcf, max_fmiss_site)
    } else {
        ch_fmiss_site = ch_fmiss_ind
    }

    def ch_filt_vcf = ch_fmiss_site

    emit:
    filt_vcf = ch_filt_vcf.filt_vcf
    var_count_min_maf = ch_maf.var_count
    var_count_max_fmiss_ind = ch_fmiss_ind.var_count
    var_count_max_fmiss_site = ch_fmiss_site.var_count

}
