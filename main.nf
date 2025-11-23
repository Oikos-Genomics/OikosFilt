#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { PRINT_HELP } from './modules/print_help'
include { SPLIT_VCF } from './modules/split_vcf'
include { GET_BI_SNPS } from './modules/filt_bi_snps'
include { DEPTH_FILTER } from './modules/depth_filter'
include { QUALITY_FILTER } from './modules/qual_filter'
include { GROUP_FILTER } from './modules/group_filter'
include { MERGE_VCFS } from './modules/merge_vcf'
include { COMPRESS_VCF } from './modules/compress_vcf'

workflow {
    main:

    if (params.help) {
        PRINT_HELP()
        exit 10
    }

    log.info """\
        O I K O S F I L T - N F   P I P E L I N E
        ===================================
        vcf: ${params.vcf}
        prefix: ${params.prefix}
        threads = ${params.threads}
        bi_snps: ${params.bi_snps}
        min_dp: ${params.min_dp}
        dp_95ile = ${params.dp_95ile}
        min_dp = ${params.min_dp}
        min_qual = ${params.min_qual}
        min_gq = ${params.min_gq}
        min_maf = ${params.min_maf}
        max_fmiss_ind = ${params.max_fmiss_ind}
        max_fmiss_site = ${params.max_fmiss_ind}
    """.stripIndent()
    
    if (!params.vcf) {
        error "Please provide a VCF file with --vcf"
    }

    // Input channel
    ch_vcf = channel.fromPath(params.vcf)

    // Split VCF into individual files
    SPLIT_VCF(ch_vcf)

    // filter bi-allelic SNPs only, if desired
    def ch_bi_snps
    if (params.bi_snps) {
        GET_BI_SNPS(SPLIT_VCF.out.individual_vcfs.flatten())
        ch_bi_snps = GET_BI_SNPS.out.filt_vcf.flatten()
    } else {
        ch_bi_snps = SPLIT_VCF.out.individual_vcfs.flatten()
    }

    // apply depth and quality filters, if desired (by default both are applied)
    def ch_dp
    if (params.dp_95ile || params.min_dp) {
        DEPTH_FILTER(ch_bi_snps, params.min_dp)
        ch_dp = DEPTH_FILTER.out.filt_vcf
    } else {
        ch_dp = ch_bi_snps
    }

    def ch_qual
    if (params.min_qual || params.min_gq) {
        QUALITY_FILTER(ch_dp, params.min_qual, params.min_gq)
        ch_qual = QUALITY_FILTER.out.filt_vcf
    } else {
        ch_qual = ch_dp
    }

    // Merge VCFs and apply group-level filters, if desired
    MERGE_VCFS(ch_qual.collect(), params.prefix)

    def ch_group_filt
    if (params.min_maf || params.max_fmiss_ind || params.max_fmiss_site ) {
        GROUP_FILTER(MERGE_VCFS.out.filt_vcf, params.min_maf, params.max_fmiss_ind, params.max_fmiss_site)
        ch_group_filt = GROUP_FILTER.out.filt_vcf
    } else {
        ch_group_filt = MERGE_VCFS.out.filt_vcf
    }

    //FIXME add an optional linkage-filtering step here with plink

    publish:
    final_vcf = COMPRESS_VCF(ch_group_filt)
}

output {
    final_vcf {}
}
