include { bcf_filterIDs } from '../modules/bcftools'
include { bcf_selectFields } from '../modules/bcftools'
include { gatk_table } from '../modules/gatk'
include { clean_cnvTable } from '../modules/rscripts'
include { clean_fusionTable } from '../modules/rscripts'
include { clean_translocTable } from '../modules/rscripts'

workflow VCF_TO_TABLE {

    take:
    vcf_raw
    vcf_filt
    fields_bcf
    fields_gatk

    main:
    bcf_filterIDs(vcf_raw, vcf_filt)
    bcf_selectFields(bcf_filterIDs.out, fields_bcf)
    table_ch = gatk_table(bcf_selectFields.out, fields_gatk)
    clean_cnvTable(table_ch)
    clean_fusionTable(table_ch)
    clean_translocTable(table_ch)

}

