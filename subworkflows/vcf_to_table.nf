include { bcf_filterIDs } from '../modules/bcftools'
include { gatk_table } from '../modules/gatk'
include { clean_cnvTable } from '../modules/rscripts'
include { clean_fusionTable } from '../modules/rscripts'
include { clean_translocTable } from '../modules/rscripts'

workflow VCF_TO_TABLE {

    take:
    vcf_raw
    vcf_filt
    fields_gatk
    gene_file

    main:
    vcf_ch = vcf_raw
        .join(vcf_filt, by:0)
    bcf_filt_ch = bcf_filterIDs(vcf_ch)
        .join(fields_gatk, by:0)
    table_ch = gatk_table(bcf_filt_ch)

    cnv_ch = table_ch.filter { tuple -> tuple[0] == 'cnv' }
    fusion_ch = table_ch.filter { tuple -> tuple[0] == 'fusion' }
    transloc_ch = table_ch.filter { tuple -> tuple[0] == 'translocation' }

    clean_cnvTable(cnv_ch, gene_file)
    clean_fusionTable(fusion_ch, gene_file)
    clean_translocTable(transloc_ch, gene_file)

    emit:
    table = table_ch

}

