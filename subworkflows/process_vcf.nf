include { bcf_selectFields } from '../modules/bcftools'
include { grep_vcfGenes } from '../modules/ingress'
include { gzipd_vcf } from '../modules/ingress'
include { gatk_table } from '../modules/gatk'
include { clean_cnvTable } from '../modules/rscripts'
include { clean_svTable } from '../modules/rscripts'
include { clean_snpTable } from '../modules/rscripts'

workflow PROCESS_VCF {
    take: 
    vcf
    fields
    cols

    main:
    gzipd_vcf(vcf)
    vcf_ch = gzipd_vcf.out
        .join(fields)

    bcf_selectFields(vcf_ch)

    bcf_ch = bcf_selectFields.out
        .join(cols)

    gatk_table(bcf_ch)
    
    table_ch = grep_vcfGenes(gatk_table.out)

    cnv_ch = table_ch.filter { tuple -> tuple[0] == 'cnv' }
    sv_ch = table_ch.filter { tuple -> tuple[0] == 'sv' }
    snp_ch = table_ch.filter { tuple -> tuple[0] == 'snp' }

    clean_cnvTable(cnv_ch)
    clean_svTable(sv_ch)
    clean_snpTable(snp_ch)
    
    emit:
    table = table_ch
}