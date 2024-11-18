include { bcf_selectFields } from '../modules/bcftools'
include { grep_vcfGenes as grep_gene_snp } from '../modules/ingress'
include { grep_vcfGenes as grep_gene_sv } from '../modules/ingress'
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
    genes

    main:
    gzipd_vcf(vcf)
    vcf_ch = gzipd_vcf.out
        .join(fields)

    bcf_selectFields(vcf_ch)

    bcf_ch = bcf_selectFields.out
        .join(cols)

    gatk_table(bcf_ch)

    // Make separate channels for each type to run gene filtering:

    table_ch = gatk_table.out
        .branch {
            cnv: it[0] == 'cnv'
            sv: it[0] == 'sv'
            snp: it[0] == 'snp'
        }

    // Filter genes only in sv and snp tables
    grep_gene_snp(table_ch.snp, genes)
    grep_gene_sv(table_ch.sv, genes)
    

    //cnv_ch = final_table_ch.filter { tuple -> tuple[0] == 'cnv' }
    //sv_ch = final_table_ch.filter { tuple -> tuple[0] == 'sv' }
    //snp_ch = final_table_ch.filter { tuple -> tuple[0] == 'snp' }

    clean_cnvTable(table_ch.cnv)
    clean_svTable(grep_gene_sv.out)
    clean_snpTable(grep_gene_snp.out)

    final_ch = clean_cnvTable.out
        .mix(clean_snpTable.out, clean_svTable.out)
    
    emit:
    table = final_ch
}