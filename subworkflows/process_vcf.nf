include { gzipd_vcf } from '../modules/ingress'
include { vcf_to_table } from '../modules/ingress'
include { gatk_table } from '../modules/gatk'
include { clean_vcf } from '../modules/rscripts'

workflow PROCESS_VCF {
    take: 
    gzip
    raw
    fields
    stjude
    cancer

    main:
    gzipd_vcf(gzip)

    // Prepare channel for converting to table with gatk
    vcf_ch = gzipd_vcf.out
        .mix(raw)
        .filter { it[0] != 'qdnaseq' }
        .join(fields)

    gatk_table(vcf_ch)

    // Prepare channel for converting qdnaseq vcf to table with grep
    qdnaseq_ch = gzipd_vcf.out
        .mix(raw)
        .filter { it[0] == 'qdnaseq' }

    vcf_to_table(qdnaseq_ch)

    tables_ch = gatk_table.out
        .mix(vcf_to_table.out)
        .combine(stjude)
        .combine(cancer)


    clean_vcf(tables_ch)
}