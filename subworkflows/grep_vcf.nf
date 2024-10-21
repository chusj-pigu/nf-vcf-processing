include { grep_fusion } from '../modules/ingress'
include { grep_vcfGenes } from '../modules/ingress'
include { grep_vcfIDs } from '../modules/ingress'
include { gzipd_vcf } from '../modules/ingress'

workflow GREP_VCF {
    take: 
    vcf
    pattern

    main:
    gzipd_vcf(vcf)
    
    grep_genes_ch = grep_vcfGenes(gzipd_vcf.out)
        .join(pattern, by: 0)

    grep_ids_ch = grep_vcfIDs(grep_genes_ch)
    
    // Will not execute on vcf where there is no match on either previous step
    final_vcf = grep_ids_ch.filter { tuple ->
        def ids_file = tuple[1]  // path to the ids file
        ids_file.size() > 0  // Keep only those where the ids file is non-empty
        }
    
    emit:
    raw_vcf = gzipd_vcf.out
    final_vcf = final_vcf
}