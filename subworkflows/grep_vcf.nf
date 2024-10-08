include { grep_fusion } from '../modules/ingress'
include { grep_vcfGenes } from '../modules/ingress'
include { grep_vcf } from '../modules/ingress'
include { gzipd_vcf } from '../modules/ingress'

workflow GREP_VCF {
    take: 
    vcf
    pattern

    main:
    gzipd_vcf(vcf)

    // Will only execute on fusion vcf
    grep_fusion_ch = grep_fusion(gzipd_vcf.out)

    // Replace the fusion raw fusion vcf with the vcf containing only gene fusion
    grep_fusion_final = gzipd_vcf.out.filter { type, vcf -> type != 'fusion' }
        .mix(grep_fusion_ch)
    
    grep_genes_ch = grep_vcfGenes(grep_fusion_final)

    grep_ids_ch = grep_vcf(grep_genes_ch, pattern)
    
    // Will not execute on vcf where there is no match on either previous step
    final_vcf = grep_ids_ch.filter { tuple ->
        def ids_file = tuple[1]  // path to the ids file
        ids_file.size() > 0  // Keep only those where the ids file is non-empty
        }
    
    emit:
    raw_vcf = gzipd_vcf.out
    final_vcf = final_vcf
}