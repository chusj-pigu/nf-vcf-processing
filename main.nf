include { VCF_TO_TABLE } from './subworkflows/vcf_to_table.nf'
include { GREP_VCF } from './subworkflows/grep_vcf'


workflow {
    // vcf channel
    cnv_ch = params.cnv ? Channel.of(tuple('cnv', params.cnv)) : Channel.empty()
    fusion_ch = params.fusion ? Channel.of(tuple('fusion', params.sv)) : Channel.empty()
    transloc_ch = params.translocation ? Channel.of(tuple('translocation', params.sv)) : Channel.empty()
    all_vcf_ch = cnv_ch.mix(fusion_ch,transloc_ch)

    // Pattern channel
    pattern_ch = Channel.fromPath(params.patterns)
        .splitCsv()
        .map { row -> tuple(row[0],row[1]) }

    // Field channels
    bcf_fields_ch = Channel.fromPath(params.fields)
        .splitCsv()
        .map { row -> 
            def type = row[0]
            def fields = row[1..-1].findAll { it }.collect { "%INFO/${it}" }
            tuple(type, fields.join(' '))
        }

    gatk_fields_ch = bcf_fields_ch
        .map { type, fields ->
            def new_fields = fields.replaceAll("%INFO/", "-F ")
            tuple(type, new_fields)
        }

    GREP_VCF(all_vcf_ch, pattern_ch) 
    VCF_TO_TABLE(GREP_VCF.out.raw_vcf,GREP_VCF.out.final_vcf,bcf_fields_ch,gatk_fields_ch)

}
