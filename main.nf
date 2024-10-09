def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/nf-vcf-processing --cnv FILE --sv FILE --gene_list FILE

        Mandatory arguments:
         --sv and/or --cnv              Path to the SV and/or CNV vcf file(s)

         Optional arguments:
         --gene_list                    Path to a file with list of gene names to filter VCF with, each gene must be separated by a newline [default: EMPTY]
         --sample_id                    String of name with which to prefix output files [default: "variants"]
         --out_dir                      Output directory [default: variants]
         --patterns                     Path to a csv file containing type of event (cnv, fusion or translocation) as first column and pattern with which to filter ID field of VCF with as second column. Default is stored in asstets of this project.
         --fields                       Path to a csv file containing type of event (cnv, fusion or translocation) as first column and fields to keep in VCF as subsequent elements. Default is stored in asstets of this project.
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

include { VCF_TO_TABLE } from './subworkflows/vcf_to_table.nf'
include { GREP_VCF } from './subworkflows/grep_vcf'


workflow {

    // Automatically turn on fusion and translocation is a --sv file is provided
    def fusion = params.sv != null ? true : false
    def translocation = params.sv != null ? true : false

    // vcf channel
    cnv_ch = params.cnv != null ? Channel.of(tuple('cnv', params.cnv)) : Channel.empty()
    fusion_ch = fusion ? Channel.of(tuple('fusion', params.sv)) : Channel.empty()
    transloc_ch = translocation ? Channel.of(tuple('translocation', params.sv)) : Channel.empty()
    all_vcf_ch = cnv_ch.mix(fusion_ch,transloc_ch)

    // Create a channel that will filter types of mutation that are enabled
   // all_types = all_vcf_ch
     //   .map{ tuple -> tuple[0] }

    // Pattern channel
    pattern_ch = Channel.fromPath(params.patterns)
        .splitCsv()
        .map { row -> tuple(row[0],row[1]) }

    // Field channels
    gatk_fields_ch = Channel.fromPath(params.fields)
        .splitCsv()
        .map { row -> 
            def type = row[0]
            def fields = row[1..-1].findAll { it }.collect { "-F ${it}" }
            tuple(type, fields.join(' '))
        }
    
    genes_ch = Channel.fromPath(params.gene_list)

    GREP_VCF(all_vcf_ch, pattern_ch) 
    VCF_TO_TABLE(GREP_VCF.out.raw_vcf,GREP_VCF.out.final_vcf,gatk_fields_ch,genes_ch)

    // Prints message to indicates which tables were processed :

    types_processed = VCF_TO_TABLE.out.table
        .map { tuple -> tuple[0] }   // Extract type field (assuming type is the first element in the tuple)
        .unique()
        .toList()
    
    if (!params.gene_list.endsWith('NO_FILE')) {
    types_processed
        .view { types ->
            println "Selected genes were found for event type: ${types.join(', ')}"
        }
    } else {
    println "No gene list provided, summary tables include all fusion, translocation and cnv events."
    }

}
