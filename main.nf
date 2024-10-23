def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/nf-vcf-processing --cnv FILE --sv FILE --gene_list FILE

        Mandatory arguments:
         --sv and/or --cnv              Path to the SV and/or CNV vcf file(s)

         Optional arguments:
         --fusion, --translocation      If --sv is used, at least one of these options has to be used.
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

include { PROCESS_VCF } from './subworkflows/process_vcf'


workflow {

    // vcf channel
    cnv_ch = params.cnv != null ? Channel.fromPath(params.cnv).map { file -> tuple('cnv', file) } : Channel.empty()
    sv_ch = params.sv != null ? Channel.fromPath(params.sv).map { file -> tuple('sv', file) } : Channel.empty()
    snp_ch = params.snp != null ? Channel.fromPath(params.snp).map { file -> tuple('snp', file) } : Channel.empty()
    all_vcf_ch = cnv_ch.mix(sv_ch, snp_ch)

    // Prepare channels for fields to keep in vcf (bcftool step)
    info_ch = Channel.fromPath(params.info)   // Fields to keep in the INFO category of VCF
        .splitCsv()
        .map { row -> 
            def type = row[0]
            def fields = row[1..-1].findAll { it }.collect { "%INFO/${it}" }
            tuple(type, fields.join(' '))
        }

    format_ch = Channel.fromPath(params.format)  // Fields to keep in FORMAT category of VCF, they have to be called differently in the function
        .splitCsv()
        .map { row -> 
            def type = row[0]
            def fields = row[1..-1].findAll { it }
            def formatted_fields = fields ? "[${fields.collect { "%${it}" }.join(':')}]" : ""
            tuple(type, formatted_fields)
        }
        
    fields_ch = info_ch
        .join(format_ch, by:0, remainder:true)
        .map { type, info, format ->
            def columns = format ? "${info} ${format}".trim() : "${info}"
            tuple(type, columns)
        }

    // Prepare channel for columns to keep in table (gatk step)
    gatk_info = info_ch
        .map { type, fields ->
            def new_fields = fields.replaceAll("%INFO/", "-F ")
            tuple(type, new_fields)
        }
    
    gatk_format = format_ch 
        .map { type, cols ->
            def gt_cols = cols.replaceAll("%", "-GF ")
            def new_cols = gt_cols.replaceAll("[\\[\\]]", "")
            def spaced_cols = new_cols.replaceAll(":", " ")
            tuple(type, spaced_cols)
        }
    
    columns_ch = gatk_info
        .join(gatk_format, by:0, remainder:true)
        .map { type, info, format ->
            def columns = format ? "${info} ${format}".trim() : "${info}"
            tuple(type, columns)
        }

    PROCESS_VCF(all_vcf_ch, fields_ch, columns_ch) 

    // Prints message to indicates which tables were processed :

    types_processed = PROCESS_VCF.out.table
        .map { tuple -> tuple[0] }  
        .unique()
        .toList()
    
    types_processed
        .view { types ->
            println "Selected genes were found for event type: ${types.join(', ')}"
        }

}
