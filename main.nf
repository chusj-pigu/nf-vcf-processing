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

    // Create a channel for gene lists:
    stjude_ch = Channel.fromPath(params.stjude_list)
    cancer_ch = Channel.fromPath(params.cancer_genes)

    // Create a channel from all files in the specified directory
    all_vcf_ch = Channel.fromPath("${params.in_dir}/*")

    // Create cnv_ch for files that contain 'cnv' in their name
    cnv_spectre_ch = all_vcf_ch
        .filter { file -> file.getName().contains('cnv_spectre') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'cnv_spectre' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('spectre', file)  // Create a tuple with type 'spectre' and file path
        }
    cnv_qdnaseq_ch = all_vcf_ch
        .filter { file -> file.getName().contains('cnv_qdnaseq') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'cnv_qdnaseq' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('qdnaseq', file)  // Create a tuple with type 'qdnaseq' and file path
        }
    snp_clin_ch = all_vcf_ch
        .filter { file -> file.getName().contains('clinvar') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'clinvar' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('clinvar', file)  // Create a tuple with type 'clinvar' and file path
        }
    snp_ch = all_vcf_ch
        .filter { file -> file.getName().contains('snp.vcf') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'snp' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('snp', file)  // Create a tuple with type 'snp' and file path
        }
    sv_ch = all_vcf_ch
        .filter { file -> file.getName().contains('sv.vcf') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'sv' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('sv', file)  // Create a tuple with type 'sv' and file path
        }
    
    str_ch = all_vcf_ch
        .filter { file -> file.getName().contains('str.vcf') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'sv' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('str', file)  // Create a tuple with type 'sv' and file path
        }
    
    tuples_vcf = cnv_spectre_ch
        .mix(cnv_qdnaseq_ch,snp_clin_ch, snp_ch, sv_ch)
        .branch {                               // Separate .gzip files to be decompressed
            gzip: it[1].name.endsWith('.gz') 
            raw: !it[1].name.endsWith('.gz')
        }

    // Prepare channel for columns to keep in table (gatk step)
    gatk_info = Channel.fromPath(params.info)   // Fields to keep in the INFO category of VCF
        .splitCsv()
        .map { row -> 
            def type = row[0]
            def fields = row[1..-1].findAll { it }.collect { "-F ${it}" }
            tuple(type, fields.join(' '))
        }
    
    gatk_format = Channel.fromPath(params.format)  // Fields to keep in FORMAT category of VCF, they have to be called differently in the function
        .splitCsv()
        .map { row -> 
            def type = row[0]
            def fields = row[1..-1].findAll { it }.collect { "-GF ${it}" }
            tuple(type, fields.join(' '))
        }
    
    fields_ch = gatk_info
        .join(gatk_format, by:0, remainder:true)
        .map { type, info, format ->
            def columns = "${info} ${format}".trim()
            tuple(type, columns)
        }

    PROCESS_VCF(tuples_vcf.gzip,tuples_vcf.raw,fields_ch,stjude_ch,cancer_ch)

}
