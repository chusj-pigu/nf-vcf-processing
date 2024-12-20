def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run chusj-pigu/nf-vcf-processing --in_dir PATH/TO/VCF

        Mandatory arguments:
         --in_dir                       Path to the SV and/or CNV vcf file(s)

         Optional arguments:
         --stjude_list                  Path to a file containing a list of gene to filter SNP and SV VCF with, each gene must be separated by a newline [default: St-Jude panel]
         --cancer_gene                  Path to a file containing a list of genes to filter CNV files with when they are more than 100 rows long [default: path to a list of 1000s genes involved in cancers]
         --out_dir                      Output directory [default: variants]
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
    
    all_vcf_sources = [
        cnv_spectre_ch, 
        cnv_qdnaseq_ch, 
        snp_clin_ch, 
        snp_ch, 
        sv_ch, 
        str_ch
    ]

    // Combine all VCF channels
    combined_vcf_ch = Channel.empty().mix(*all_vcf_sources)

    // Validation: Ensure at least one VCF file is provided
    combined_vcf_ch
        .ifEmpty {
            log.error "No valid VCF files were found. Please provide at least one VCF file (e.g., --sv, --cnv, or other VCF sources)."
            exit 1
        }
    
    tuples_vcf = combined_vcf_ch
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
