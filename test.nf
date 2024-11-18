workflow {

    // Create a channel from all files in the specified directory
    all_vcf_ch = Channel.fromPath("${params.in_dir}/*")

    // Create cnv_ch for files that contain 'cnv' in their name
    cnv_spectre_ch = all_vcf_ch
        .filter { file -> 
            file.getName().contains('cnv_spectre') && (file.getName() =~ /.*\.vcf(\.gz)?$/) 
            }  // Filter files with 'cnv_spectre' and ending with '.vcf' or '.vcf.gz'
        .map { file -> 
            tuple('spectre', file)  // Create a tuple with type 'cnv' and file path
        }
        .view()
}