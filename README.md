# nf-vcf-processing

[Nextflow] workflow for generating summary tables for STR, CNV, SV and SNP VCF files obtained from [wf-human-variation]. Uses [gatk] and [R].

## Dependencies

Requires either [Docker] or [Apptainer] installed.

## Usage

``` sh
nextflow run chusj-pigu/nf-vcf-processing -r main [--in_dir DIR] 
```

## Inputs

This workflow takes as the only input the directory within which vcf files to be processed are found. At least one vcf containing type of variant in the name (ex. 'cnv') is required.

## Outputs

This workflow outputs one `*_summary.tsv` per vcf file, and outputs two tsv files for Spectre VCF (one that is more detailed and one more summarized with the suffix `_short.tsv`).

### Parameters

- `--in_dir`: Path to the directory containing VCF files to be processed.
- `--stjude_list`: Path to a file containing a list of gene to filter SNP and SV VCF with, each gene must be separated by a newline [default: St-Jude panel].
- `--cancer_genes`: Path to a file containing a list of genes to filter CNV files with when they are more than 100 rows long [default: path to a list of 1000s genes involved in cancers]
- `-profile`: Use drac if running on Narval [default: standard].
- `--out_dir`: Path to the directory that will store output files [default: same as `--in_dir`]
- `--info`: Usually not required. Path to a csv file containing type of vcf as first column and vcf INFO fields to keep. Default is stored in assets of this project.
- `--format`: Usually not required. Path to a csv file containing type of vcf as first column and vcf FORMAT fields to keep. Default is stored in assets of this project.

[Docker]: https://www.docker.com
[Apptainer]: https://apptainer.org
[Nextflow]: https://www.nextflow.io/docs/latest/index.html
[gatk]: https://gatk.broadinstitute.org/hc/en-us
[R]: https://www.r-project.org
[wf-human-variation]: https://github.com/epi2me-labs/wf-human-variation