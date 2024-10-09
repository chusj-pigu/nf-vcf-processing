# nf-vcf-processing

[Nextflow] workflow for generating summary tables for gene fusions, translocations and copy number variants from VCF files obtained from [wf-human-variation].

## Dependencies

Requires either [Docker] or [Apptainer] installed.

## Usage

``` sh
nextflow run chusj-pigu/nf-vcf-processing -r main [--cnv FILE] [--sv FILE] [--gene_list FILE] [--sample_id STR]
```

Overview:

This workflow does not require extensive memory or cpu usage, and can be run locally or on the login node of Compute Canada Narval (with `-profile drac`). Requires at least one of `--cnv FILE` or `--sv FILE`.

Filtering VCF files is done with [bcftools], converting to a table is done with [gatk] and arranging the table is done with [R].

### Parameters

- `--cnv`: Path to the VCF file containing CNVs
- `--sv`: Path to the VCF file containing SVs
- `--gene_list`: Path to a file with list of gene names to filter VCF with, each gene must be separated by a newline [default: None]
- `--sample_id`: Name with which to prefix output files [default: "variants"].
- `-profile`: Use drac if running on Narval [default: standard].
- `--out_dir`: Path to the directory that will store output files [default: variants/]
- `--pattern`: Usually not reqiuired. Path to a csv file containing type of event (cnv, fusion or translocation) as first column and pattern with which to filter ID field of VCF with as second column. Default is stored in asstets of this project.
- `--fields`: Usually not reqiuired. Path to a csv file containing type of event (cnv, fusion or translocation) as first column and fields to keep in VCF as subsequent elements. Default is stored in asstets of this project.

## Outputs

Fusions and translocations summary tsv tables will be generated if `--sv` is used, and copy number variants summary tsv if `--cnv` is used.

[Docker]: https://www.docker.com
[Apptainer]: https://apptainer.org
[Nextflow]: https://www.nextflow.io/docs/latest/index.html
[bcftools]: https://samtools.github.io/bcftools/bcftools.html
[gatk]: https://gatk.broadinstitute.org/hc/en-us
[R]: https://www.r-project.org
[wf-human-variation]: https://github.com/epi2me-labs/wf-human-variation