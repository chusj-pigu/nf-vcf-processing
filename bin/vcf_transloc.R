#!/usr/bin/env Rscript

# Load necessary library
library(tidyverse)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of command line arguments are provided
if(length(args) < 2){
  print("Please include input and output file paths as arguments!")
  stop("Requires command line arguments: <input_file> <output_file>")
}

# Assign input and output file paths
input_file <- args[1]
output_file <- args[2]
if(length(args) > 2) {
  genes <- args[3]
}

# Separate the ANN annotation field into columns and keep only useful columns
vcf_translocation <- read.delim(input_file) %>% 
  separate(ANN, into = c("Allele", "Annotation", "Annotation_impact", "Gene_name", "Gene_id", "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos-length", "CDS.pos-length", "AA.pos-length", "Distance", "Message"), sep = "\\|") %>%
  select(CHROM:ID, Annotation, Gene_name, Gene_id, `HGVS.c`, `HGVS.p`, Annotation_impact, Allele)

# Replace %3 html characters to ;
vcf_translocation$HGVS.c <- gsub("%3", ";", vcf_translocation$HGVS.c)
vcf_translocation$HGVS.p <- gsub("%3", ";", vcf_translocation$HGVS.p)

# Write the processed data to the output file
write_tsv(vcf_translocation, output_file)

if (exists("genes")) {
  gene_fusion <- readLines(genes)
  gene_regex <- paste(gene_fusion, collapse = "|")
  vcf_translocation <- vcf_translocation %>% filter(str_detect(Gene_name, gene_regex))
  write_tsv(vcf_translocation, output_file)
} else {
  write_tsv(vcf_translocation, output_file)
}

